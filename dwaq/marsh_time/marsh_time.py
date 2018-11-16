"""
Integrate time water spends in the Calhoun Cut marshes.

Goal is to evaluate spatial influence of restored connections to
marsh.

For usual "age", water from the "source" gets unit concentration,
and zero time.
 time is then incremented by the source concentration * dt

What I want to know is where does water go that spends time in
the marsh.  so I could just do the time integration part.

"""
import six
import os
import shutil
import datetime
from stompy.spatial import wkb2shp

from stompy.model.delft import waq_scenario as waq
from stompy.model.delft import dfm_grid

##

Sub=waq.Substance

class Scen(waq.Scenario):
    # csc_sed_9500: map output interval too small
    # csc_sed_9501: output every 30 minutes
    # 001: atmospheric deposition
    # 002: 0th order nitrification
    # 003: run length of hydro, just 10 days for this run.
    base_path="csc_marsh_time_003"
    map_formats=['binary']

    # With integration_option="15.60"
    # 3000 was unstable.
    # 1000 was unstable.
    #  500 stable at least to 5%, running ~ 600x realtime
    # Also tried 1000 with integration_option="21.70" and "16.70", no better.
    # 21.70 appears to be about 15% faster than 15.60, a little over 700x.
    #  that was outputting a lot.
    integration_option="21.70"
    time_step=300
    map_time_step=3000

    def init_parameters(self):
        params=super(Scen,self).init_parameters()

        params['MaxIter']= 100
        params['Tolerance']= 1.0e-7
        params['Iteration Report']= 0
        params['ONLY_ACTIVE']= 1.0

        params['ACTIVE_DynDepth']  =1 
        params['ACTIVE_TotDepth']  =1
        #params['ACTIVE_ATMDEP_NO3']=1
        params['ACTIVE_Nitrif_NH4']=1
        # disable regular nitrification
        params['RcNit']=0.0
        params['RcNit20']=0.0

        # params['fAtmDepNO3']=0.01 # g m-2 d-1
        polys=wkb2shp.shp2geom('../../dflowfm/gis/poly-features.shp')
        for poly in polys:
            if poly['name']=='calhoun_marsh':
                g=self.hydro.grid()
                cells=g.select_cells_intersecting(poly['geom'])

                cell_rates=np.zeros(g.Ncells())
                # for atm dep: g m-2 d-1
                # for nitrif: mg/l d-1
                # hmm - I'm getting a max value of 11.6 over a two day run.
                #  growth in one of those vernal ponds is linear, so that
                #  much is good.
                cell_rates[cells] = 0.01 
                self.hydro.infer_2d_elements()
                seg_rates=cell_rates[self.hydro.seg_to_2d_element]

                # 001:
                # params['fAtmDepNO3']=waq.ParameterSpatial(per_segment=seg_rates)
                # 002:
                params['ZNit']=waq.ParameterSpatial(per_segment=seg_rates)
                params['NH4']=10.0 # always plenty of NH4
                print("Set rate in %d/%d elements"%(cells.sum(),g.Ncells()))

        if self.hydro.n_exch_z>0:
            # if we don't have vertical diffusion specified, this will
            # raise an error, but that's good.  3D runs should have it.
            params['ACTIVE_VertDisp']= 1.0

        return params

    def init_substances(self):
        subs=super(Scen,self).init_substances()
        subs['NO3']=Sub(0)

        return subs

    #def init_bcs(self):
    #def init_loads(self):

hydro=waq.HydroFiles("../../dflowfm/runs/20180807_grid98_16_single/DFM_DELWAQ_flowfm/flowfm.hyd")

# if it runs all the way to the end of the hydro, the bad fluxes on the last step
# appear to cause a gmres error.
scen=Scen(hydro)  # this is still setup

scen.delft_bin="/home/rusty/src/dfm/1.4.6/bin"
os.environ['LD_LIBRARY_PATH']="/home/rusty/src/dfm/1.4.6/lib"
scen.share_path="/home/rusty/src/dfm/1.4.6/share/delft3d"
scen.map_output += ('TotalDepth','Depth','LocalDepth')

os.path.exists(scen.base_path) and shutil.rmtree(scen.base_path) # for dev

scen.cmd_default() # actually write delwaq's input files
scen.cmd_delwaq1() # delwaq1 is the preprocessor - compiles the water quality processes
scen.cmd_delwaq2() # actually runs the simulation
scen.write_binary_map_nc() # translate simulation output from proprietary binary=>nc


##

# what's a good tracer that I can set a 0th order accumulation?
# what NO3, and use atmospheric deposition for the source?

