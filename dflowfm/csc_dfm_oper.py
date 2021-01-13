#!/usr/bin/env python
"""
Script to configure and generate mdu ready for running dflow fm
for Cache Slough Complex
 - *_decker version switches to a grid that extends down to decker island,
   necessitating different forcing config.
"""

import os, sys
import numpy as np
import logging
import xarray as xr
import six
import shutil

from stompy import utils, filters
import stompy.model.delft.io as dio
from stompy.io.local import usgs_nwis
from stompy.model.delft import dfm_grid
from stompy.grid import unstructured_grid
from stompy.spatial import wkb2shp, field

import gen_polygons
        
if __name__=='__main__':
    logging.basicConfig(level=logging.INFO)

log=logging.getLogger('csc_dfm')

try:
    import local_config
except ImportError:
    log.error("local_config.py not found. May need to copy from local_config.py.in")
    sys.exit(1)

import barker_data

import stompy.model.delft.dflow_model as dfm
six.moves.reload_module(local_config)
import stompy.model.hydro_model as hm
import csc_dfm_decker_roughsearch as rs

try:
    here=os.path.dirname(__file__)
except NameError:
    here="."
    log.info("Assuming script is in %s"%here)

class CscDeckerModel(dfm.DFlowModel):
    cache_dir='cache'
    z_datum='NAVD88' # not really used right now.
    projection='EPSG:26910'
    # Forcing data is fetched as UTC, and adjusted according to this offset
    utc_offset=np.timedelta64(-8,'h') # PST.  

    src_grid_fn=os.path.join(here,'../grid/CacheSloughComplex_v111-edit21.nc')

    def load_default_mdu(self):
        self.load_mdu('template.mdu')
        
    def ccs_pre_restoration(self):
        return self.run_start < np.datetime64("2014-11-01")

    @property
    def tidal_bc_location(self):
        if self.run_start < np.datetime64("2016-05-01"):
            return 'riovista'
        else:
            return 'decker'
    
    def get_grid(self):
        """ 
        Get a grid with bathy. Depends on run_start in order to choose
        the right bathymetry.
        """
        if self.ccs_pre_restoration():
            bathy_fn=os.path.join(here,"../bathy/merged-20190530-pre_calhoun.tif")
            dst_grid=os.path.basename(self.src_grid_fn).replace('.nc','-pre-bathy.nc')
        else:    
            bathy_fn=os.path.join(here,"../bathy/merged_2m-20190122.tif")
            dst_grid=os.path.basename(self.src_grid_fn).replace('.nc','-bathy.nc')

        if utils.is_stale(dst_grid,[self.src_grid_fn,bathy_fn],ignore_missing=True):
            g=unstructured_grid.UnstructuredGrid.from_ugrid(self.src_grid_fn)
            dem=field.GdalGrid(bathy_fn)
            import dem_cell_node_bathy
            node_depths=dem_cell_node_bathy.dem_to_cell_node_bathy(dem,g)
            g.add_node_field('depth',node_depths,on_exists='overwrite')
            g.write_ugrid(dst_grid,overwrite=True)
        else:
            g=unstructured_grid.UnstructuredGrid.from_ugrid(dst_grid)

        if self.tidal_bc_location=='riovista':
            # truncate domain:
            srv_line=model.get_geometry(name='SRV',geom_type='LineString')
            to_keep=g.select_cells_by_cut(srv_line)
            for c in np.nonzero(~to_keep)[0]:
                g.delete_cell(c)
            g.delete_orphan_edges()
            g.delete_orphan_nodes()
            g.renumber()
            
        return g

    def setup_bcs(self):
        dredge_depth=-1

        # check_bspp.py has the code that converted original tim to csv.
        barker=barker_data.BarkerPumpsBC(name='Barker_Pumping_Plant',dredge_depth=dredge_depth)
        self.add_bcs(barker)

        if self.tidal_bc_location=='decker':
            # Decker only exists post-2015
            self.add_bcs(hm.NwisStageBC(name='decker',station=11455478,cache_dir=self.cache_dir,
                                        filters=[hm.Lowpass(cutoff_hours=1.0)]))
        elif self.tidal_bc_location=='riovista':
            self.add_bcs(hm.NwisStageBC(name='SRV',station=11455420,cache_dir=self.cache_dir,
                                        filters=[hm.Lowpass(cutoff_hours=1.0)]))
        else:
            raise Exception("Bad value for tidal_bc_location: %s"%self.tidal_bc_location)

        if self.tidal_bc_location=='decker':
            # unclear whether threemile should also be flipped.  mean flows typically Sac->SJ,
            # and the test period shows slightly negative means, so it's possible that it is
            # correct.
            # flipping this did improve a lot of phases, but stage at TSL is much worse, and
            # my best reading of the metadata is that it should not be flipped.
            self.add_bcs(hm.NwisFlowBC(name='threemile',station=11337080,cache_dir=self.cache_dir,
                                       dredge_depth=dredge_depth,
                                       filters=[hm.Lowpass(cutoff_hours=1.0)]))

        # GSS: from compare_flows and lit, must be flipped.
        self.add_bcs(hm.NwisFlowBC(name='Georgiana',station=11447903,cache_dir=self.cache_dir,
                                   dredge_depth=dredge_depth,
                                   filters=[hm.Transform(lambda x: -x),
                                            hm.Lowpass(cutoff_hours=1.0)] ))

        self.add_bcs(hm.NwisFlowBC(name='dcc',station=11336600,cache_dir=self.cache_dir,
                                   default=0.0,dredge_depth=dredge_depth,
                                   filters=[hm.FillGaps(large_gap_value=0.0),
                                            hm.Lowpass(cutoff_hours=1.0),
                                            hm.Transform(lambda x: -x)] ) )

        # moving Sac flows upstream and removing tides.
        sac=hm.NwisFlowBC(name="SacramentoRiver",station=11447650,
                          pad=np.timedelta64(5,'D'),cache_dir=self.cache_dir,
                          dredge_depth=dredge_depth,
                          filters=[hm.LowpassGodin(),
                                   hm.Lag(np.timedelta64(-2*3600,'s'))])
        self.add_bcs(sac)

        if 1:
            lisbon_bc=hm.CdecFlowBC(name='lis',station="LIS",pad=np.timedelta64(5,'D'),
                                    default=0.0,cache_dir=self.cache_dir,
                                    dredge_depth=dredge_depth,
                                    filters=[hm.FillGaps(),
                                             hm.Lowpass(cutoff_hours=1.0)])
            self.add_bcs(lisbon_bc)
        else:
            log.warning("TEMPORARILY Disabling Lisbon flow due to CDEC issues")

        # Ulatis inflow
        # There are probably timezone issues here - they are coming in PST, but
        # this code probably assumes UTC.
        ulatis_ds=xr.open_dataset(os.path.join(here,"../bcs/ulatis/ulatis_hwy113.nc"))
        ulatis_ds['flow']=ulatis_ds.flow_cfs*0.02832
        ulatis_ds['flow'].attrs['units']='m3 s-1'

        def pad_with_zero(ds,pad=np.timedelta64(1,'D')):
        # if ( (self.run_start>=ulatis_ds.time.min()) and
            if self.run_stop+pad>ds.time[-1]:
                log.warning("Will extend flow with 0 flow! (data end %s)"%(ds.time.values[-1]))
                # with xarray, easier to just overwrite the last sample.  lazy lazy.
                ds.time.values[-1] = self.run_stop+pad
                ds.flow.values[-1] = 0.0
            if self.run_start-pad<ds.time[0]:
                log.warning("Will prepend flow with 0 flow! (data starts %s)"%(ds.time.values[0]))
                # with xarray, easier to just overwrite the last sample.  lazy lazy.
                ds.time.values[0] = self.run_start - pad
                ds.flow.values[0] = 0.0
            return ds

        ulatis_ds=pad_with_zero(ulatis_ds)
        self.add_bcs(hm.FlowBC(name='ULATIS',flow=ulatis_ds.flow,dredge_depth=dredge_depth))

        # Campbell Lake
        campbell_ds=xr.open_dataset(os.path.join(here,"../bcs/ulatis/campbell_lake.nc"))
        campbell_ds['flow']=campbell_ds.flow_cfs*0.02832
        campbell_ds['flow'].attrs['units']='m3 s-1'
        campbell_ds=pad_with_zero(campbell_ds)
        self.add_bcs(hm.FlowBC(name='CAMPBELL',flow=campbell_ds.flow,dredge_depth=dredge_depth))

        if 0: # not including wind right now
            # Try file pass through for forcing data:
            windxy=self.read_tim('forcing-data/windxy.tim',columns=['wind_x','wind_y'])
            windxy['wind_xy']=('time','xy'),np.c_[ windxy['wind_x'].values, windxy['wind_y'].values]
            self.add_WindBC(wind=windxy['wind_xy'])

        self.setup_roughness()
        
    def setup_roughness(self):
        # Roughness 
        if 1:
            # These are setting that came out of the optimization
            settings={}
            settings['cache']=0.04
            settings['dws']=0.025
            settings['elk']=0.015
            settings['fpt_to_dcc']=0.030
            settings['lindsey']=0.04
            settings['miner']=0.035
            settings['rio_vista']=0.025
            settings['sac_below_ges']=0.0225
            settings['steamboat']=0.025
            settings['toe']=0.0375
            settings['upper_sac']=0.0175
            xyn=rs.settings_to_roughness_xyz(self,settings)
            # Turn that into a DataArray
            da=xr.DataArray( xyn[:,2],dims=['location'],name='n' )
            da=da.assign_coords(x=xr.DataArray(xyn[:,0],dims='location'),
                                y=xr.DataArray(xyn[:,1],dims='location'))
            da.attrs['long_name']='Manning n'
            rough_bc=hm.RoughnessBC(data_array=da)
            self.add_bcs(rough_bc)

    def setup_structures(self):
        # Culvert at CCS
        # rough accounting of restoration
        if self.ccs_pre_restoration:
            # name => id
            # polylinefile is filled in from shapefile and name

            # these structure parameters get an amplitude ratio at CCS of 1.241
            # the shape is not great, and suggests that as the stage gets close to 1.7
            # or so, wetted area expands rapidly. That is probably a shortcoming of the
            # current bathymetry.
            self.add_Structure(name='ccs_breach',
                               type='gate',
                               door_height=15, # no overtopping?
                               lower_edge_level=0.9,
                               # in release 1.5.2, this needs to be nonzero.  Could use
                               # width of 0.0 in some older versions, but no longer.
                               opening_width=0.1,
                               sill_level=0.8,
                               horizontal_opening_direction = 'symmetric')
            # these are really just closed
            for levee_name in ['ccs_west','ccs_east']:
                self.add_Structure(name=levee_name,
                                   type='gate',
                                   door_height=15,
                                   lower_edge_level=3.5,
                                   opening_width=0.0,
                                   sill_level=3.5,
                                   horizontal_opening_direction = 'symmetric')

    def setup_monitoring(self):
        # -- Extract locations for sections and monitor points
        mon_sections=self.match_gazetteer(monitor=1,geom_type='LineString')
        mon_points  =self.match_gazetteer(geom_type='Point')
        self.add_monitor_sections(mon_sections)
        self.add_monitor_points(mon_points)

    def __init__(self,*a,**kw):
        super(CscDeckerModel,self).__init__(*a,**kw)
        
        # Register linear and point feature shapefiles
        self.add_gazetteer('gis/model-features.shp')
        self.add_gazetteer('gis/point-features.shp')

        self.set_grid(self.get_grid())

        # 6 is maybe better for getting good edges
        self.mdu['geometry','BedlevType']=6
        self.mdu['geometry','Conveyance2D']=-1
        self.mdu['geometry','Nonlin2D']=0

        self.mdu['physics','UnifFrictCoef']= 0.023
        # fail out when it goes unstable.
        self.mdu['numerics','MinTimestepBreak']=0.05
        if 0: # For PTM usage
            self.mdu['output','WaqInterval']=1800 
            self.mdu['output','MapInterval']=1800
        else:
            self.mdu['output','MapInterval']=3600

        self.mdu['output','RstInterval']=86400
            
        self.create_with_mode(self.cache_dir,mode='create')

        # -- Boundary Conditions --
        # All sources and flow BCs will be dredged to this depth to make
        # they remain wet and active.
        self.setup_bcs()
        
        self.setup_structures() # 1 call
        self.setup_monitoring()

    def write(self):
        super(CscDeckerModel,self).write()
        
        gen_polygons.gen_polygons(self.run_dir)
        self.write_bc_plots()

        try:
            script=__file__
        except NameError:
            script=None
        if script:
            shutil.copyfile(script,
                            os.path.join( os.path.join(model.run_dir,
                                                       os.path.basename(script) ) ))
        
    def write_bc_plots(self):
        # experimental saving of BC data to html plots
        for bc in self.bcs:
            bc.write_bokeh(path=self.run_dir)

## 
# if not invoked directly, just set up the model and let
# the importer decide what to do with model.

if __name__=='__main__':
    import argparse

    # v00: pre-restoration period with post-restoration grid
    # v01: use pre-restoration bathy
    # v02: tuned CCS culvert
    # v03: Add Ulatis and Campbell Lake inflows
    
    parser = argparse.ArgumentParser(description='Set up and run Cache Slough Complex simulation.')
    parser.add_argument("-s", "--start", help="Date of simulation start",
                        default="2014-02-25")
    parser.add_argument("-e", "--end", help="Date of simulation stop",
                        default="2014-05-15")
    parser.add_argument("-d", "--dir", help="Run directory",
                        default=None,required=True)
    parser.add_argument("-r", "--resume", help="Resume from run",
                        default=None)
    parser.add_argument("-n", "--dryrun", help="Do not actually partition or run the simulation",
                        action='store_true')
    parser.add_argument("-N", "--no-run", help="Write and partition but do not run the simulation",
                        action='store_true')
    parser.add_argument("--interval",help="Interval for multiple shorter runs, e.g. 1D for 1 day",
                        default=None)

    args=None
    # Uncomment this to hardwire in a command line, useful for testing from within
    # a python console session.
    # args="-d runs/v03/v03 -s 2019-01-15T00:00 -e 2019-07-05T00:00 --interval 30D".split()
    args=parser.parse_args(args=args)
    
    if args.resume is not None:
        # Go ahead and get the restart time, so that intervals and directory names
        # can be chosen below
        last_run_dir=args.resume
        last_model=drv.SuntansModel.load(last_run_dir)
        multi_run_start=last_model.restartable_time()
        print("Will resume run in %s from %s"%(last_run_dir,multi_run_start))
    else:
        multi_run_start=np.datetime64(args.start)
        print("Run start: ",multi_run_start)
        last_model=None
        last_run_dir=None
        
    run_start=multi_run_start
    multi_run_stop=np.datetime64(args.end)

    # series of 1-day runs
    # run_interval=np.timedelta64(1,'D')
    # But restarts are not good with average output.
    # Kludge it and run the whole period in one go.
    if args.interval is not None:
        # annoying that this can only be integer values.
        # possible inputs:
        # 5D for 5 days
        # 12h for 12 hours
        # 1800s for half an hour
        run_interval=np.timedelta64(int(args.interval[:-1]),args.interval[-1])
    else:
        # in one go.
        run_interval=multi_run_stop-run_start

    assert args.dir is not None
        
    run_count=0
    while True:
        if (last_run_dir is not None) and (last_model is None):
            # this would happen if a restart wasn't specifically
            # requested, it's a multi-run setup (interval was set)
            # and some of the earlier runs have already completed.
            last_model=CscDeckerModel.load(last_run_dir)

        if last_model is not None:
            run_start=last_model.restartable_time()

        if run_start >= multi_run_stop:
            break
        
        run_stop=run_start+run_interval
        print("Simulation period: %s -- %s"%(run_start,run_stop))
        date_str=utils.to_datetime(run_start).strftime('%Y%m%d')
        
        run_dir=f"{args.dir}_{date_str}"

        if not CscDeckerModel.run_completed(run_dir):
            model=CscDeckerModel(run_dir=run_dir,
                                 run_start=run_start,
                                 run_stop=run_stop,
                                 restart_from=last_model)
            if args.dryrun:
                print("Dry run - dropping out of loop")
                break
            
            model.write()
            model.partition()
            if args.no_run:
                print("No run - dropping out of loop")
                break
            model.run_model()
            if not model.is_completed():
                log.error("Breaking out of loop -- run %s did not complete"%run_dir)
                break
            last_model=model
        else:
            print("  -- appears to have already been run")
            last_model=CscDeckerModel.load(run_dir)
            
        run_count+=1
        last_run_dir=run_dir
        # not precisely. The restartable time may be different than
        # run_stop, based on frequency of writing restart data.
        run_start=run_stop

