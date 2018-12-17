#!/usr/bin/env python
"""
Script to configure and generate mdu ready for running dflow fm
for Cache Slough Complex
 - *_decker version switches to a grid that extends down to decker island,
   necessitating different forcing config.
"""

import os
import shutil
import numpy as np
import logging
import pandas as pd
import xarray as xr
import six

from stompy import utils, filters
import stompy.model.delft.io as dio
from stompy.io.local import usgs_nwis
from stompy.model.delft import dfm_grid
from stompy.grid import unstructured_grid
from stompy.spatial import wkb2shp, field

if __name__=='__main__':
    logging.basicConfig(level=logging.INFO)

log=logging.getLogger('csc_dfm')

import barker_data
import nwis_bc
import stompy.model.delft.dflow_model as dfm
cache_dir='cache'

## --------------------------------------------------

six.moves.reload_module(dfm_grid)
six.moves.reload_module(dfm)
six.moves.reload_module(dio)
six.moves.reload_module(barker_data)
six.moves.reload_module(nwis_bc)
import local_config
local_config.install()

model=dfm.DFlowModel()


# Having issues with 53925-opt, and 52184-dbg, both
# with MPI.  Looks like dwaq output is not compatible
# with ugrid+mpi.
model.num_procs=4
model.z_datum='NAVD88'
model.projection='EPSG:26910'
model.utc_offset=np.timedelta64(-8,'h') # PST

# Parameters to control more specific aspects of the run
# grid100_00: First run with post-restoration grid, new features in HydroModel.
#   Fix timezone, previous run looks about 8h ahead.  Running again.

# grid100_01: bedlevtyp=3, conveyance2d=3, nonlin2d=1
# grid100_02: fix sign for DCC and Georgiana flows.
#        _03: conveyanced2d=2.  Basically same results, and 80% faster.
# grid100_04: bedlevtyp=3, conveyance2d=1, nonlin2d=0
# grid100_05: bedlevtyp=3, conveyance2d=0, nonlin2d=0
# grid100_06: bedlevtyp=3, conveyance2d=2, nonlin2d=0
#          drop the roughness shapefile, set global 0.025
# grid100_07: back to same settings as _03, but with new grid that
#      slims DWS and parts of DWB (Miner slough?)
#        _08: fix some bathy, add sections.
#        _09: more sections, dial down DWS friction to 0.02, dial up friction
#             north of Sutter
#        _10: new bit of stairstep grid, more sections, increase friction in
#             SUT, HWB, SSS
#        _11: even lower friction on DWS
#        _12: revert friction on DWS, bump up friction on dead-end north of SUT,
#             and drop friction on lower sac to match rest of sac.
#        _13: return to Thomas' roughness for a reality check.
#        _14: run a December 2017 period, see if phasing is better when DCC is closed.
#        _15: same, but flip TSL to see if it makes a difference in cal.
model.set_run_dir("runs/grid100_15", mode='askclobber')

#model.run_start=np.datetime64('2018-08-01')
#model.run_stop=np.datetime64('2018-09-01')
model.run_start=np.datetime64('2017-12-01')
model.run_stop=np.datetime64('2017-12-10')

model.load_mdu('template.mdu')

src_grid='../grid/CacheSloughComplex_v111-edit01.nc'
dst_grid=os.path.basename(src_grid).replace('.nc','-bathy.nc')
bathy_fn="../bathy/merged_2m-20181113.tif"
if utils.is_stale(dst_grid,[src_grid,bathy_fn]):
    g=unstructured_grid.UnstructuredGrid.from_ugrid(src_grid)
    dem=field.GdalGrid(bathy_fn)
    node_depths=dem(g.nodes['x'])
    while np.any(np.isnan(node_depths)):
        missing=np.nonzero(np.isnan(node_depths))[0]
        print("Looping to fill in %d missing node depths"%len(missing))
        for n in missing:
            nbrs=g.node_to_nodes(n)
            node_depths[n]=np.nanmean(node_depths[nbrs])
    g.add_node_field('depth',node_depths,on_exists='overwrite')
    g.write_ugrid(dst_grid,overwrite=True)
else:
    g=unstructured_grid.UnstructuredGrid.from_ugrid(dst_grid)
model.set_grid(g)

# I think this doesn't matter with conveyance2d>0
model.mdu['geometry','BedlevType']=3
# ostensibly analytic 2D conveyance with 3.
# 2 is faster, and appears less twitchy while showing very similar
#  calibration results.
# 0 blocks some flow, 1 was a little twitchy.
model.mdu['geometry','Conveyance2D']=2
# enabling this seems to cause a lot of oscillation.
model.mdu['geometry','Nonlin2D']=0

model.mdu['output','MapInterval']=7200
#model.mdu['output','WaqInterval']=1800
model.mdu['physics','UnifFrictCoef']= 0.023

model.set_cache_dir('cache')

model.add_gazetteer('gis/model-features.shp')
model.add_gazetteer('gis/point-features.shp')

dfm.SourceSinkBC.dredge_depth=-1
dfm.FlowBC.dredge_depth=-1

# check_bspp.py has the code that converted original tim to csv.
# awkward reloading of that, but at least it's independent of the
# simulation period.
model.add_bcs(barker_data.BarkerPumpsBC(name='Barker_Pumping_Plant'))

six.moves.reload_module(nwis_bc)

# Decker only exists post-2015
if model.run_start>np.datetime64("2015-11-16"):
    model.add_bcs(nwis_bc.NwisStageBC(name='decker',station=11455478,cache_dir='cache'))
else:
    # maybe fall back to Rio Vista, or some adjustment thereof
    raise Exception("Decker tidal data starts 2015-11-16, too late for this simulation period")
# unclear whether threemile should also be flipped.  mean flows typically Sac->SJ,
# and the test period shows slightly negative means, so it's possible that it is
# correct.
# flipping this did improve a lot of phases, but stage at TSL is much worse, and
# my best reading of the metadata is that it should not be flipped.
model.add_bcs(nwis_bc.NwisFlowBC(name='threemile',station=11337080,cache_dir='cache'))
# GSS: from compare_flows and lit, must be flipped.
model.add_bcs(nwis_bc.NwisFlowBC(name='Georgiana',station=11447903,cache_dir='cache',
                                 filters=[dfm.Transform(lambda x: -x)] ))

class FillGaps(dfm.BCFilter):
    max_gap_interp_s=2*60*60
    large_gap_value=0.0
    
    def transform_output(self,da):
        # have self.bc, self.bc.model
        # self.bc.data_start, self.bc.data_stop
        if len(da)==0:
            log.warning("FillGaps called with no input data")
            da=xr.DataArray(self.large_gap_value)
            return da
        log.warning("FillGaps code incomplete")
        return da

model.add_bcs(nwis_bc.NwisFlowBC(name='dcc',station=11336600,cache_dir='cache',
                                 filters=[FillGaps(),
                                          dfm.Transform(lambda x: -x)] ) )

sac=nwis_bc.NwisFlowBC(name="SacramentoRiver",station=11447650,
                       pad=np.timedelta64(5,'D'),cache_dir='cache',
                       filters=[dfm.LowpassGodin(),
                                dfm.Lag(np.timedelta64(-2*3600,'s'))])
model.add_bcs(sac)

## 
from stompy.io.local import cdec
pad=np.timedelta64(5,'D')
lisbon_ds=cdec.cdec_dataset(station='LIS',
                            start_date=model.run_start-pad, end_date=model.run_stop+pad,
                            sensor=20, cache_dir=cache_dir)
# to m3/s
lisbon_ds['Q']=lisbon_ds['sensor0020'] * 0.028316847
# to hourly average
lisbon_ds.groupby( lisbon_ds.time.astype('M8[h]')).mean()
lisbon_bc=dfm.FlowBC(name='lis',Q=lisbon_ds.Q)

## 
if 0: # don't have wind data for newer period
    # Try file pass through for forcing data:
    windxy=model.read_tim('forcing-data/windxy.tim',columns=['wind_x','wind_y'])
    windxy['wind_xy']=('time','xy'),np.c_[ windxy['wind_x'].values, windxy['wind_y'].values]
    model.add_WindBC(wind=windxy['wind_xy'])
    
# Roughness
# model.add_RoughnessBC(shapefile='forcing-data/manning_slick_sac.shp')
model.add_RoughnessBC(shapefile='forcing-data/manning_n.shp')

# Culvert at CCS
if 1:
    # rough accounting of restoration
    if model.run_start < np.datetime64("2014-11-01"):
        # name => id
        # polylinefile is filled in from shapefile and name
        model.add_Structure(name='ccs_breach',
                            type='gate',
                            door_height=15, # no overtopping?
                            lower_edge_level=0.89,
                            opening_width=0.0, # pretty sure this is ignored.
                            sill_level=0.85, # gives us a 0.05m opening?
                            horizontal_opening_direction = 'symmetric')
        # these are really just closed
        for levee_name in ['ccs_west','ccs_east']:
            model.add_Structure(name=levee_name,
                                type='gate',
                                door_height=15,
                                lower_edge_level=3.5,
                                opening_width=0.0,
                                sill_level=3.5,
                                horizontal_opening_direction = 'symmetric')

mon_sections=model.match_gazetteer(monitor=1,geom_type='LineString')
mon_points  =model.match_gazetteer(geom_type='Point')
model.add_monitor_sections(mon_sections)
model.add_monitor_points(mon_points)

## 
if 1:
    for bc in model.bcs:
        bc.write_bokeh(path=model.run_dir)

##

if __name__=='__main__':
    model.write()
    model.partition()
    try:
        script=__file__
    except NameError:
        script=None
    if script:
        shutil.copyfile(script,
                        os.path.join( os.path.join(model.run_dir,
                                                   os.path.basename(script) ) ))

    model.run_model()

# 2018-11-30
# trying to understand whether the GES/DCC reach is correct
#  GES: biased negative (floodward), in particular about 50 m3/s too strong on flood
#        and slightly weak on ebb.  mean delta 24m3/s, 7 minute lag.
#  the Georgiana flows oscillate between 50 and 100, presumably drawing water out.
