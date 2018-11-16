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
from stompy.spatial import wkb2shp

if __name__=='__main__':
    logging.basicConfig(level=logging.INFO)

log=logging.getLogger('csc_dfm')

import barker_data
import nwis_bc
import stompy.model.delft.dflow_model as dfm

## --------------------------------------------------

six.moves.reload_module(dfm_grid)
six.moves.reload_module(dfm)
six.moves.reload_module(dio)
six.moves.reload_module(barker_data)
six.moves.reload_module(nwis_bc)

model=dfm.DFlowModel()

# Having issues with 53925-opt, and 52184-dbg, both
# with MPI.  Looks like dwaq output is not compatible
# with ugrid+mpi.
model.dfm_bin_dir="/home/rusty/src/dfm/r53925-opt/bin"
model.num_procs=4
model.z_datum='NAVD88'
model.projection='EPSG:26910'
model.utc_offset=np.timedelta64(-8,'h') # PST

# Parameters to control more specific aspects of the run
# grid100_00: First run with post-restoration grid, new features in HydroModel.
#   Fix timezone, previous run looks about 8h ahead.  Running again.
model.set_run_dir("runs/grid100_00", mode='pristine')

model.run_start=np.datetime64('2018-08-01')
model.run_stop=np.datetime64('2018-09-01')

if 0:
    # this switches to low-biased edge depths, still optimized, and
    # uses a master DEM with CCS cut down a bit.
    model.set_grid("CacheSloughComplex_v100_bathy2_sparse_net.nc")
if 1: # Instead, try the bedlevtype=2 with straight node bathy
    src_grid='../grid/CacheSloughComplex_v108.nc'
    dst_grid="CacheSloughComplex_v108-bathy.nc"
    bathy_fn="../bathy/merged_2m-20181113.tif"
    if utils.is_stale(dst_grid,[src_grid,bathy_fn]):
        g=unstructured_grid.UnstructuredGrid.from_ugrid(src_grid)
        dem=field.GdalGrid(bathy_fn)
        g.add_node_field('depth',dem(g.nodes['x']),on_exists='overwrite')
        g.write_ugrid(dst_grid,overwrite=True)
    else:
        g=unstructured_grid.UnstructuredGrid.from_ugrid(dst_grid)
    model.set_grid(g)

model.load_mdu('template.mdu')
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
model.add_bcs(nwis_bc.NwisFlowBC(name='threemile',station=11337080,cache_dir='cache'))
model.add_bcs(nwis_bc.NwisFlowBC(name='Georgiana',station=11447903,cache_dir='cache'))
model.add_bcs(nwis_bc.NwisFlowBC(name='dcc',station=11336600,cache_dir='cache'))

sac=nwis_bc.NwisFlowBC(name="SacramentoRiver",station=11447650,
                       pad=np.timedelta64(5,'D'),cache_dir='cache',
                       filters=[dfm.LowpassGodin(),
                                dfm.Lag(np.timedelta64(-2*3600,'s'))])
model.add_bcs(sac)

if 0: # don't have wind data for newer period
    # Try file pass through for forcing data:
    windxy=model.read_tim('forcing-data/windxy.tim',columns=['wind_x','wind_y'])
    windxy['wind_xy']=('time','xy'),np.c_[ windxy['wind_x'].values, windxy['wind_y'].values]
    model.add_WindBC(wind=windxy['wind_xy'])
    
# Roughness
# 0.020 on the Sac, 0.023 in CSC
# also 0.2(!) in the marsh off CCS
model.add_RoughnessBC(shapefile='forcing-data/manning_slick_sac.shp')

# Culvert at CCS
if 1:
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

