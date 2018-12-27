#!/usr/bin/env python
"""
Script to configure and generate mdu ready for running dflow fm
for Cache Slough Complex
 - *_decker version switches to a grid that extends down to decker island,
   necessitating different forcing config.
 - this version slowly fills the domain to extract effective hypsometry
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

import stompy.model.delft.dflow_model as dfm
cache_dir='cache'

## --------------------------------------------------

six.moves.reload_module(dfm_grid)
six.moves.reload_module(dfm)
six.moves.reload_module(dio)
dfm_bin_dir=os.path.join(os.environ['HOME'],"src/dfm/1.5.0/bin")
os.environ['LD_LIBRARY_PATH']=dfm_bin_dir.replace('bin','lib')
dfm.DFlowModel.dfm_bin_dir=dfm_bin_dir
dfm.DFlowModel.mpi_bin_dir=dfm_bin_dir

model=dfm.DFlowModel()

model.num_procs=1 
model.z_datum='NAVD88'
model.projection='EPSG:26910'
model.utc_offset=np.timedelta64(-8,'h') # PST

# 00: initial
# 01: fix grid error
# 02: fix location of sections
# 03: fix initial condition
# 04: conveyance2d=-1
# 05: bedlevtype=6
model.set_run_dir("runs/hypso_05", mode='askclobber')

model.run_start=np.datetime64('2017-12-01')
model.run_stop=np.datetime64('2017-12-10')

model.load_mdu('template.mdu')

src_grid='../grid/CacheSloughComplex_v111-edit19.nc'
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
# 6 means cells=mean(nodes), edges=shallower(cells)
model.mdu['geometry','BedlevType']=6
model.mdu['geometry','Conveyance2D']=-1
# enabling this seems to cause a lot of oscillation.
# but it may be necessary to get good prism on the narrow channels.
model.mdu['geometry','Nonlin2D']=0
model.mdu['output','MapInterval']=7200
# low friction to let it fill quickly and evenly
model.mdu['physics','UnifFrictCoef']= 0.01
# fail out when it goes unstable.
model.mdu['numerics','MinTimestepBreak']=0.05

model.set_cache_dir('cache')

model.add_gazetteer('gis/model-features.shp')
model.add_gazetteer('gis/point-features.shp')

pad=np.timedelta64(15*60,'s')
bc_time=np.arange(model.run_start,model.run_stop+pad,np.timedelta64(1,'h'))
bc_h=np.linspace(-2,4,len(bc_time))

da=xr.DataArray(bc_h,dims=['time'], coords={'time':bc_time})

model.add_bcs(dfm.StageBC( name='decker',
                           z=da ))

model.mdu['geometry','WaterLevIni']=bc_h[0]

mon_sections=model.match_gazetteer(monitor=1,geom_type='LineString')
mon_points  =model.match_gazetteer(geom_type='Point')
model.add_monitor_sections(mon_sections)
model.add_monitor_points(mon_points)


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
