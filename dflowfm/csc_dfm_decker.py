#!/usr/bin/env python
"""
Script to configure and generate mdu ready for running dflow fm
for Cache Slough Complex
 - *_decker version switches to a grid that extends down to decker island,
   necessitating different forcing config.
"""

import os
import numpy as np
import logging
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

import stompy.model.delft.dflow_model as dfm
cache_dir='cache'

try:
    here=os.path.dirname(__file__)
except NameError:
    here="."
    log.info("Assuming script is in %s"%here)

## --------------------------------------------------

import local_config
local_config.install()

model=dfm.DFlowModel()

model.num_procs=12
model.z_datum='NAVD88' # not really used right now.
model.projection='EPSG:26910'
# Forcing data is fetched as UTC, and adjusted according to this offset
model.utc_offset=np.timedelta64(-8,'h') # PST.  

# v00: pre-restoration period with post-restoration grid
# v01: use pre-restoration bathy
# v02: tuned CCS culvert
# v03: Add Ulatis and Campbell Lake inflows
model.set_run_dir("runs/val-201404_v03", mode='pristine')

model.run_start=np.datetime64('2014-02-25')
model.run_stop=np.datetime64('2014-05-15')

model.load_mdu('template.mdu')

# Register linear and point feature shapefiles
model.add_gazetteer('gis/model-features.shp')
model.add_gazetteer('gis/point-features.shp')

ccs_pre_restoration=model.run_start < np.datetime64("2014-11-01")

src_grid='../grid/CacheSloughComplex_v111-edit19fix.nc'
if ccs_pre_restoration:
    bathy_fn="../bathy/merged-20190530-pre_calhoun.tif"
    dst_grid=os.path.basename(src_grid).replace('.nc','-pre-bathy.nc')
else:    
    bathy_fn="../bathy/merged_2m-20190122.tif"
    dst_grid=os.path.basename(src_grid).replace('.nc','-bathy.nc')
    
if utils.is_stale(dst_grid,[src_grid,bathy_fn],ignore_missing=True):
    g=unstructured_grid.UnstructuredGrid.from_ugrid(src_grid)
    dem=field.GdalGrid(bathy_fn)
    import dem_cell_node_bathy
    node_depths=dem_cell_node_bathy.dem_to_cell_node_bathy(dem,g)
    g.add_node_field('depth',node_depths,on_exists='overwrite')
    g.write_ugrid(dst_grid,overwrite=True)
else:
    g=unstructured_grid.UnstructuredGrid.from_ugrid(dst_grid)

tidal_bc_location='decker'

if model.run_start < np.datetime64("2016-05-01"):
    tidal_bc_location='riovista'

if tidal_bc_location=='riovista':
    # truncate domain:
    srv_line=model.get_geometry(name='SRV',geom_type='LineString')
    to_keep=g.select_cells_by_cut(srv_line)
    for c in np.nonzero(~to_keep)[0]:
        g.delete_cell(c)
    g.delete_orphan_edges()
    g.delete_orphan_nodes()
    g.renumber()


##    
    
model.set_grid(g)

# 6 is maybe better for getting good edges
model.mdu['geometry','BedlevType']=6
model.mdu['geometry','Conveyance2D']=-1
model.mdu['geometry','Nonlin2D']=0

model.mdu['physics','UnifFrictCoef']= 0.023
# fail out when it goes unstable.
model.mdu['numerics','MinTimestepBreak']=0.05
# For PTM usage
model.mdu['output','WaqInterval']=1800 
model.mdu['output','MapInterval']=1800

# Make sure cache dir exists
os.path.exists(cache_dir) or os.makedirs(cache_dir)


# -- Boundary Conditions --

# All sources and flow BCs will be dredged to this depth to make
# they remain wet and active.
dfm.SourceSinkBC.dredge_depth=-1
dfm.FlowBC.dredge_depth=-1

# check_bspp.py has the code that converted original tim to csv.
model.add_bcs(barker_data.BarkerPumpsBC(name='Barker_Pumping_Plant'))

if tidal_bc_location=='decker':
    # Decker only exists post-2015
    model.add_bcs(dfm.NwisStageBC(name='decker',station=11455478,cache_dir='cache'))
elif tidal_bc_location=='riovista':
    model.add_bcs(dfm.NwisStageBC(name='SRV',station=11455420,cache_dir='cache'))
else:
    raise Exception("Bad value for tidal_bc_location: %s"%tidal_bc_location)

if tidal_bc_location=='decker':
    # unclear whether threemile should also be flipped.  mean flows typically Sac->SJ,
    # and the test period shows slightly negative means, so it's possible that it is
    # correct.
    # flipping this did improve a lot of phases, but stage at TSL is much worse, and
    # my best reading of the metadata is that it should not be flipped.
    model.add_bcs(dfm.NwisFlowBC(name='threemile',station=11337080,cache_dir='cache'))
    
# GSS: from compare_flows and lit, must be flipped.
model.add_bcs(dfm.NwisFlowBC(name='Georgiana',station=11447903,cache_dir='cache',
                                 filters=[dfm.Transform(lambda x: -x)] ))

model.add_bcs(dfm.NwisFlowBC(name='dcc',station=11336600,cache_dir='cache',
                             default=0.0,
                             filters=[dfm.FillGaps(large_gap_value=0.0),
                                      dfm.Transform(lambda x: -x)] ) )

# moving Sac flows upstream and removing tides.
sac=dfm.NwisFlowBC(name="SacramentoRiver",station=11447650,
                   pad=np.timedelta64(5,'D'),cache_dir='cache',
                   filters=[dfm.LowpassGodin(),
                            dfm.Lag(np.timedelta64(-2*3600,'s'))])
model.add_bcs(sac)


if 1:
    lisbon_bc=dfm.CdecFlowBC(name='lis',station="LIS",pad=np.timedelta64(5,'D'),
                             default=0.0)
    model.add_bcs(lisbon_bc)
else:
    log.warning("TEMPORARILY Disabling Lisbon flow due to CDEC issues")

# Ulatis inflow
# There are probably timezone issues here - they are coming in PST, but
# this code probably assumes UTC.
ulatis_ds=xr.open_dataset(os.path.join(here,"../bcs/ulatis/ulatis_hwy113.nc"))
ulatis_ds['flow']=ulatis_ds.flow_cfs*0.02832
ulatis_ds['flow'].attrs['units']='m3 s-1'
model.add_bcs(dfm.FlowBC(name='ULATIS',Q=ulatis_ds.flow))

# Campbell Lake
campbell_ds=xr.open_dataset(os.path.join(here,"../bcs/ulatis/campbell_lake.nc"))
campbell_ds['flow']=campbell_ds.flow_cfs*0.02832
campbell_ds['flow'].attrs['units']='m3 s-1'
model.add_bcs(dfm.FlowBC(name='CAMPBELL',Q=campbell_ds.flow))

    
if 0: # not including wind right now
    # Try file pass through for forcing data:
    windxy=model.read_tim('forcing-data/windxy.tim',columns=['wind_x','wind_y'])
    windxy['wind_xy']=('time','xy'),np.c_[ windxy['wind_x'].values, windxy['wind_y'].values]
    model.add_WindBC(wind=windxy['wind_xy'])
    

# Roughness 
# model.add_RoughnessBC(shapefile='forcing-data/manning_slick_sac.shp')
# model.add_RoughnessBC(shapefile='forcing-data/manning_n.shp')
if 1:
    import csc_dfm_decker_roughsearch as rs
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
    xyn=rs.settings_to_roughness_xyz(model,settings)
    # Turn that into a DataArray
    da=xr.DataArray( xyn[:,2],dims=['location'],name='n' )
    da=da.assign_coords(x=xr.DataArray(xyn[:,0],dims='location'),
                        y=xr.DataArray(xyn[:,1],dims='location'))
    da.attrs['long_name']='Manning n'
    rough_bc=dfm.RoughnessBC(data_array=da)
    model.add_bcs(rough_bc)

# Culvert at CCS
# rough accounting of restoration
if ccs_pre_restoration:
    # name => id
    # polylinefile is filled in from shapefile and name

    # these structure parameters get an amplitude ratio at CCS of 1.241
    # the shape is not great, and suggests that as the stage gets close to 1.7
    # or so, wetted area expands rapidly. That is probably a shortcoming of the
    # current bathymetry.
    model.add_Structure(name='ccs_breach',
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
        model.add_Structure(name=levee_name,
                            type='gate',
                            door_height=15,
                            lower_edge_level=3.5,
                            opening_width=0.0,
                            sill_level=3.5,
                            horizontal_opening_direction = 'symmetric')


        
# -- Extract locations for sections and monitor points
mon_sections=model.match_gazetteer(monitor=1,geom_type='LineString')
mon_points  =model.match_gazetteer(geom_type='Point')
model.add_monitor_sections(mon_sections)
model.add_monitor_points(mon_points)

if 1:
    # experimental saving of BC data to html plots
    for bc in model.bcs:
        bc.write_bokeh(path=model.run_dir)

##

import gen_polygons
gen_polygons.gen_polygons(model.run_dir)

## 
# if not invoked directly, just set up the model and let
# the importer decide what to do with model.

if __name__=='__main__':
    import shutil
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

