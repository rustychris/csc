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

import local_config
local_config.install()

model=dfm.DFlowModel()

model.num_procs=4
model.z_datum='NAVD88' # not really used right now.
model.projection='EPSG:26910'
# Forcing data is fetched as UTC, and adjusted according to this offset
model.utc_offset=np.timedelta64(-8,'h') # PST.  

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
#        _16: try nonlin2d, and revert to constant n=0.023 roughness
#        _17: for apples-to-apples, same as 16 but nonlin2d=0
#        _18: keep nonlin2d=0, and bring in updated grid with doubled shorelines for
#             Sac and DWS
#        _19: hmm - step back to conveyance2d=-1
#             amplitudes are down, and lags are crazy bad at fpt. can looker closer, but this
#             seems much worse.
#        _20: switch to bedlevtype=6
#        _21: and use adjusted node elevations
#        _22: simpler adjustment just to get nodes to reflect means
#        _23: bringing back the settings from roughsearch
#        _24: add extra section on Liberty Cut, fix LIS BC, and add channel
#             up to LIY in the DEM.
model.set_run_dir("runs/grid100_24", mode='askclobber')

model.run_start=np.datetime64('2017-12-01')
model.run_stop=np.datetime64('2017-12-10')

model.load_mdu('template.mdu')

src_grid='../grid/CacheSloughComplex_v111-edit19fix.nc'
dst_grid=os.path.basename(src_grid).replace('.nc','-bathy.nc')
bathy_fn="../bathy/merged_2m-20190122.tif"
if utils.is_stale(dst_grid,[src_grid,bathy_fn]):
    g=unstructured_grid.UnstructuredGrid.from_ugrid(src_grid)
    dem=field.GdalGrid(bathy_fn)
    if 0:
        node_depths=dem(g.nodes['x'])
        while np.any(np.isnan(node_depths)):
            missing=np.nonzero(np.isnan(node_depths))[0]
            print("Looping to fill in %d missing node depths"%len(missing))
            for n in missing:
                nbrs=g.node_to_nodes(n)
                node_depths[n]=np.nanmean(node_depths[nbrs])
    else:
        import dem_cell_node_bathy
        node_depths=dem_cell_node_bathy.dem_to_cell_node_bathy(dem,g)
    g.add_node_field('depth',node_depths,on_exists='overwrite')
    g.write_ugrid(dst_grid,overwrite=True)
else:
    g=unstructured_grid.UnstructuredGrid.from_ugrid(dst_grid)
model.set_grid(g)

# I think this doesn't matter with conveyance2d>0
# 6 is maybe better for getting good edges
model.mdu['geometry','BedlevType']=6
# ostensibly analytic 2D conveyance with 3.
# 2 is faster, and appears less twitchy while showing very similar
#  calibration results.
# 0 blocks some flow, 1 was a little twitchy.
model.mdu['geometry','Conveyance2D']=-1
# enabling this seems to cause a lot of oscillation.
# but it may be necessary to get good prism on the narrow channels.
model.mdu['geometry','Nonlin2D']=0

model.mdu['output','MapInterval']=7200
model.mdu['physics','UnifFrictCoef']= 0.023
# fail out when it goes unstable.
model.mdu['numerics','MinTimestepBreak']=0.05

# Make sure cache dir exists
os.path.exists(cache_dir) or os.makedirs(cache_dir)


# -- Boundary Conditions --

# Register linear and point feature shapefiles
model.add_gazetteer('gis/model-features.shp')
model.add_gazetteer('gis/point-features.shp')

# All sources and flow BCs will be dredge to this depth to make
# they remain wet and active.
dfm.SourceSinkBC.dredge_depth=-1
dfm.FlowBC.dredge_depth=-1

# check_bspp.py has the code that converted original tim to csv.
model.add_bcs(barker_data.BarkerPumpsBC(name='Barker_Pumping_Plant'))

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

model.add_bcs(nwis_bc.NwisFlowBC(name='dcc',station=11336600,cache_dir='cache',
                                 filters=[dfm.FillGaps(large_gap_value=0.0),
                                          dfm.Transform(lambda x: -x)] ) )

# moving Sac flows upstream and removing tides.
sac=nwis_bc.NwisFlowBC(name="SacramentoRiver",station=11447650,
                       pad=np.timedelta64(5,'D'),cache_dir='cache',
                       filters=[dfm.LowpassGodin(),
                                dfm.Lag(np.timedelta64(-2*3600,'s'))])
model.add_bcs(sac)

lisbon_bc=dfm.CdecFlowBC(name='lis',station="LIS",pad=np.timedelta64(5,'D'))
model.add_bcs(lisbon_bc)

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

# if not invoked directly, just set up the model and let
# the importer decide what to do with model.

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

