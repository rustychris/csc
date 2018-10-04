#!/usr/bin/env python
"""
Script to configure and generate mdu ready for running dflow fm
for Cache Slough Complex
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


## --------------------------------------------------

from stompy.model.delft import dflow_model
six.moves.reload_module(dfm_grid)
six.moves.reload_module(dflow_model)
six.moves.reload_module(dio)

model=dflow_model.DFlowModel()

# Having issues with 53925-opt, and 52184-dbg, both
# with MPI.  Looks like dwaq output is not compatible
# with ugrid+mpi.
model.dfm_bin_dir="/home/rusty/src/dfm/r53925-opt/bin"
model.num_procs=1
model.z_datum='NAVD88'
model.projection='EPSG:26910'

# Parameters to control more specific aspects of the run
# 20180807_grid97: download FPX flows and shift -2h to use as
# sac flow.
# 20180807_grid97_02: use verona, no lag.
# 20180807_grid97_03: use verona, no lag, back to low friction settings.
# 20180807_grid98_00: extended sac grid, combine verona and american river, lp
#                     for upstream flow.
# 20180807_grid98_01: back to uniform friction settings as in baseline.
#      failed early.
#  early results appeared to be no better, so dial friction down to 0.02
#  everywhere, try again.
# runs/20180807_grid98_02: fancy setting of depths, targetting conveyance area.
#    return friction to 0.023
# runs/20180807_grid98_03: split friction 0.02 on sac, 0.023 in CSC, turn off DXC
#   as the gates were closed.
# 04: add structure
# 05: more outputs
# 06: more outputs, fix sign of LSHB section
# 07: edge depths are connectivity-based, not mean, and CCS
#   is dredge out.
# 08: conn-based was too much - back to edge-means
# 09: same as 08 but make the structure deeper
# 10: structure sill up to 0.47, for a gap of 0.03m, and zero out
#    opening width.
# 11: limit CCS culvert flow to center edge
# 12: increase friction in the marsh, and raise the gate to 0.85--0.88m
# 13: decrease the opening height down to 0.02m ??
# 14: based on report, make that 0.04, and crank up friction even more.
#     note that the high friction area now covers a small chunk of the channel, too.
# 15: 14 was too frictional -- back off n to 0.12.
# 16: a short, single core DWAQ-enabled run.
# 17: month long DWAQ-enabled run, back to MPI.
model.set_run_dir("runs/20180807_grid98_17", mode='noclobber')

model.run_start=np.datetime64('2014-04-01')
model.run_stop=np.datetime64('2014-05-01')

# this switches to low-biased edge depths, still optimized, and
# uses a master DEM with CCS cut down a bit.
model.set_grid("CacheSloughComplex_v98_bathy2_sparse_net.nc")

model.load_mdu('template.mdu')
model.mdu['output','MapInterval']=1800 # 7200
model.mdu['output','WaqInterval']=1800
model.mdu['physics','UnifFrictCoef']= 0.023

model.set_cache_dir('cache')

model.add_gazetteer('gis/model-features.shp')
model.add_gazetteer('gis/point-features.shp')

# Default to no dredging for flow and discharge BCs.
# 2018-08-19: why did I disable dredging?
dflow_model.SourceSinkBC.dredge_depth=-1
dflow_model.FlowBC.dredge_depth=-1

# check_bspp.py has the code that converted original tim to csv.
# awkward reloading of that, but at least it's independent of the
# simulation period.
barker=xr.Dataset.from_dataframe(pd.read_csv('forcing-data/Barker_Pumping_Plant.csv',
                                             parse_dates=['time']))
barker=barker.set_coords('time')
model.add_FlowBC(name='Barker_Pumping_Plant',Q=barker['Q'])

rio_vista=model.read_bc('forcing-data/WaterLevel.bc')['SRV_0001']
model.add_StageBC(name='SRV',z=rio_vista['waterlevelbnd'])

Qshared=model.read_bc('forcing-data/Discharge.bc')
model.add_FlowBC(name='Georgiana',Q=Qshared['Georgiana_0001']['dischargebnd'])
# DXC was closed during this period.
# model.add_FlowBC(name='DXC',Q=Qshared['DXC_0001']['dischargebnd'])

if 1:  # freeport flows, lowpass
    ds_fpx=usgs_nwis.nwis_dataset(station=11447650,
                                  start_date=model.run_start - np.timedelta64(5,'D'),
                                  end_date=model.run_stop + np.timedelta64(5,'D'),
                                  products=[60,65],cache_dir="cache")
    ds_fpx.time.values -= np.timedelta64(8,'h') # convert UTC to PST.
    # original data had flows shifted 14.5h, try just 2h? it's possible that they really should
    # be lagged, not shifted back in time, since the signal is mostly tidal and we're talking
    # about propagation of the tides.
    ds_fpx.time.values -= np.timedelta64(int(2.*3600),'s')
    ds_fpx.stream_flow_mean_daily.values *= 0.028316847 # cfs => m3/s
    # that's 15 minute data.
    flow=ds_fpx.stream_flow_mean_daily
    flow.values[:]=filters.lowpass_godin(flow.values,
                                         utils.to_dnum(flow.time))
    model.add_FlowBC(name="SacramentoRiver",Q=flow)

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


##

mon_sections=model.match_gazetteer(monitor=1,geom_type='LineString')
mon_points  =model.match_gazetteer(geom_type='Point')
model.add_monitor_sections(mon_sections)
model.add_monitor_points(mon_points)

#
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

