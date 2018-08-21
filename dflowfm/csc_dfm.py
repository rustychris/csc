#!/usr/bin/env python
"""
Script to configure and generate mdu ready for running dflow fm
for Cache Slough Complex
"""

import os
import shutil
import numpy as np
import logging
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
model.num_procs=4
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
model.set_run_dir("runs/20180807_grid98_02", mode='clean')

model.run_start=np.datetime64('2014-04-01')
model.run_stop=np.datetime64('2014-05-01')

# model.set_grid("CacheSloughComplex_v95_bathy01_net.nc")
# model.set_grid("CacheSloughComplex_v97_bathy_net.nc")
# model.set_grid("CacheSloughComplex_v98_bathy_net.nc")
model.set_grid("CacheSloughComplex_v98_bathy_sparse_net.nc")

model.load_mdu('template.mdu')
model.mdu['output','MapInterval']=7200
model.mdu['physics','UnifFrictCoef']= 0.023

model.set_cache_dir('cache')

model.add_gazetteer('gis/model-features.shp')

# Default to no dredging for flow and discharge BCs.
# 2018-08-19: why did I disable dredging?
dflow_model.SourceSinkBC.dredge_depth=-1
dflow_model.FlowBC.dredge_depth=-1

# original data was in seconds -- so use that explicitly
barker=model.read_tim('forcing-data/Barker_Pumping_Plant.tim',columns=['Q','s','T'],
                      time_unit='S')
model.add_FlowBC(name='Barker_Pumping_Plant',Q=barker['Q'])

rio_vista=model.read_bc('forcing-data/WaterLevel.bc')['SRV_0001']
model.add_StageBC(name='SRV',z=rio_vista['waterlevelbnd'])

Qshared=model.read_bc('forcing-data/Discharge.bc')
model.add_FlowBC(name='Georgiana',Q=Qshared['Georgiana_0001']['dischargebnd'])
model.add_FlowBC(name='DXC',Q=Qshared['DXC_0001']['dischargebnd'])

if 0:
    # original flows, maybe shifted in time too much?
    model.add_FlowBC(name="SacramentoRiver",Q=Qshared['SacramentoRiver_0001']['dischargebnd'])
elif 0:
    ds_fpx=usgs_nwis.nwis_dataset(station=11447650,
                                  start_date=model.run_start - np.timedelta64(5,'D'),
                                  end_date=model.run_stop + np.timedelta64(5,'D'),
                                  products=[60,65],cache_dir="cache")
    ds_fpx.time.values -= np.timedelta64(8,'h') # convert UTC to PST.
    # original data had flows shifted 14.5h, try just 2h? it's possible that they really should
    # be lagged, not shifted back in time, since the signal is mostly tidal and we're talking
    # about propagation of the tides.
    ds_fpx.time.values -= np.timedelta64(int(2.*3600),'s')
    ds_fpx.stream_flow_mean_daily.values *= 0.9*0.028316847 # cfs => m3/s, and scale down a bit
    model.add_FlowBC(name="SacramentoRiver",Q=ds_fpx.stream_flow_mean_daily)
elif 0:  # verona flows+american at fair oaks, no lag.
    dss=[]
    for stn in [11425500,11446500]:
        ds=usgs_nwis.nwis_dataset(station=stn,
                                  start_date=model.run_start - np.timedelta64(5,'D'),
                                  end_date=model.run_stop + np.timedelta64(5,'D'),
                                  products=[60,65],cache_dir="cache")
        ds.time.values -= np.timedelta64(8,'h') # convert UTC to PST.
        ds.stream_flow_mean_daily.values *= 0.028316847 # cfs => m3/s, and scale down a bit
        dss.append(ds)
    flow_sum=dss[0].stream_flow_mean_daily + dss[1].stream_flow_mean_daily

    # that's 15 minute data.
    lp_values=filters.lowpass_godin(flow_sum.values,
                                    utils.to_dnum(flow_sum.time))
    flow_sum.values[:]=lp_values
    model.add_FlowBC(name="SacramentoRiver",Q=flow_sum)
elif 1:  # freeport flows, lowpass
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
if 0:
    model.add_RoughnessBC(shapefile='forcing-data/manning_n.shp')
elif 0:
    rough=wkb2shp.shp2geom('forcing-data/manning_n.shp')
    rough['n']=rough['n'].clip(0.027,np.inf)
    wkb2shp.wkb2shp('forcing-data/manning_n_027.shp',rough['geom'],fields={'n':rough['n']},
                    overwrite=True)
    model.add_RoughnessBC(shapefile='forcing-data/manning_n_027.shp')

# TODO: shift these to come in from GIS
model.add_extra_file('ND_stations.xyn')
model.add_extra_file('FlowFMcrs.pli')


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

##

# This is failing on mpi when using WAQ output.

# There are several hundred messages about ** INFO   : Removed link     754, because of tiny angles at endpoints.
# not sure which domain that is coming from.
# is it not properly loading metis?  on the command line it looks okay, no need even to set
# LD_LIBRARY_PATH
# Also weird that I don't get those messages when runnig 53925-opt from the command line.
# is it related to forcing?
# pare it down: remove obs points, cross-sections, all entries from flowfm.ext => no help
# drop waq output => runs??
# Try that from the start with the script, but omit WAQ output.

# Now running on 4 cores.
# BLT 4 run had no tides.  revert to BLT and make sure I can get back to the previous state.
# is it possible that tim files are always in minutes, regardless of Tunit?

# I'd like to replicate the old run, but it's going to be annoying.
# switch the code to assume tim is always minutes, but then temporarily in csc_dfm.py
# muck with the tim data to get back to the old behavior.

##

# Problem with forcing - seems all of the timestamps are off??
# 
