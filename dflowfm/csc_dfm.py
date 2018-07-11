#!/usr/bin/env python
"""
Script to configure and generate mdu ready for running dflow fm
for Cache Slough Complex
"""

import os
import numpy as np
import logging
import xarray as xr
import six

from stompy import utils
import stompy.model.delft.io as dio
from stompy.model.delft import dfm_grid
from stompy.grid import unstructured_grid
from stompy.spatial import wkb2shp

if __name__=='__main__':
    logging.basicConfig(level=logging.INFO)

log=logging.getLogger('csc_dfm')


## --------------------------------------------------

from stompy.model.delft import dflow_model
six.moves.reload_module(dflow_model)
six.moves.reload_module(dio)

model=dflow_model.DFlowModel()

# Having issues with 53925-opt, and 52184-dbg, both
# with MPI.  Runs single core okay.
model.dfm_bin_dir="/home/rusty/src/dfm/r52184-dbg/bin"
model.num_procs=4
model.z_datum='NAVD88'
model.projection='EPSG:26910'

# Parameters to control more specific aspects of the run
model.set_run_dir("runs/base20180710", mode='clean')
model.run_start=np.datetime64('2014-04-01')
model.run_stop=np.datetime64('2014-05-01')

model.set_grid("CacheSloughComplex_v95_net.nc")
model.load_mdu('template.mdu')

model.set_cache_dir('cache')

model.add_gazetteer('gis/model-features.shp')

# Default to no dredging for flow and discharge BCs.
dflow_model.SourceSinkBC.dredge_depth=None
dflow_model.FlowBC.dredge_depth=None

barker=model.read_tim('forcing-data/Barker_Pumping_Plant.tim',columns=['Q','s','T'])
model.add_SourceSinkBC(name='Barker_Pumping_Plant',Q=barker['Q'])

rio_vista=model.read_bc('forcing-data/WaterLevel.bc')['SRV_0001']
model.add_StageBC(name='SRV',z=rio_vista['waterlevelbnd'])

# Is it weird that Georgiana, DXC, American and Sac all have the same forcing?
Qshared=model.read_bc('forcing-data/Discharge.bc')

model.add_FlowBC(name='Georgiana',Q=Qshared['Georgiana_0001']['dischargebnd'])
model.add_FlowBC(name='DXC',Q=Qshared['DXC_0001']['dischargebnd'])
model.add_FlowBC(name="SacramentoRiver",Q=Qshared['SacramentoRiver_0001']['dischargebnd'])
model.add_FlowBC(name="AmericanRiver",Q=Qshared['AmericanRiver_0001']['dischargebnd'])

# Try file pass through for forcing data:
windxy=model.read_tim('forcing-data/windxy.tim',columns=['wind_x','wind_y'])
windxy['wind_xy']=('time','xy'),np.c_[ windxy['wind_x'].values, windxy['wind_y'].values]
model.add_WindBC(wind=windxy['wind_xy'])

# TODO: shift these to come in from GIS
model.add_extra_file('ND_stations.xyn')
model.add_extra_file('FlowFMcrs.pli')

if __name__=='__main__':
    model.write()
    model.partition()
    model.run_model()



##
