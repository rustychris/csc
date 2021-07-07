import stompy.model.delft.waq_scenario as dwaq
import stompy.model.delft.dflow_model as dfm
import stompy.model.hydro_model as hm

import csc_dfm_oper
import barker_data
import numpy as np
import os, shutil
import six

six.moves.reload_module(hm)
six.moves.reload_module(dfm)
six.moves.reload_module(dwaq)
six.moves.reload_module(barker_data)
six.moves.reload_module(csc_dfm_oper)
    
target_date=np.datetime64("2014-10-01")
spinup=np.timedelta64(150,'D')
post  =np.timedelta64(1,'D')

model=csc_dfm_oper.CscDeckerModel(run_start=target_date-spinup,
                                  run_stop=target_date +post,
                                  dcd=True,
                                  run_dir="runs/age2014_v04")

# Not too much output for starters
# v00: wedged on a bum node
# v01: good, but short.
# v02: longer.  75 day spinup.
# v03: Considerably longer: 150 day spinup, in case DWS is lagging
# v04: Including DCD
if model.ref_date is None:
    model.ref_date=model.run_start

map_start=((target_date-np.timedelta64(1,'D')) - model.ref_date) / np.timedelta64(1,'s')
map_stop =(model.run_stop-model.ref_date)/np.timedelta64(1,'s')
model.mdu['Output','MapInterval']="900 %.0f %.0f"%(map_start,map_stop)

model.write()
shutil.copyfile(__file__,os.path.join(model.run_dir,"script.py"))
model.partition()
model.run_simulation()
