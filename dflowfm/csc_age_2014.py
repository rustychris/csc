import csc_dfm_oper
import numpy as np
import os, shutil

target_date=np.datetime64("2014-10-01")
spinup=np.timedelta64(150,'D')
post  =np.timedelta64(1,'D')

model=csc_dfm_oper.CscDeckerModel(run_start=target_date-spinup,
                                  run_stop=target_date +post,
                                  run_dir="runs/age2014_v03")

# Not too much output for starters
# v00: wedged on a bum node
# v01: good, but short.
# v02: longer.  75 day spinup.
# v03: Considerably longer: 150 day spinup, in case DWS is lagging
if model.ref_date is None:
    model.ref_date=model.run_start

map_start=((target_date-np.timedelta64(1,'D')) - model.ref_date) / np.timedelta64(1,'s')
map_stop =(model.run_stop-model.ref_date)/np.timedelta64(1,'s')
model.mdu['Output','MapInterval']="900 %.0f %.0f"%(map_start,map_stop)

model.write()
shutil.copyfile(__file__,os.path.join(model.run_dir,"script.py"))
model.partition()
model.run_simulation()
