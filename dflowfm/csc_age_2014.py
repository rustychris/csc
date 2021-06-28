import csc_dfm_oper
import numpy as np
import os, shutil

target_date=np.datetime64("2014-10-01")
spinup=np.timedelta64(30,'D')
post  =np.timedelta64(1,'D')

model=csc_dfm_oper.CscDeckerModel(run_start=target_date-spinup,
                                  run_stop=target_date +post,
                                  run_dir="runs/age2014_v00")

# Not too much output for starters
model.mdu['Output','MapInterval']=86400//2

model.write()
shutil.copyfile(__file__,os.path.join(model.run_dir,"script.py"))
model.partition()
model.run_simulation()
