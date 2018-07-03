"""
Copy and adjust data in preparation for run_cal_plotter.py

This should not be needed routinely, but is here to document
the provenance of the calibratino data in the repository.
"""
from __future__ import print_function

import glob
import os
import pandas as pd
import shutil

##

csc_stations_dir="/home/rusty/mirrors/Ed/Sync/UCD/Projects/CDFW_Arc/dflowfm/27-comp_bedlevtyp_2_3_4/stations"

def any_x_in_y(x,y):
    """ true if any of the elements of list x are found in y"""
    for an_x in x:
        if an_x in y:
            return True
    return False

for fn in glob.glob(os.path.join( csc_stations_dir,"wsel","*.csv")):
    base_fn=os.path.basename(fn)
    print(fn)
    if base_fn in ['GES_STAGE_april2014.csv']:
        print("  Copying without modification")
        shutil.copyfile(fn,base_fn)
    elif base_fn.startswith('HecDssExcel'):
        if any_x_in_y( ('UL1','DOP','HAAS'),base_fn):  # no correction
            # unsure of status of UL1
            print("  CWS station with no time offset")
            shutil.copyfile(fn,base_fn)
        elif any_x_in_y( ('HS1','LN2','SG1','CCS'), base_fn):
            # Not positive about CCS, though.  Will apply 1 hour shift anyway
            print("  CWS station with assumed 1 hour offset")
            df=pd.read_csv(fn,parse_dates=['Time'],infer_datetime_format=True)
            df.Time -= np.timedelta64(3600,'s')
            # Transition to ISO date formats
            df.to_csv(base_fn,index=False,date_format="%Y-%m-%d %H:%M")
        else:
            raise Exception("No match for %s"%base_fn)




