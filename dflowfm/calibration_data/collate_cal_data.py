"""
Copy and adjust data in preparation for run_cal_plotter.py

This should not be needed routinely, but is here to document
the provenance of the calibratino data in the repository.
"""
from __future__ import print_function
import logging
import glob
import os
import pandas as pd
import shutil

cache_dir="../cache"


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


##

# USGS stations:
download_period=[np.datetime64("2014-03-01"),
                 np.datetime64("2014-06-01")]
from stompy.io.local import usgs_nwis


# USGS gauges with Flow and Stage:
for usgs_name,usgs_station in [ ("SRV","11455420"),  # Sac River at Rio Vista
                                ("FPX","11447650"),  # Sac River at Freeport
                                ("RYI", "11455350"), # Cache Slough at Ryer Island
                                ("HWB","11455165"),  # Miner Slough at HWY 84 Bridge
                                ("SSS","11447850"),  # Steamboat Slough Btw Sac R And Sutter Sl, aka Steamboat Slough nr Walnut
                                ("SUT","11447830"),  # Sutter Slough at Courtland
                                ("DWS","11455335"),  # Sacramento R Deep Water Ship Channel Nr Rio Vista
                                # no physical data until 2015-07:
                                # ("LIB","11455315"),  # Cache Slough A S Liberty Island Nr Rio Vista CA
]:
    ds=usgs_nwis.nwis_dataset(usgs_station,
                              download_period[0],download_period[1],
                              [60,65], # Discharge and Stage
                              days_per_request='M',cache_dir=cache_dir)
    # nwis_dataset() returns UTC data.  Convert to PST:
    ds['time'] = ds.time - np.timedelta64(8,'h')

    # Match the names up with existing csv files:
    df=ds.rename(
        {'time':'Time',
         'stream_flow_mean_daily':'Flow',
         'height_gage':'Stage'}
    ).to_dataframe()

    df.Stage.to_csv('%s-2014-04-stage.csv'%usgs_name,index=True,date_format="%Y-%m-%d %H:%M",header=True)
    df.Flow.to_csv('%s-2014-04-flow.csv'%usgs_name,index=True,date_format="%Y-%m-%d %H:%M",header=True)

##
logging.warning("No timezone adjustments applied down here!  May need to go UTC=>PST")


from stompy import utils
# Fetch a year of stage data for Yolo Bypass near Lisbon.
# the original file is put in cache_dir, from which
lis_fn="LIS-stage-WY2014.csv"
lis_orig_fn=os.path.join(cache_dir,'LIS-WY2014-STAGE_15-MINUTE_DATA_DATA.CSV')
url="http://wdl.water.ca.gov/waterdatalibrary/docs/Hydstra/docs/B91560/2014/STAGE_15-MINUTE_DATA_DATA.CSV"

# only fetch when necessary.
if not os.path.exists(lis_fn):
    if not os.path.exists(lis_orig_fn):
        utils.download_url(url,lis_orig_fn,log=logging)
    df=pd.read_csv(lis_orig_fn,skiprows=3,parse_dates=['Time'],infer_datetime_format=True,
                   names=["Time","Stage","Quality","Comment"])
    df.to_csv(lis_fn,columns=["Time","Stage"],date_format="%Y-%m-%d %H:%M",index=False)





