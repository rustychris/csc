"""
Read timeseries from SCWA spreadsheet, write out simple CSV.
"""

import os
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
## 
# Stage data:
# Arc_Hydro="/media/cws/Arc_Hydro"
Arc_Hydro="/home/rusty/mirrors/ucd-X/Arc_Hydro"

ulatis_data=os.path.join(Arc_Hydro,"CSC_Project/Common_Source_Data/N Delta hydro data/SCWA/Ulatis/Ulatis_Data")
#ulatis_stage=os.path.join(ulatis_data,"UlatisCreek_2010-Present.xlsx")

##

easterly_mgd=pd.read_excel(os.path.join( ulatis_data, "Vacaville WWTP Effluent Flow 0110-1017 tp.xlsx"),
                           header=None,
                           names=["date","_","flow_mgd"],
                           skiprows=10,
                           usecols=[0,1,2])
valid=np.isfinite(easterly_mgd['flow_mgd'].values)

## 
ds=xr.Dataset()
ds['time'] = ('time',), easterly_mgd['date'].iloc[valid]
ds['flow'] = ('time',), 0.043812636 * easterly_mgd['flow_mgd'].iloc[valid]
ds.flow.attrs['units']="m3 s-1"
ds.to_netcdf("easterly_flow.nc")

##

# This reads in a spreadsheet for which we are not sure of the distribution
# rights.  
ulatis_hwy113=pd.read_excel(os.path.join("scwa-nonpublic","N20.Stage-Flow Data.CLake, Ulatis, Hass.xlsx"),
                            skiprows=3,sheet_name=3,usecols=[0,3,4],
                            names=['time','stage_ft','flow_cfs'])
print(ulatis_hwy113.head())


campbell_lake=pd.read_excel(os.path.join("scwa-nonpublic","N20.Stage-Flow Data.CLake, Ulatis, Hass.xlsx"),
                            skiprows=3,sheet_name=3,usecols=[0,1,2],
                            names=['time','stage_ft','flow_cfs'])
print(campbell_lake.head())


## 
ulatis_ds=xr.Dataset.from_dataframe(ulatis_hwy113.set_index('time'))
ulatis_ds['flow_cfs'].attrs['description']="Flow from Hwy 113 stage rating curve, with nominal 5 cfs when check dams in place"
ulatis_ds.to_netcdf('ulatis_hwy113.nc')

campbell_ds=xr.Dataset.from_dataframe(campbell_lake.set_index('time'))
campbell_ds.to_netcdf('campbell_lake.nc')

