# Synthesize Ulatis flow from some rating information
# and a long time series of stage.
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
ulatis_stage=os.path.join(ulatis_data,"UlatisCreek_2010-Present.xlsx")

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

# Rating curves:
# flow/stage data is from Ulatis Creek at Hwy 113.

ratings=["SC13-2003 WY.CSV",
         "SC13-2004 WY.CSV",
         "SC13-2005 WY.CSV"]

rating_fns=[os.path.join(ulatis_data,"Ulatis_Ratings",fn) for fn in ratings]

dfs=[pd.read_csv(fn,skiprows=3,names=['date','stage','cfs'],parse_dates=['date'])
     for fn in rating_fns]
df=pd.concat(dfs)


##

plt.figure(1).clf()
fig,axs=plt.subplots(2,1,sharex=True,num=1)

axs[0].plot(df.date, df.stage, label='stage')
axs[1].plot(df.date, df.cfs, label='cfs')

##


stage=df.stage.values
cfs  =df.cfs.values
valid=np.isfinite(stage)&np.isfinite(cfs)

coeffs=np.polyfit(stage[valid],cfs[valid],3)

ordered=stage.copy()
ordered.sort()

fit=np.polyval(coeffs,ordered)

best=None
def ulatis_stage_to_flow(x,coeffs=None):
    if coeffs is None:
        coeffs=best
    s=x.clip(coeffs[1],np.inf)
    return coeffs[0]*(s-coeffs[1])**coeffs[2]        

def cost(coeffs):
    pred=ulatis_stage_to_flow(stage[valid],coeffs)
    return np.mean( (pred-cfs[valid])**2 )

from scipy.optimize import fmin

best=fmin(cost,[80,0.2,1.7])

plt.figure(2).clf()
plt.plot( df.stage, df.cfs, 'g.')
plt.plot(ordered,fit,'k-',label='Fit')

plt.plot(ordered,ulatis_stage_to_flow(ordered),label="opt power")

plt.xlabel('Stage')
plt.ylabel('CFS')

##

# cubic fit is okay, but goes negative at low values.
# the power law is better, and non-negative.
# data max out around 8000.

# not sure if Alamo Creek flows need to be added in.
# though Easterly is part of that.
ulatis_hwy113=pd.read_excel(os.path.join(ulatis_data,"Ulatis @ Hwy 113.xlsx"),
                            skiprows=2)
ulatis_hwy113.rename({'Date (PST)':'time',
                      'Stage (ft)':'stage_ft'},
                     axis=1,
                     inplace=True)
                           
##

ulatis_hwy113['flow_cfs']=ulatis_stage_to_flow(ulatis_hwy113.stage_ft)

##

plt.figure(3).clf()
fig,axs=plt.subplots(2,1,sharex=True,num=3)

axs[0].plot(df.date, df.stage, label='stage')
axs[0].plot(ulatis_hwy113.time, ulatis_hwy113.stage_ft, label='pred.')

axs[1].plot(df.date, df.cfs, label='cfs')
axs[1].plot(ulatis_hwy113.time, ulatis_hwy113.flow_cfs, label='pred')

axs[0].set_ylabel('Stage (ft)')
axs[0].legend()
axs[1].set_ylabel('Flow (cfs)')
axs[1].legend()

# Weird -
#  the recent data shows high stage each year from April through October.
#  the winter period is punctuating by short flood flows with higher
#  stage and flow, but it looks like there is some sort of
#  impounding going on in the summer. 

# To really use this, then we need some way of knowing when flows are
# real, and when the check dams are in place.

# For starters, see how DCD net and gross inflow compare.

##
# node 320 is from the junction of Ulatis and Cache Ck.
dcd_node=pd.read_excel("../dcd/node320.xls",
                       header=0,index_col=1,skiprows=[0,1,3,4,5,6])
print(dcd_node.head())
sel=dcd_node.index.values>np.datetime64("2000-01-01")

dcd_node_sel=dcd_node.iloc[sel,:]

##



plt.figure(4).clf()
fig,axs=plt.subplots(2,1,sharex=True,num=4)

axs[0].plot(df.date, df.stage, label='stage')
axs[0].plot(ulatis_hwy113.time, ulatis_hwy113.stage_ft, label='pred.')

axs[1].plot(df.date, df.cfs, label='cfs')
axs[1].plot(ulatis_hwy113.time, ulatis_hwy113.flow_cfs, label='pred')
for field in ['DIV-FLOW','DRAIN-FLOW','SEEP-FLOW']:
    axs[1].plot(dcd_node_sel.index.values,dcd_node_sel[field],label=field)

axs[0].set_ylabel('Stage (ft)')
axs[0].legend()
axs[1].set_ylabel('Flow (cfs)')
axs[1].legend()
axs[0].axis( (731914.0764334609, 732117.9989572018, -0.680919476, 14.531472356) )
axs[1].axis( ymin=-194.05608590185568, ymax=2160.107367229529)

fig.savefig('scwa-vs-dcd.png')
