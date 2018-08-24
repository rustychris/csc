import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

##
X_path="/home/rusty/mirrors/ucd-X"
src_dir=os.path.join(X_path,"Arc_Hydro/CSC_Project/Common_Source_Data/N Delta hydro data/SCWA/ARabidoux")

cceh_path=os.path.join(src_dir,'CCEH.txt')
ccs1_path=os.path.join(src_dir,'CCS1.txt')
ccs_path="HecDssExcel1634846602065331917_CCS_wsel.csv"
ccs_orig_path=os.path.join(src_dir,"CCS.txt")
lshb_path=os.path.join(src_dir,'LSHB.txt')
dop_path =os.path.join(src_dir,'DOP.txt')

##
opts=dict(skiprows=2,sep="\t",parse_dates=['Date (PST)'],
          infer_datetime_format=True)

cceh=pd.read_csv(cceh_path,**opts).rename(columns={'Date (PST)':'Time'})
ccs1=pd.read_csv(ccs1_path,**opts).rename(columns={'Date (PST)':'Time'})
ccs_orig=pd.read_csv(ccs_orig_path,**opts).rename(columns={'Date (PST)':'Time'})
lshb=pd.read_csv(lshb_path,**opts).rename(columns={'Date (PST)':'Time'})
dop=pd.read_csv(dop_path,**opts).rename(columns={'Date (PST)':'Time'})
ccs=pd.read_csv(ccs_path,parse_dates=['Time'],infer_datetime_format=True)

##

for df in cceh,ccs1,ccs_orig,lshb,dop:
    bad=df['WSEL (NAVD88)']<-5
    df.loc[bad,'WSEL (NAVD88)']=np.nan
    df['WSEL_m']=df['WSEL (NAVD88)']*0.3048

##

plt.figure(1).clf()
fig,ax=plt.subplots(num=1)

ax.plot(cceh.Time,cceh['WSEL_m'],label='CCEH')
ax.plot(ccs1.Time,ccs1['WSEL_m'],label='CCS1')
ax.plot(lshb.Time,lshb['WSEL_m'],label='LSHB')

#ax.plot(ccs_orig.Time,ccs_orig['WSEL_m'],label='CCS orig.')
ax.plot(ccs_orig.Time + np.timedelta64(1100,'s'),
        ccs_orig['WSEL_m']-0.07,label='CCS orig lag')
ax.plot(ccs.Time,ccs['Stage'],label='CCS')

ax.legend()
ax.axis((734988.6442548756, 734988.9485486166, 0.6038336064456505, 1.430241285739577))

##

# At the very least, write LSHB and DOP out in formats ready for hydro_cal_plotter.
for df,name in [ (lshb,'LSHB'),
                 (dop,'DOP') ]:
    df=df.rename(columns={'WSEL_m':'Stage',
                          'Flow Discharge (cfs)':'Flow',
                          'Mean Discharge (cfs)':'Flow'})
    df.loc[:,['Time','Stage']].to_csv('%s_stage_m.csv'%name,index=False)
    df.loc[:,['Time','Flow']].to_csv('%s_flow_cfs.csv'%name,index=False)

