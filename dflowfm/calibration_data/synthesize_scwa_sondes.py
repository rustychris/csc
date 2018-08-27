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
# Also write this version of the CCS data, as the DSS one I think has an erroneous
# 1 hour shift.
for df,name in [ (lshb,'LSHB'),
                 (dop,'DOP'),
                 (ccs_orig,'CCS_orig') ]:
    df=df.rename(columns={'WSEL_m':'Stage',
                          'Flow Discharge (cfs)':'Flow',
                          'Mean Discharge (cfs)':'Flow'})

    df.loc[:,['Time','Stage']].to_csv('%s_stage_m.csv'%name,index=False)
    if 'Flow' in df:
        df.loc[:,['Time','Flow']].to_csv('%s_flow_cfs.csv'%name,index=False)

##

# And estimate longer period of CCEH based on LSHB

plt.figure(2).clf()
fig,ax=plt.subplots(num=2)

ax.plot(cceh.Time,cceh['WSEL_m'],label='CCEH')
ax.plot(lshb.Time,lshb['WSEL_m'],label='LSHB')

ax.legend()
# ax.axis((734988.6442548756, 734988.9485486166, 0.6038336064456505, 1.430241285739577))

##

# time dimensions are compatible, 15-minute
cceh_dn=utils.to_dnum(cceh.Time)
lshb_dn=utils.to_dnum(lshb.Time)
cceh_h=cceh['WSEL_m']
lshb_h=lshb['WSEL_m']

# more lag at lower water
# around LSHB=2.0, 2700s lag is good
# at 0.75, a lag of 4000 is okay
# sharp-ish break at 0.9m
# this is just hand-tuned, but looks decent enough
# for this application
lag_lshb_s=4000
lag_sec_per_m=0

lag=np.interp( lshb_h,
               [0.8,1.0],
               [4000,2700] )
lshb_h_res=np.interp(cceh_dn,
                     lshb_dn+lag/86400.,
                     lshb_h)

# Turn that scatter into a monotonic function

fig=plt.figure(3)
fig.clf()
ax=fig.add_subplot(2,1,1)
ax_scat=fig.add_subplot(2,2,3)

ax.plot_date(cceh_dn,cceh_h,'b-',label='CCEH')
ax.plot_date(lshb_dn,lshb_h,'g-',label='LSHB')
ax_scat.plot(lshb_h_res, cceh_h,'k.',alpha=0.3)

valid=np.isfinite(lshb_h_res+cceh_h)

# quartic fit not bad
coeffs=np.polyfit(lshb_h_res[valid], cceh_h[valid],4)
ax_scat.plot(lshb_h_res, np.polyval(coeffs,lshb_h_res),label='Quartic fit')

# Reconstruct!
def lshb_to_cceh(target_dn,lshb_h,lshb_dn):
    lag=np.interp( lshb_h,
                   [0.8,1.0],
                   [4000,2700] )
    lshb_h_res=np.interp(target_dn,
                         lshb_dn+lag/86400.,
                         lshb_h)
    cceh_est=np.polyval(coeffs,lshb_h_res)
    return cceh_est

cceh_est=lshb_to_cceh(lshb_dn,lshb_h,lshb_dn)

ax.plot(lshb_dn,cceh_est,'k--',label='Reconstructed')

ax_scat.set_xlabel('LSHB resampled')
ax_scat.set_ylabel('CCEH')

ax.legend()
ax_scat.legend()
# ax.axis( (734947, 735003, -0.7707635966445383, 3.0148056777744747))
ax.axis( (734960.91856, 734963.0595, 0.40, 2.1342))

fig.savefig('CCEH_synthesise.png')

##
df=pd.DataFrame()
df['Time']=utils.to_dt64(lshb_dn)
df['Stage']=cceh_est

df.to_csv('CCEH_syn_stage_m.csv',index=False)
