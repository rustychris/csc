# Look at station flows and stage on a common axis trying
# to understand what drives bad phase in DWS.

import pandas as pd
import matplotlib.pyplot as plt

## 
# First, the data:
Qdws=pd.read_csv('calibration_data/DWS-flow.csv',parse_dates=['Time'])
Qryi=pd.read_csv('calibration_data/RYI-flow.csv',parse_dates=['Time'])
Qhwb=pd.read_csv('calibration_data/HWB-flow.csv',parse_dates=['Time'])

##

##
import stompy.model.delft.dflow_model as dfm
six.moves.reload_module(dfm)

model=dfm.DFlowModel.load('runs/grid100_13/flowfm.mdu')
his=xr.open_dataset(model.his_output())

##


def get_flow(name):
    name=name.encode()
    xs_idx=list(his.cross_section_name.values).index(name)
    return his.cross_section_discharge.isel(cross_section=xs_idx)
modQdws=get_flow('DWS')
modQryi=get_flow('RYI')
modQhwb=get_flow('HWB')

##
fig=plt.figure(1)
fig.clf()

fig,axs=plt.subplots(2,1,sharex=True,num=1)

axQ=axs[0]

def add_Q(label,df):
    ft3_to_m3=0.028316847
    tsel=(df.Time.values>=model.run_start) & (df.Time.values<=model.run_stop)
    axQ.plot( df.Time.values[tsel], ft3_to_m3*df.Flow.values[tsel], label=label)
    
add_Q('DWS',Qdws)
add_Q('RYI',Qryi)
add_Q('HWB',Qhwb)

axQ.plot( modQdws.time, modQdws.values, label='mod DWS')
axQ.plot( modQryi.time, modQryi.values, label='mod RYI')
axQ.plot( modQhwb.time, modQhwb.values, label='mod HWB')

axQ.legend()


##

# Confirming signs around the DCC

Qsdc=pd.read_csv('calibration_data/SDC-flow.csv',parse_dates=['Time'])
Qges=pd.read_csv('calibration_data/GES-flow.csv',parse_dates=['Time'])
Qgss=pd.read_csv('calibration_data/GSS-flow.csv',parse_dates=['Time'])
Qdlc=pd.read_csv('calibration_data/DCC-flow.csv',parse_dates=['Time'])

##

zoom=ax.axis()
if zoom[0]<70000:
    zoom=(736162.8298482717, 736166.1012803828, -20455.38228094729, 29577.86075025231) 

start=zoom[0]
stop=zoom[1]

fig=plt.figure(1)
fig.clf()
ax=fig.gca()

sel_sdc=utils.within(utils.to_dnum(Qsdc.Time.values),[start,stop])
sel_ges=utils.within(utils.to_dnum(Qges.Time.values),[start,stop])
sel_gss=utils.within(utils.to_dnum(Qgss.Time.values),[start,stop])
sel_dlc=utils.within(utils.to_dnum(Qdlc.Time.values),[start,stop])


ax.plot(Qsdc.Time.values[sel_sdc],
        Qsdc.Flow.values[sel_sdc],label='SDC')
ax.plot(Qges.Time.values[sel_ges],
        -Qges.Flow.values[sel_ges],label='GES')
ax.plot(Qgss.Time.values[sel_gss],
        Qgss.Flow.values[sel_gss],label='GSS')
ax.plot(Qdlc.Time.values[sel_dlc],
        -Qdlc.Flow.values[sel_dlc],label='DCC')

# this is positive for flow into the control volume.
# to get the best balance, only SDC retains a positive sign.
# this makes sense for GES.
# it's arbitrary for DCC and GSS, and this suggests that both are
# positive *out* of the domain
flow_sum=Qsdc.Flow.values[sel_sdc] - Qges.Flow.values[sel_ges] - Qgss.Flow.values[sel_gss] - Qdlc.Flow.values[sel_dlc]
# mean is -80 to -700. 
print("Mean balance: %.1f cfs"%flow_sum.mean())
ax.plot(Qsdc.Time.values[sel_sdc],flow_sum,label='sum')

ax.legend()
ax.axis(zoom)

ax.axhline(0.0,color='k')
