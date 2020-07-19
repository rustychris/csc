import six
import os
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib import ticker

from stompy import utils
import matplotlib.pyplot as plt


from stompy.model import data_comparison
import stompy.model.hydro_model as hm
from stompy.model.delft import dflow_model

import pandas as pd

##

# model=dflow_model.DFlowModel.load("runs/20180807_grid98_15")
six.moves.reload_module(dflow_model)

model=dflow_model.DFlowModel.load('runs/v03/v03_20190415/flowfm.mdu')

##

fig_dir=os.path.join(os.path.dirname(model.his_output()),"figs-20200718")
os.path.exists(fig_dir) or os.mkdir(fig_dir)

## 
# CSV datasources:
base_dir="."
csv_dir=os.path.join(base_dir,"calibration_data")

##
def add_flow(name,csv_fn,**extra):
    fn=os.path.join(csv_dir,csv_fn)
    df=pd.read_csv(fn,parse_dates=['Time'])
    obs=xr.Dataset.from_dataframe(df)
    # standardize and add metadata
    obs['time']=('time',),df.Time.values
    obs.time.attrs['timezone']='PST'
    obs['flow']=('time',),obs.Flow*0.028316847
    obs.flow.attrs['units']='m3 s-1'
    obs_da=obs.flow

    mod=model.extract_section(name=name)

    model_da=mod.cross_section_discharge
    model_da.name='flow'
    model_da=model_da.isel(time=model_da.time-model_da.time.values[0]>np.timedelta64(24,'h'))

    model_da.attrs['label']="Model"
    obs_da.attrs['label']="Obs."

    model_da=model_da.assign_coords(label="Model")
    obs_da=obs_da.assign_coords(label="Obs.")

    plot_def=dict(sources=[obs_da,model_da],
                  station_name=name)
    plot_def.update(extra)
    plot_defs.append(plot_def)

def add_stage(name,csv_fn,**extra):
    fn=os.path.join(csv_dir,csv_fn)
    df=pd.read_csv(fn,parse_dates=['Time'])
    obs=xr.Dataset.from_dataframe(df)
    # standardize and add metadata
    obs['time']=('time',),df.Time.values
    obs.time.attrs['timezone']='PST'
    obs['water_level']=('time',),obs.Stage
    obs.water_level.attrs['units']='m'
    
    # some source files have missing records
    invalid=utils.isnat(obs.time.values)
    obs=obs.isel(time=~invalid)
    
    obs_da=obs.water_level

    mod=model.extract_station(name=name)

    model_da=mod.waterlevel
    model_da.name='water_level'
    model_da=model_da.isel(time=model_da.time-model_da.time.values[0]>np.timedelta64(24,'h'))

    model_da.attrs['label']="Model"
    obs_da.attrs['label']="Obs."

    model_da=model_da.assign_coords(label="Model")
    obs_da=obs_da.assign_coords(label="Obs.")

    plot_def=dict(sources=[obs_da,model_da],
                  station_name=name)
    plot_def.update(extra)
    plot_defs.append(plot_def)
        
##

plot_defs=[]
add_flow("DOP","DOP_flow_cfs.csv")
add_stage("CCS",'CCS_orig_stage_m.csv')
add_stage("LN2",'HecDssExcel6060592544570966047_LN2_wsel.csv',
          zoom_period=[735336.386, 735347.])

# And some on-demand downloaded data:
# Or chained..
model_period=[model.run_start,model.run_stop]

# What USGS and WDL stations are worth looking at?

# def add_usgs_stage(label,station,parameter=65):
#     fn=os.path.join(csv_dir,csv_fn)
#     df=pd.read_csv(fn,parse_dates=['Time'])
#     obs=xr.Dataset.from_dataframe(df)
#     # standardize and add metadata
#     obs['time']=('time',),df.Time.values
#     obs.time.attrs['timezone']='PST'
#     obs['water_level']=('time',),obs.Stage
#     obs.water_level.attrs['units']='m'
#     
#     # some source files have missing records
#     invalid=utils.isnat(obs.time.values)
#     obs=obs.isel(time=~invalid)
#     
#     obs_da=obs.water_level
# 
#     mod=model.extract_station(name=name)
# 
#     model_da=mod.waterlevel
#     model_da.name='water_level'
#     model_da=model_da.isel(time=model_da.time-model_da.time.values[0]>np.timedelta64(24,'h'))
# 
#     model_da.attrs['label']="Model"
#     obs_da.attrs['label']="Obs."
# 
#     model_da=model_da.assign_coords(label="Model")
#     obs_da=obs_da.assign_coords(label="Obs.")
# 
#     plot_def=dict(sources=[obs_da,model_da],
#                   station_name=name)
#     plot_def.update(extra)
#     plot_defs.append(plot_def)
    
##
six.moves.reload_module(hm)
six.moves.reload_module(data_comparison)
# 
rv_stage=hm.NwisStageBC(name='SRV',station=11455420)
# Already have a framework in place to handle data like NWIS stage,
# and it's tucked into the BC class hierarchy.
# What does extracting a model DataArray look like then?

# When loading a model from a run directory, can put together
# a basic gazetteer from BC and output selections

HERE -- move this logic into data_comparison.assemble_comparison_data.

bc_data=rv_stage
if isinstance( bc_data, hm.StageBC):
    print("okie -- will get stage")
    # In an ideal world, the model would have a gazetteer.
    # feat=model.match_gazetteer( name=rv_stage.name )
    # if feat is None:
    #     print("Hmm - will have to be more clever about figuring out location")
    # else:
    model_ds=model.extract_station(name=rv_stage.name)
    model_da=model_ds['waterlevel']

## 

plot_def=dict(sources=[rv_stage, model_da],
              station_name=rv_stage.name)
plot_defs=[plot_def]

# add_usgs_stage('Decker', 11455478, parameter=65)
# add_usgs_velocity('Decker', 11455478, parameter=72254) # sensor velocity
# add_usgs_stage('Rio Vista', 11455420)
# add_usgs_discharge('Rio Vista', 11455420, parameter=60)
# add_usgs_velocity('Rio Vista', 11455420, parameter=72255) # mean water velocity
# 
# add_usgs_discharge('Cache above Ryer',11455385)
# add_usgs_stage('Cache above Ryer',11455385)
# 
# add_usgs_discharge('Cache at Liberty',11455315)
# add_usgs_stage('Cache at Liberty',11455315)


defaults={}
# defaults['zoom_period']=[735336.386, 735347.]
# zoom_period=[735326.6,735342.00] )

for plot_def in plot_defs:
    settings=dict(defaults)# copy
    settings.update(plot_def)
    srcs=settings['sources']
    fig=data_comparison.calibration_figure_3panel(srcs,trim_time=True,
                                                  styles=[dict(color='k'),
                                                          dict(color='g')])
    if fig is None:
        continue
    param=srcs[0].name
    
    if param=='flow':
        fig.axes[0].set_ylabel('Flow (m$^3$s$^{-1}$)')
        fig.axes[1].set_ylabel('Residual flow (m$^3$s$^{-1}$)')
        fig.axes[0].set_title("Flow Calibration: %s"%plot_def['station_name'])
    elif param=='water_level':
        fig.axes[0].set_ylabel('Stage (m)')
        fig.axes[1].set_ylabel('Lowpass stage (m)')
        fig.axes[0].set_title("Stage Calibration: %s"%plot_def['station_name'])
        
    fig.axes[0].axis(xmin=settings['zoom_period'][0],
                     xmax=settings['zoom_period'][1])

    # this really should be automatic, but RRuleLocator does some weird stuff
    # with months.
    fig.axes[0].xaxis.set_major_locator(ticker.MultipleLocator(3))
    fig.axes[1].xaxis.set_major_locator(ticker.MultipleLocator(7))
    fig.subplots_adjust(wspace=0.35)

    img_fn=os.path.join(fig_dir,"%s-%s.png"%(plot_def['station_name'],
                                             param))
    fig.savefig(img_fn,dpi=200)

##     


# Older code

# get DOP data
fn=os.path.join(csv_dir,'DOP_flow_cfs.csv')
df=pd.read_csv(fn,parse_dates=['Time'])
dop_obs=xr.Dataset.from_dataframe(df)
# standardize and add metadata
dop_obs['time']=('time',),df.Time.values
dop_obs.time.attrs['timezone']='PST'
dop_obs['flow']=('time',),dop.Flow*0.028316847
dop_obs.flow.attrs['units']='m3 s-1'

##




dop_mod=model.extract_section(name='DOP')

##

model_da=dop_mod.cross_section_discharge
model_da.name='flow'
obs_da=dop_obs.flow


model_da.attrs['label']="Model"
obs_da.attrs['label']="Obs."

model_da=model_da.assign_coords(label="Model")
obs_da=obs_da.assign_coords(label="Obs.")
##
six.moves.reload_module(data_comparison)

all_sources=[model_da,obs_da]

combined=data_comparison.combine_sources(all_sources)
data_comparison.calibration_figure_3panel(all_sources,combined,num=1)


# plt.figure(1).clf()
# 
# plt.plot(model_da.time,model_da,label='model')
# plt.plot(obs_da.time,obs_da,label='obs')


