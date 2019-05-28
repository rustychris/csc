import six
from stompy.model import data_comparison
from stompy.model.delft import dflow_model
import pandas as pd
##
six.moves.reload_module(dflow_model)
## 


model=dflow_model.DFlowModel.load("runs/20180807_grid98_15")

##

# CSV datasources:
base_dir="."
csv_dir=os.path.join(base_dir,"calibration_data")

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


