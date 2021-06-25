import six
import os
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib import ticker

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from stompy import utils
from stompy.model import data_comparison
import stompy.model.hydro_model as hm
from stompy.model.delft import dflow_model

## 
# CSV datasources:
#base_dir="."
#csv_dir=os.path.join(base_dir,"calibration_data")

six.moves.reload_module(hm)
six.moves.reload_module(data_comparison)
six.moves.reload_module(dflow_model)

plot_defs=[]

model=dflow_model.DFlowModel.load('runs/sb_rv4_20190601')
model.utc_offset=np.timedelta64(-8,'h') # model runs in PST

# And some on-demand downloaded data:
# Or chained..
model_period=[model.run_start,model.run_stop]

# change because run didn't finish
# model_period[1] = np.datetime64('2019-08-08')

# Passing the model to the data source enables time zone correction

plot_defs=[
    # dict(sources=[hm.NwisScalarBC(name='SDI', station=11455478, scalar='salinity',
    #                               cache_dir='cache', model=model)],
    #      station_name='Sac R. at Decker Island'),
    # Decker Island no longer in model domain
    # dict(sources=[hm.NwisStageBC(name='SDI', station=11455478,
    #                              cache_dir='cache', model=model)],
    #      station_name='Sac R. at Decker Island'),
    dict(sources=[hm.NwisStageBC(name='SRV',station=11455420,
                                 cache_dir='cache',model=model)],
         station_name='Rio Vista'),

    dict(sources=[hm.NwisFlowBC(name='RioVista',station=11455420,
                                cache_dir='cache',model=model)],
         station_name='Rio Vista'),

    dict(sources=[hm.NwisFlowBC(name='RYI',station=11455385,
                                cache_dir='cache',model=model)],
         station_name='Cache above Ryer'),

    dict(sources=[hm.NwisStageBC(name='RYI',station=11455385,
                                 cache_dir='cache',model=model)],
         station_name='Cache above Ryer'),
    
    dict(sources=[hm.NwisStageBC(name='LIB',station=11455315,
                                 cache_dir='cache',model=model)],
         station_name='Cache at Liberty'),

    dict(sources=[hm.NwisFlowBC(name='LIB',station=11455315,
                                cache_dir='cache',model=model)],
         station_name='Cache at Liberty'),

    dict(sources=[hm.NwisFlowBC(name='HST',station=11455280,
                                cache_dir='cache',model=model)],
         station_name='Cache at Hastings'),
    
    dict(sources=[hm.NwisStageBC(name='HST',station=11455280,
                                 cache_dir='cache',model=model)],
         station_name='Cache at Hastings'),

    dict(sources=[hm.NwisStageBC(name='SG1',station=11455276,
                                 cache_dir='cache',model=model)],
         station_name='Shag, near Courtland'),
    dict(sources=[hm.NwisFlowBC(name='SG1',station=11455276,
                                cache_dir='cache',model=model)],
         station_name='Shag, near Courtland'),

    # These just have BGC parameters-- 
    # dict(sources=[hm.NwisScalarBC(name='DWSCTL',station=11455142, scalar='salinity',
    #                               cache_dir='cache',model=model)],
    #      station_name="DWS, Courtland"),
    # dict(sources=[hm.NwisFlowBC(name='DWSCTL',station=11455142,
    #                             cache_dir='cache',model=model)],
    #      station_name="DWS, Courtland"),
    
    dict(sources=[hm.NwisStageBC(name='DWS',station=11455335,
                                cache_dir='cache',model=model)],
         station_name="DWS"),
    dict(sources=[hm.NwisFlowBC(name='DWS',station=11455335,
                                cache_dir='cache',model=model)],
         station_name="DWS"),
]    
    
defaults=dict( zoom_period=[model.run_stop-np.timedelta64(5,'D'),
                            model.run_stop],
               models=[model] )

# also update because run didn't finish
# defaults=dict( zoom_period=[model_period[1]-np.timedelta64(5,'D'),
#                             model_period[1]],
#                models=[model] )


fig_dir=os.path.join(os.path.dirname(model.his_output()),"cal_figs")
os.path.exists(fig_dir) or os.mkdir(fig_dir)

force=True

for plot_def in plot_defs:
    settings=dict(defaults) # copy
    settings.update(plot_def)

    sources,combined=data_comparison.assemble_comparison_data(settings['models'],
                                                              settings['sources'])
    param=combined.name
    if len(sources)<2:
        print("%s has %d sources (maybe model missing output here for %s)"%(plot_def['station_name'],
                                                                            len(sources),
                                                                            param))
        continue

    img_fn=os.path.join(fig_dir,"%s-%s.png"%(plot_def['station_name'],
                                             param))
    if (not force) and os.path.exists(img_fn):
        print('%s exists'%img_fn)
        continue
    
    fig=data_comparison.calibration_figure_3panel(sources,trim_time=True,
                                                  styles=[dict(color='k'),
                                                          dict(color='g')])
    if fig is None:
        continue
    if param=='flow':
        fig.axes[0].set_ylabel('Flow (m$^3$s$^{-1}$)')
        fig.axes[1].set_ylabel('Residual flow (m$^3$s$^{-1}$)')
        fig.axes[0].set_title("Flow Calibration: %s"%plot_def['station_name'])
    elif param=='water_level':
        fig.axes[0].set_ylabel('Stage (m)')
        fig.axes[1].set_ylabel('Lowpass stage (m)')
        fig.axes[0].set_title("Stage Calibration: %s"%plot_def['station_name'])
    elif param=='salinity':
        fig.axes[0].set_ylabel('Salinity (ppt)')
        fig.axes[1].set_ylabel('Lowpass salinity (ppt)')
        fig.axes[0].set_title("Salinity Calibration: %s" % plot_def['station_name'])
        
    fig.axes[0].axis(xmin=settings['zoom_period'][0],
                     xmax=settings['zoom_period'][1])

    # this really should be automatic, but RRuleLocator does some weird stuff
    # with months.
    fig.axes[0].xaxis.set_major_locator(ticker.MultipleLocator(1))
    fig.axes[1].xaxis.set_major_locator(mdates.DayLocator(bymonthday=[1, 15]))
    fig.axes[1].tick_params(axis='x', labelrotation=-15)
    fig.subplots_adjust(wspace=0.35)
    fig.set_size_inches(10, 8)
    fig.savefig(img_fn,dpi=200, bbox_inches='tight')
    print(img_fn)

plt.show()
