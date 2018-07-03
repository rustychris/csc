his_file='DFM_OUTPUT_FlowFM/FlowFM_his.nc'

ds=xr.open_dataset(his_file)

##

fig=plt.figure(10)
fig.clf()
ax=fig.add_subplot(1,1,1)

srv_idx=list(ds.station_id.values).index(b'SRV')
ges_idx=list(ds.station_id.values).index(b'GES')
fpx_idx=list(ds.station_id.values).index(b'FPX')

ax.plot(ds.time.values,
        ds.waterlevel.isel(stations=srv_idx))

# 15 minutes into the run (after IC)
# we're at 1.31m here, just as WaterLevel.bc dictates

# From NWIS, 2014-04-01 00:15 PST has gage height of ..
#   4.10 ft
# max following that is at 2014-04-01 04:30 or 04:45
#   at 6.42 ft

# High tide from WaterLevel.bc is ... between 16200
# and 17100.  This also shows the same symmetry (ABBA)
# as nwis, so good.  And agrees with the
# 4:30 to 4:45 timing.
# Interestingly, the 6.42ft from USGS would be 1.9568m,
# but WaterLevel.bc has 2.0051,
# which represents a 1.9in adjustment??

# So far: Waterlevel.bc and USGS NWIS agree on the time zone of
#   Rio Vista stage, and USGS reports that as PST.
#  WaterLevel.bc includes an adjustment of +1.9in to stage.

# NOAA reports observed high tide at 02:36 PST at Port Chicago
# or 3:36 PDT.

# End of February of 2014, before PDT, (Mar 9)
#  NOAA gives first high tide of Mar 2014 observed at
#  01:48 / 01:54
#  USGS that day gives high tide at 04:00 PST
# That means that, at least for Feb 28 conditions, tidal
# propagation from Port Chicago to Rio Vista takes 2h10m
#

# While that seems long, it may be slowed down by Delta flows
# Based on rough analysis of CASCADE output, they got
# 124 min propagation time from Port Chicago to Rio Vista.
# Excellent.  We have decent agreement then.
# this was not during significant flows.

# All of this suggests that the WaterLevel.bc for Rio Vista
# is in fact properly in PST.

##

from stompy import utils

lag=utils.find_lag_xr(ds.waterlevel.isel(stations=ges_idx),
                      ds.waterlevel.isel(stations=srv_idx))
print("Lag SRV=>GES (min): %.2f"% (lag/np.timedelta64(60,'s')))

lag=utils.find_lag_xr(ds.waterlevel.isel(stations=fpx_idx),
                      ds.waterlevel.isel(stations=ges_idx))
print("Lag GES=>FPX (min): %.2f"% (lag/np.timedelta64(60,'s')))


##
plt.figure(11).clf()

for idx in [srv_idx,ges_idx,fpx_idx]:
    plt.plot( ds.time.values,ds.waterlevel.isel(stations=idx),label=ds.station_id.isel(stations=idx).values)
plt.legend()

##

# Pasted in from HPC.sfei.org
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from stompy.spatial import proj_utils
from stompy.grid import unstructured_grid
import stompy.model.delft.waq_scenario as waq
from stompy.model.delft import hydro_utils
from stompy import utils

srv_xy=[615117,4224383]
pc_ll=[-(122+2.4/60.),(38+3.3/60.)]
pc_xy=proj_utils.mapper('WGS84','EPSG:26910')(pc_ll)

ges_xy=[629223,4233353]
fpx_xy=[630728,4257433]

hydro=waq.HydroFiles('/hpcvol1/cascade/WY2011/DFM_DELWAQ_sal+temp/sal+temp.hyd')

g=hydro.grid()

plt.figure()
g.plot_edges(lw=0.3,color='k')
plt.axis('equal')

t0=utils.to_dt64(hydro.time0)

t1=t0+np.timedelta64(170,'D')
tn=t1+np.timedelta64(10,'D')

srv=hydro_utils.extract_water_level(hydro,srv_xy,t1,tn).isel(station=0)
pc=hydro_utils.extract_water_level(hydro,pc_xy,t1,tn).isel(station=0)
ges=hydro_utils.extract_water_level(hydro,ges_xy,t1,tn).isel(station=0)
fpx=hydro_utils.extract_water_level(hydro,fpx_xy,t1,tn).isel(station=0)

plt.figure()

plt.plot(srv.time,srv.water_level,label="SRV")
plt.plot(pc.time,pc.water_level,label="PC")
plt.plot(ges.time,ges.water_level,label="GES")
plt.plot(ges.time,fpx.water_level,label="FPX")

lag= utils.find_lag_xr( srv.water_level,pc.water_level)
print("Lag PC=>SRV (min): %.2f"% (lag/np.timedelta64(60,'s')))

lag= utils.find_lag_xr( ges.water_level,srv.water_level)
print("Lag SRV=>GES (min): %.2f"% (lag/np.timedelta64(60,'s')))

lag= utils.find_lag_xr( fpx.water_level, ges.water_level)
print("Lag GES=>FPX (min): %.2f"% (lag/np.timedelta64(60,'s')))

##
import pandas as pd
# Look at two copies of GES stage:
ges_mwtract=pd.read_csv( ('/home/rusty/mirrors/ucd-X/mwtract/TASK2_Modeling/Data/'
                          'Time_Series/INTERNAL_CALIBRATION/Int_Cal_05182018/'
                          'GES_stage.csv'),
                         infer_datetime_format=True, # way faster
                         parse_dates=['Time'])
ges_april=pd.read_csv( ('/home/rusty/mirrors/Ed/Sync/UCD/Projects/CDFW_Arc/dflowfm/'
                        '27-comp_bedlevtyp_2_3_4/stations/wsel/'
                        'GES_STAGE_april2014.csv'),
                       parse_dates=['Time'])

sel=(ges_mwtract.Time>=ges_april.Time.iloc[0])&(ges_mwtract.Time<=ges_april.Time.iloc[-1])
##

from stompy.io.local import usgs_nwis

cache_dir='cache'
os.path.exists(cache_dir) or os.mkdir(cache_dir)

ges_ds=usgs_nwis.nwis_dataset(station="11447905",
                              start_date=np.datetime64("2014-04-01"),
                              end_date=np.datetime64("2014-05-01"),
                              cache_dir=cache_dir,
                              products=[60,65], # discharge and stage
                              days_per_request='M')

# that always returns UTC, so add a PST time:
ges_ds['time_pst']=('time',), ges_ds.time + np.timedelta64(-8,'h')

##
import matplotlib.pyplot as plt
plt.ion()

plt.figure(12).clf()

# These are all in agreement now.
plt.plot( ges_april.Time, ges_april.Stage,label='april')
plt.plot( ges_mwtract.Time[sel], ges_mwtract.Stage[sel],label='MWT')
plt.plot( ges_ds.time_pst, ges_ds.height_gage, label='NWIS')
plt.legend()

