# Get a sense of how to fabricate Decker Island tides when only Rio Vista
# data is available.
# choosing not to do this -- better to truncate the domain at Rio Vista.

import numpy as np
from stompy.io.local import usgs_nwis
from stompy import utils
import matplotlib.pyplot as plt

## 
period=[np.datetime64("2016-01-01"),
        np.datetime64("2019-01-01")]

decker=usgs_nwis.nwis_dataset(station=11455478,
                              start_date=period[0], end_date=period[1],
                              products=[60,65], cache_dir='cache')

riovista=usgs_nwis.nwis_dataset(station=11455420,
                                start_date=period[0], end_date=period[1],
                                products=[60,65],cache_dir='cache')

##


##
from stompy import filters

# separate into tidal, subtidal
for ds in [decker,riovista]:
    da_fill=utils.fill_tidal_data(ds['height_gage'])
    ds['ftime']=('ftime',), da_fill.time
    ds['stage_fill']=('ftime',),da_fill.values
    ds['stage_lp']=('ftime',),filters.lowpass(ds['stage_fill'].values,
                                              (ds.ftime.values-ds.ftime.values[0])/np.timedelta64(1,'s'),
                                              cutoff=40*3600)
    ds['stage_hp']=ds.stage_fill - ds.stage_lp
    
## 
# Find the tidal lag:

lag_hp=utils.find_lag_xr(decker.stage_hp,riovista.stage_hp)

# Decker leads Rio Vista by 1738s
lag_hp_s = lag/np.timedelta64(1,'s')

# and subtidal lag is almost exactly 2 hours.  Weird.
lag_lp=utils.find_lag_xr(decker.stage_lp,riovista.stage_lp)
lag_lp_s=lag_lp/np.timedelta64(1,'s')

## 


plt.figure(1).clf()

plt.plot(decker.time,decker.height_gage,label='Decker')
plt.plot(riovista.time,riovista.height_gage,label='Rio Vista')

#plt.plot(decker.ftime,decker.stage_lp,label='Decker')
#plt.plot(riovista.ftime,riovista.stage_lp,label='Rio Vista')

plt.plot(riovista.ftime+lag_hp,
         riovista.stage_fill,
         label="RV lagged")


plt.legend()

##
