from stompy.io.local import usgs_nwis
from stompy import utils

import os
import matplotlib.pyplot as plt
##

cache_dir="cache"

if not os.path.exists(cache_dir):
    os.makedirs(cache_dir)

    
ds=usgs_nwis.nwis_dataset(11455420,
                          np.datetime64("2019-08-01"),
                          np.datetime64("2019-09-01"),
                          products=[65],
                          cache_dir=cache_dir)

ds_long=usgs_nwis.nwis_dataset(11455420,
                               np.datetime64("2019-07-01"),
                               np.datetime64("2019-10-01"),
                               products=[65],
                               cache_dir=cache_dir)

##

# Try basic tidal filling:
height_filled=utils.fill_tidal_data(ds.height_gage)

height_filled_long=utils.fill_tidal_data(ds_long.height_gage)



##
plt.figure(1).clf()

plt.plot(ds.time, ds.height_gage,label="original" )
plt.plot(height_filled.time, height_filled,label="fill, 1 mo")
plt.plot(height_filled_long.time, height_filled_long,label="fill, 3 mos")
plt.legend()

# Stitching is not great -- but a lowpass on tidal data is a good idea
# regardless, so just lowpass it after the stitch.

