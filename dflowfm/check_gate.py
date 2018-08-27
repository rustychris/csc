import xarray as xr
import numpy as np
from stompy.grid import unstructured_grid

import matplotlib.pyplot as plt

##

ds3=xr.open_dataset('runs/20180807_grid98_13/DFM_OUTPUT_flowfm/flowfm_0003_map.nc')

g3=unstructured_grid.UnstructuredGrid.from_ugrid(ds3)

##
zoom0=(604981.1101668236, 606079.0483167086, 4234753.578142932, 4235571.719151394)
zoom1=(605562.8808368467, 605697.4585567462, 4235272.981963599, 4235373.264071008)
plt.figure(10).clf()
ax=plt.gca()

ti=146
ecoll=g3.plot_edges(ax=ax,values=np.abs(ds3.mesh2d_u1.isel(time=ti).values),
                    lw=4,cmap='jet')
ccoll=g3.plot_cells(ax=ax,values=ds3.mesh2d_s1.isel(time=ti).values,
                    cmap='jet')

ecoll.set_clim([0,0.5])
ccoll.set_clim([0,2])

ax.axis('equal')
ax.axis(zoom1)
plt.colorbar(ecoll,label='Edge |u|')
plt.colorbar(ccoll,label='Cell eta')

##

#mesh2d_q1
#mesh2d_u1
# mesh2d_hu: water depth at velocity points.
