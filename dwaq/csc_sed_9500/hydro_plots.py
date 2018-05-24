# show difference in bed stress when including waves
# if any...
import xarray as xr
import matplotlib.pyplot as plt

from stompy.grid import unstructured_grid
from stompy.model.delft import dfm_grid
##

nc_nowaves=xr.open_dataset("../dflowfm/csc_95/run-no_waves/DFM_OUTPUT_FlowFM/FlowFM_map.nc")
nc_waves=xr.open_dataset("../dflowfm/csc_95/DFM_OUTPUT_FlowFM/FlowFM_map.nc")


##

g=unstructured_grid.UnstructuredGrid.from_ugrid(nc_nowaves)

##
tau_nowaves=nc_nowaves.mesh2d_taus.isel(time=200).values

tau_waves=nc_waves.mesh2d_taus.isel(time=200).values

##

plt.figure(1).clf()
fig,ax=plt.subplots(num=1)

ccoll=g.plot_cells(values=tau_nowaves,ax=ax,cmap='jet')
ccoll.set_clim([0,0.2])
ax.axis('equal')

# Yeah - those are exactly the same.
