"""
more debugging of stuck particles
now see if a short run with 36s time steps and output
also gets stuck so quickly.

now with a streamlined debug run
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

from stompy.grid import unstructured_grid
from stompy.model.fish_ptm import ptm_tools

##

hyd=xr.open_dataset('../../dflowfm/runs/20180807_grid98_17/ptm_hydro.nc')
g=unstructured_grid.UnstructuredGrid.from_ugrid(hyd)

##
run_dir="run_debug_stuck"
init=ptm_tools.PtmBin(run_dir+'/TEST_bin.out')
ntimes=init.count_timesteps()
##

# zoom=(605889.6569457075, 638002.2586920519, 4217801.158715993, 4241730.226468915)
# zoom=(597913.7274775933, 648118.8262812896, 4217179.54644355, 4301202.344200377)
# zoom=(611280.377359663, 632614.9072567355, 4222938.787804629, 4248182.140275016)
# zoom=(626037.7515578158, 626228.6109768279, 4232804.050163795, 4233029.878040465)
zoom=(626114.3567959621, 626237.6438243476, 4232894.817537309, 4233010.536791836)

plt.figure(1).clf()
fig,ax=plt.subplots(num=1)

ti=20

init.plot(ti,ax=ax,zoom=zoom,update=False,ms=4)

g.plot_edges(color='k',lw=0.4,ax=ax,clip=zoom,labeler='id')
g.plot_cells(centers=True,labeler=lambda i,r:str(i),clip=zoom,ax=ax)
plt.setp(ax.texts,fontsize=6)

ax.axis(zoom)

##

# Looks like velocity is a little crazy in those cells, or maybe it's
# just the dispersion.  but the particles stick when they hit the boundary.
for ti in range(100):
    init.plot(ti,ax=ax,zoom=zoom,update=True,ms=4)
    plt.pause(0.05)

##

# For example,
j=25175
c_deep=51176
c_shallow=20356

##

t=hyd.nMesh2_data_time
hyd_short=hyd.isel(nMesh2_data_time=slice(48,48+10))

# Flow on this edge is 0 for all time.
Qj=hyd_short.h_flow_avg.isel(nMesh2_edge=j,nMesh2_layer_3d=0)

# 1 for all time
j_bot=hyd_short.Mesh2_edge_bottom_layer.isel(nMesh2_edge=j)
# 0 for all time.
j_top=hyd_short.Mesh2_edge_top_layer.isel(nMesh2_edge=j)

# 0 for all time
Aj=hyd_short.Mesh2_edge_wet_area.isel(nMesh2_edge=j,nMesh2_layer_3d=0)

# 1.8238m positive down
hyd_short.Mesh2_face_depth.isel(nMesh2_face=51176)

# eta 1.39, 1.365, 1.364
hyd_short.Mesh2_sea_surface_elevation.isel(nMesh2_face=51176)

##

# So there is a lot of diffusive random walk going on.
# What particle is getting stuck there?
# Could there be an issue with the in-polygon testing due to
# some issue with node order?
# Or maybe on these edges, where one side is wet, k_top should
# be equal to k_bottom, not k_bottom-1?
# can check untrim output for that

# t,part=init.read_timestep(20)
part_idx=0

traj= [ init.read_timestep(ti)[1]['x'][part_idx]
        for ti in range(20) ]

traj=np.array(traj)

ax.plot(traj[:,0],traj[:,1],'g-o')

# Would increasing the horizontal diffusion substeps help?
# NORMAL_VELOCITY_GRADIENT?
# HORIZONTAL_DIFFUSION_METHOD
# CONSTANT_HORIZONTAL_EDDY_DIFFUSIVITY?

# How is it supposed to deal with bouncing into a boundary?
# based on reading the code, it should be able to diffuse
# off of the boundary.
##

plt.figure(2).clf()
plt.plot(t,Qj)


##

untrim_hyd=xr.open_dataset('/home/rusty/mirrors/ucd-X/mwtract/'
                           'TASK2_Modeling/Scratch/dfm2ptm/untrim2ptm/'
                           'untrim_run/untrim_hydro_Jan2017.nc')
untrim_g=unstructured_grid.UnstructuredGrid.from_ugrid(untrim_hyd)

# 34:
# untrim_hyd.Mesh2_edge_bottom_layer.min()
# 0
# untrim_hyd.Mesh2_edge_top_layer.min()

# That's a bit confusing
# untrim_hyd.Mesh2_edge_bottom_layer.isel(nMesh2_data_time=2).min()

# there are times that j_top!=0, but edge wet area IS 0.
# maybe j_top should be zero only in cases where both cells
# are dry?

# okay - edge_tops are allowed to be set to 0, but cell tops
# never go below cell_bot.

# status at time 1:
ti=1
j_top=untrim_hyd.Mesh2_edge_top_layer.isel(nMesh2_data_time=ti).values
j_bot=untrim_hyd.Mesh2_edge_bottom_layer.isel(nMesh2_data_time=ti).values
i_top=untrim_hyd.Mesh2_face_top_layer.isel(nMesh2_data_time=ti).values
i_bot=untrim_hyd.Mesh2_face_bottom_layer.isel(nMesh2_data_time=ti).values

fig=plt.figure(10)
fig.clf()
fig,ax=plt.subplots(num=10)
ecoll=untrim_g.plot_edges(ax=ax,values=j_top-j_bot,lw=3,cmap='jet')
ecoll.set_clim([-2,2])
ccoll=untrim_g.plot_cells(ax=ax,values=i_top-i_bot,cmap='jet')
ccoll.set_clim([-2,2])

ax.axis('equal')

i=9552

##

# Something weird is going on with cell indices
hyd=xr.open_dataset('../../dflowfm/runs/20180807_grid98_17/ptm_hydro.nc')
g=unstructured_grid.UnstructuredGrid.from_ugrid(hyd)

g2=unstructured_grid.UnstructuredGrid.from_ugrid('../../dflowfm/runs/20180807_grid98_17/'
                                                 'DFM_OUTPUT_flowfm/flowfm_map.nc')
# half of the cell centers do not agree.
# first bad one is cell 4.  but that's just a matter of a centimeter.
# just difference in how cc is calculated.
# maybe ds was stale?  largest difference is 23m.
