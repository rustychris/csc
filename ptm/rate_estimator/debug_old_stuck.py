# Go back to 36 second output to see more clearly how particles
# get stuck.
from stompy.grid import unstructured_grid
from stompy.model.fish_ptm import ptm_tools

##

hyd=xr.open_dataset('../../dflowfm/runs/20180807_grid98_16_single/ptm_hydro.nc')
g=unstructured_grid.UnstructuredGrid.from_ugrid(hyd)

##
init=ptm_tools.PtmBin('../scenario_divergence/DIST_bin.out')

ntimes=init.count_timesteps()
##

zoom=(605889.6569457075, 638002.2586920519, 4217801.158715993, 4241730.226468915)
# zoom=(597913.7274775933, 648118.8262812896, 4217179.54644355, 4301202.344200377)
# zoom=(611280.377359663, 632614.9072567355, 4222938.787804629, 4248182.140275016)
# zoom=(626037.7515578158, 626228.6109768279, 4232804.050163795, 4233029.878040465)

plt.figure(1).clf()
fig,ax=plt.subplots(num=1)

ti=5000

# Just plot stuck particles
tA,pA=init.read_timestep(ti)
tB,pB=init.read_timestep(ti+1)

stuck=(np.all(pA['x'][:,:2]==pB['x'][:,:2],axis=1))
inactive=pA['active']!=1

# init.plot(ti,ax=ax,zoom=zoom,update=False,ms=1)
ax.plot( pA['x'][stuck,0],
         pA['x'][stuck,1],
         'r.')

g.plot_edges(color='k',lw=0.4,ax=ax,clip=zoom) # ,labeler='id')
#g.plot_cells(centers=True, labeler=lambda i,r:str(i),
#             clip=zoom,ax=ax)

ax.axis(zoom)

# np.nonzero(stuck & utils.within_2d(pA['x'][:,:2],ax.axis()))
i=1962

##

traj_i=[]

for ti in utils.progress(range(4500,len(init.time))):
    t,p=init.read_timestep(ti)
    if len(p)>i:
        traj_i.append( [ti,p['x'][i,0],p['x'][i,1]] )

##
del ax.lines[-1]
del ax.collections[-1]

traj=np.array(traj_i)

plt.plot(traj[:,1],traj[:,2],'k-')
scat=plt.scatter(traj[:,1],traj[:,2],50,traj[:,0],lw=0,zorder=3, cmap='jet')


##

# For example,
j=25170
c_deep=21111
c_shallow=51090

##


t=hyd.nMesh2_data_time
# Flow on this edge is 0 for all time.
Qj=hyd.h_flow_avg.isel(nMesh2_edge=j,nMesh2_layer_3d=0)

# 1 for all time
j_bot=hyd.Mesh2_edge_bottom_layer.isel(nMesh2_edge=j)
# 0 for all time.
j_top=hyd.Mesh2_edge_top_layer.isel(nMesh2_edge=j)

# 0 for all time
Aj=hyd.Mesh2_edge_wet_area.isel(nMesh2_edge=j,nMesh2_layer_3d=0)

# shallow cell 51090 is always bottom layer=1, top=0
# deep cell 21111 is always bottom=top=1


##

plt.figure(2).clf()
plt.plot(t,Qj)


