# Take a look at how the thacker test case is setup, in case
# 1D channel networks are a good way to go for the CSC
# model.

base="/home/rusty/src/dfm/1.5.0/validation/cases/f03_advection/c040_thacker1d_standard/"

ds_net=xr.open_dataset(os.path.join(base,'planar1d_net.nc'))

##
import stompy.grid.unstructured_grid as ugrid

g=ugrid.UnstructuredGrid.read_dfm(ds_net)

##

fig=plt.figure(2)
fig.clf()
ax=fig.gca()

# This includes the 1D area.
g.plot_edges(ax=ax,values=g.edges['NetLinkType'])
ax.axis('equal')

##

# cells (=elements) only in the grids, not the network
g.plot_cells(ax=ax,color='0.7',lw=0)

##

# An elevation coinciding with nodes of each
xyb=np.loadtxt(os.path.join(base,'inibott.xyb'))
scat=ax.scatter(xyb[:,0],xyb[:,1],20,xyb[:,2])
plt.colorbar(scat)

##
