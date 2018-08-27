from stompy.plot import plot_utils
from stompy.model.delft import dfm_grid

##

g06=dfm_grid.DFMGrid('runs/20180807_grid98_06/CacheSloughComplex_v98_bathy_sparse_net.nc')
g07=dfm_grid.DFMGrid('runs/20180807_grid98_07/CacheSloughComplex_v98_bathy2_sparse_net.nc')

##

plt.figure(1).clf()
fig,ax=plt.subplots(1,1,num=1)

g06.plot_edges(color='k',lw=0.4)

ncoll=g06.plot_nodes(values=(g07.nodes['depth']-g06.nodes['depth']),
                     ax=ax)
plot_utils.cbar(ncoll,label='Depth 7 - depth 6')
ncoll.set_clim([-5,5])
ncoll.set_cmap('seismic')
ax.axis('equal')

