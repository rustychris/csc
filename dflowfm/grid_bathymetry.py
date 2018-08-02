import numpy as np
import matplotlib.pyplot as plt

from stompy.spatial import field
import stompy.plot.cmap as scmap
from stompy.model.delft import dfm_grid
cmap=scmap.load_gradient('hot_desaturated.cpt')

##

g=dfm_grid.DFMGrid('CacheSloughComplex_v95_net.nc')

# Also compared to DWR bathy database, but that was significantly
# different.
X="/home/rusty/mirrors/ucd-X"
dem=field.GdalGrid(os.path.join(X,
                                "Arc_Hydro/CSC_Project/MODELING/1_Hydro_Model_Files",
                                "Geometry/Bathymetry_Tiles/NDelta_hydro_2m_v4",
                                "ndelta_hydro_2m_v4z.tif"))
# There is some weird stuff with the extents of that DEM.
dem.extents[0]=float(int(dem.extents[0])) # drops 0.80m
dem.extents[1]=float(int(dem.extents[1])) # same.
dem.extents[2]=float(int(dem.extents[2])) # drops 0.88m
dem.extents[3]=float(int(dem.extents[3])) # same.

##

z_dem=dem(g.nodes['x'])

##

plt.figure(1).clf()
fig,axs=plt.subplots(1,2,num=1,sharex=True,sharey=True)

for ax in axs:
    g.plot_edges(ax=ax,lw=0.5,color='k')


ncoll1=g.plot_nodes(ax=axs[0],values=g.nodes['depth'],sizes=40)
ncoll2=g.plot_nodes(ax=axs[1],values=z_dem,sizes=40)

plt.setp( [ncoll1,ncoll2],clim=[-4,3],cmap=cmap,zorder=1)

##

plt.figure(2).clf()
plt.plot(g.nodes['depth'], z_dem, 'b.',alpha=0.3,ms=1)

##

plt.figure(3).clf()
fig,ax=plt.subplots(num=3)
g.plot_edges(ax=ax,lw=0.5,color='k',zorder=-1)
ncoll=g.plot_nodes(ax=ax,values=z_dem-g.nodes['depth'],sizes=40)
ncoll.set_cmap('seismic')
ncoll.set_clim([-1,1])
ax.axis('equal')
plt.colorbar(ncoll,label='DEM-grid (m)')

##

# Run my bathy script on this data, see if that gets me back to Thomas'
from stompy.grid import depth_connectivity
edge_depths=depth_connectivity.edge_connection_depth(g,dem,edge_mask=None,centers='lowest')

node_depths=depth_connectivity.greedy_edgemin_to_node(g,z_dem,edge_depths)

##

plt.figure(4).clf()
# plt.plot(g.nodes['depth'], z_dem, 'b.',alpha=0.3,ms=1)
plt.plot(g.nodes['depth'], node_depths, 'k.',alpha=0.3,ms=1)

##

valid=np.isfinite( g.nodes['depth']+node_depths+z_dem)

# so the above is slightly closer to the grid bathy, but not quite there.
Ra=np.corrcoef(g.nodes['depth'][valid],
               node_depths[valid])[0,1]
Rb=np.corrcoef(g.nodes['depth'][valid],
               z_dem[valid])[0,1]
# z_dem has 4 nan?  This is from the southern edge of the DEM being corrupted

##

# Where/how different is the grid and the connectivity bathy?
plt.figure(5).clf()
fig,ax=plt.subplots(num=5)
g.plot_edges(ax=ax,lw=0.5,color='k',zorder=-1)
ncoll=g.plot_nodes(ax=ax,values=node_depths-g.nodes['depth'],sizes=40)
ncoll.set_cmap('seismic')
ncoll.set_clim([-1,1])
ax.axis('equal')
plt.colorbar(ncoll,label='Connect-grid (m)')

##

# Write a new version of the grid with this bathy
g.nodes['depth']=node_depths

dfm_grid.write_dfm(g,'CacheSloughComplex_v95_bathy01_net.nc')


