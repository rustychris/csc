# Extract a DEM from a node-centered grid

# be borrowed.
import matplotlib.pyplot as plt

from stompy.model.delft import dfm_grid
from stompy.grid import unstructured_grid
from stompy.plot import plot_utils

##

sch=unstructured_grid.UnstructuredGrid.from_ugrid('../bathy/grid-sources/schout_161.nc')
sch_ds=xr.open_dataset('../bathy/grid-sources/schout_161.nc')
poly=sch.boundary_polygon() # 10s

node_values=sch.nodes['depth']

##
from stompy.spatial import constrained_delaunay

g=sch
# slow!
cdf=constrained_delaunay.ConstrainedXYZField( g.nodes['x'],
                                              g.nodes['depth'],
                                              edges=g.edges['nodes'] )

##
bounds=(618784.5020805883, 619126.090221111, 4293373.8480339, 4293788.271474491)
bounds_sm=(618860.3632124562, 618991.2337552443, 4293473.6282357, 4293688.171493264)
bounds_lg=(618090., 621066, 4291272, 4294714)
bounds_nsac=(608680, 632020., 4271024, 4298020)
# pretty slow to get triangulation

# this huge area takes maybe 60s
dem0=cdf.to_grid(bounds=bounds_nsac,dx=2,dy=2)

##

mask=dem0.polygon_mask(poly)
dem0.F[~mask] = np.nan

##

# maybe another 60s
dem0.smooth_by_convolution(iterations=5)

##

clim=[-2,6]
g=sch
plt.figure(21).clf()
fig,ax=plt.subplots(num=21)

dem0.plot(ax=ax,vmin=clim[0],vmax=clim[1],cmap='jet')

g.plot_edges(ax=ax,clip=bounds,color='k',lw=0.3)
ncoll=g.plot_nodes(values=g.nodes['depth'],ax=ax,clip=bounds,
                   cmap='jet')
ncoll.set_clim(clim)


##

dem0.write_gdal("schism_gridded_2m.tif")

