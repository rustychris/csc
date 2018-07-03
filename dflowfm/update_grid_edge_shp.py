from stompy.model.delft import dfm_grid

g=dfm_grid.DFMGrid('CacheSloughComplex_v95_net.nc')

g.write_edges_shp('derived/grid-edges.shp',overwrite=True)
g.write_cells_shp('derived/grid-cells.shp',overwrite=True)
# g.write_nodes_shp('derived/grid-nodes.shp',overwrite=True)




