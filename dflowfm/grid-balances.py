g=unstructured_grid.UnstructuredGrid.read_ugrid("CacheSloughComplex_v111-edit19fix-bathy.nc")
dws_line=np.array( [(615963.7641646799, 4233147.324837522), (616387.435913091, 4233066.625456871)] )
left_cut=g.select_cells_by_cut(dws_line)
gcut=g.copy()

## 
for c in np.nonzero(~left_cut)[0]:
    print(c)
    gcut.delete_cell(c)
##

gcut.renumber()

##

gcut.plot_cells(color='orange')

##

dws_poly=g.boundary_polygon_by_union(left_cut)
##

dws_points=np.array( dws_poly.exterior )
print("%d 2"%len(dws_points) )

for i in dws_points:
    print("%.3f %.3f"%(i[0],i[1]))

##



    
