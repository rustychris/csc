"""
Look more closely around FPX to understand the differences in setting bathy.
"""
import xarray as xr
from stompy.grid import unstructured_grid

##

# Load the bathy, see how the area compares
dem=field.GdalGrid(os.path.join("/home/rusty/mirrors/ucd-X",
                                "Arc_Hydro/CSC_Project/MODELING/1_Hydro_Model_Files",
                                "Geometry/Bathymetry_Tiles/NDelta_hydro_2m_v4",
                                "ndelta_hydro_2m_v4z.tif"))
# There is some weird stuff with the extents of that DEM.
dem.extents[0]=float(int(dem.extents[0])) # drops 0.80m
dem.extents[1]=float(int(dem.extents[1])) # same.
dem.extents[2]=float(int(dem.extents[2])) # drops 0.88m
dem.extents[3]=float(int(dem.extents[3])) # same.

##

runs=[ "runs/base20180701/DFM_OUTPUT_FlowFM",
       "runs/20180712_rough027_blt4_bathy01/DFM_OUTPUT_flowfm",
       "runs/base20180701-xs_align/DFM_OUTPUT_FlowFM",
       "runs/20180725_test_xs_align/DFM_OUTPUT_flowfm",
]

hists= [ "runs/base20180701/DFM_OUTPUT_FlowFM/FlowFM_his.nc",
         "runs/20180712_rough027_blt4_bathy01/DFM_OUTPUT_flowfm/flowfm_0000_his.nc",
         "runs/base20180701-xs_align/DFM_OUTPUT_FlowFM/FlowFM_his.nc",
         "runs/20180725_test_xs_align/DFM_OUTPUT_flowfm/flowfm_0000_his.nc"]

hist_dss=[xr.open_dataset(fn) for fn in hists]

map_ds=xr.open_dataset( os.path.join(runs[0],"FlowFM_map.nc") )

##
cross_section_name="Freeport"

# 1
cs_idx=list(hist_dss[3].cross_section_name.values).index(b'Freeport')

stn_idx=list(hist_dss[3].station_name.values).index(b'FPX')

#   cross_section_x_coordinate                    (cross_section, cross_section_pts) float64 ...
#   cross_section_y_coordinate                    (cross_section, cross_section_pts) float64 ...
#   cross_section_discharge                       (time, cross_section) float64 ...
#   cross_section_cumulative_discharge            (time, cross_section) float64 ...
#   cross_section_area                            (time, cross_section) float64 ...
#   cross_section_velocity                        (time, cross_section) float64 ...
#   cross_section_salt                            (time, cross_section) float64 ...
#   cross_section_temperature                     (time, cross_section) float64 ...

cs_x= hist_dss[3].cross_section_x_coordinate.isel(cross_section=cs_idx).values
cs_y= hist_dss[3].cross_section_y_coordinate.isel(cross_section=cs_idx).values
valid=(cs_x<1e10)
cs_xy=np.c_[cs_x[valid],cs_y[valid]]


##

g=csc_grid=unstructured_grid.UnstructuredGrid.from_ugrid(map_ds)

##
# zoom=(629571.6057946421, 631493.473644821, 4256748.848504359, 4258180.950031428)

zoom=[cs_xy[:,0].min()-300,cs_xy[:,0].max()+300,
      cs_xy[:,1].min()-300,cs_xy[:,1].max()+300]

fig=plt.figure(1)
fig.clf()
ax=fig.add_subplot(1,1,1)

csc_grid.plot_edges(ax=ax)
ax.plot(cs_xy[:,0],cs_xy[:,1],'g-',label='Cross section')
ax.axis('equal')
ax.axis(zoom)

##

from stompy.spatial import linestring_utils

# Pull profile from DEM along those lines:
cs_line=linestring_utils.resample_linearring(cs_xy,2.0,closed_ring=False)
cs_dist=utils.dist_along(cs_line)
# this cross-section includes levees -- skip outside the levees
# channel_sel=(cs_dist>=106)&(cs_dist<=345)
# or with the aligned cross-section:
channel_sel=(cs_dist>=61)&(cs_dist<=299)


cs_line_z=dem(cs_line)

plt.figure(3).clf()
fig,ax=plt.subplots(num=3)

ax.plot(cs_dist,cs_line_z,'k-')

sel_z=cs_line_z[channel_sel]
sel_d=cs_dist[channel_sel]
ax.plot(sel_d,sel_z,'g-',lw=2)
dd=np.median(np.diff(sel_d))# 2.00
def h_to_A(h):
    return (dd*(h-sel_z).clip(0,np.inf)).sum()

##

plt.figure(2).clf()
fig_time,(ax_Q,ax_A,ax_h)=plt.subplots(3,1,sharex=True,num=2)

colors=['r','g','b','orange']

for i,hist_ds in enumerate(hist_dss):
    col=colors[i]

    ax_Q.plot(hist_ds.time,
              hist_ds.cross_section_discharge.isel(cross_section=cs_idx),
              color=col,
              label='Run %d'%i)
    ax_A.plot(hist_ds.time,
              hist_ds.cross_section_area.isel(cross_section=cs_idx),
              color=col,
              label='Run %d'%i)

    ax_A.plot(hist_ds.time,
              [ h_to_A(h) for h in hist_ds.waterlevel.isel(stations=stn_idx).values],
              color=col,ls='--',
              label='Run %d from DEM'%i)

    ax_h.plot(hist_ds.time,hist_ds.waterlevel.isel(stations=stn_idx),
              color=col,
              label='Run %d'%i)
ax_Q.legend()
ax_A.legend()
ax_h.legend()

ax_Q.set_ylabel('Flow (m3s)')
ax_A.set_ylabel('Area (m2)')
ax_h.set_ylabel('waterlevel (m)')

##
tile=dem.extract_tile(zoom)
##
# running this again, but with XS coordinates tweaked to better align with the grid.

# What does the grid geometry look like in cross-section?
if 0:
    g=csc_grid
    hist=hist_dss[2]
    blt=3
else:
    # Load the same details but for the aligned blt4 run.
    this_map_ds=xr.open_dataset( os.path.join(runs[3],"flowfm_0000_map.nc") )
    g=unstructured_grid.UnstructuredGrid.from_ugrid(this_map_ds)
    hist=hist_dss[3]
    blt=4

nodes=[g.select_nodes_nearest(p) for p in cs_xy]
node_path=np.concatenate( [g.shortest_path(nodes[i],nodes[i+1])
                           for i in range(len(nodes)-1)] )
node_xy=g.nodes['x'][node_path]

plt.figure(4).clf()
fig,(ax_map,ax_xs,ax_A)=plt.subplots(1,3,num=4)
tile.plot(ax=ax_map,cmap='gray')
g.plot_edges(ax=ax_map,clip=zoom,lw=0.5,color='k')
ax_map.plot(cs_xy[:,0],cs_xy[:,1],'g-o',label='Cross section as input')
ax_map.plot(node_xy[:,0],node_xy[:,1],'r-o',label='Snapped to grid')
ax_map.axis('equal')

node_xy_resamp=linestring_utils.resample_linearring(node_xy,2.0,closed_ring=0)
d_resamp=utils.dist_along(node_xy_resamp)
z_resamp=dem(node_xy_resamp)

ax_xs.plot(d_resamp,z_resamp,label='DEM elevation along grid edges')

d_nodes=utils.dist_along(node_xy)
z_path=g.nodes['mesh2d_node_z'][node_path]

ax_xs.plot(d_nodes,z_path,
           'b-o',label='Grid elevation at nodse')

if blt==3:
    # for bedlevtype=3, can I get back to DFM's cross-sectional area? yes.
    edge_z=0.5*( z_path[:-1] + z_path[1:] )
elif blt==4:
    edge_z=np.minimum(z_path[:-1],z_path[1:])

edge_L=np.diff(d_nodes)


hist_h=hist.waterlevel.isel(stations=stn_idx).values

edge_A=np.array( [ (edge_L*(h-edge_z).clip(0,np.inf)).sum()
                   for h in hist_h] )

hist_A=hist.cross_section_area.isel(cross_section=cs_idx).values

ax_A.scatter(hist_A,edge_A)
dom=[min(edge_A.min(),hist_A.min()),
     max(edge_A.max(),hist_A.max())]

ax_A.plot(dom,dom,'k-',lw=0.5)
ax_A.set_xlabel("History A")
ax_A.set_ylabel("Calc A")



##
g_run0=unstructured_grid.UnstructuredGrid.from_ugrid(map_ds)

plt.figure(10).clf()

g_run0.plot_edges(lw=0.5,color='k')
scat=g_run0.plot_nodes(values=g_run0.nodes['mesh2d_node_z'])
plt.colorbar(scat)
plt.axis('equal')




