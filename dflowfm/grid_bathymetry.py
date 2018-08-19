import numpy as np
import matplotlib.pyplot as plt

from stompy.spatial import field,linestring_utils
import stompy.plot.cmap as scmap

from stompy.model.delft import dfm_grid
cmap=scmap.load_gradient('hot_desaturated.cpt')

##

g=unstructured_grid.UnTRIM08Grid('../grid/CacheSloughComplex_v98.grd')

dem=field.GdalGrid("../bathy/merged_2m.tif")

##

z_dem=dem(g.nodes['x'])
node_depths=z_dem
##

from stompy.grid import depth_connectivity
edge_depths=depth_connectivity.edge_connection_depth(g,dem,edge_mask=None,centers='lowest')
node_depths=depth_connectivity.greedy_edgemin_to_node(g,z_dem,edge_depths)

##

# Write a new version of the grid with this bathy
g.add_node_field('depth',node_depths,on_exists='overwrite')

dfm_grid.write_dfm(g,'CacheSloughComplex_v98_bathy_net.nc',overwrite=True)

##

# some spot checks to see how that's doing.
plt.figure(10)
ax=plt.gca()
g.plot_edges(ax=ax)
ax.axis('equal')

##
# vertex locations defining ends of a transect:
pnt_pairs=[
    [ [629710, 4250426], [629927, 4250402] ],
    [ [626423, 4238373], [626524, 4238469] ]
    ]

##

def eval_pnt_pair(pnt_pair):
    metrics={}
    nodes=metrics['nodes']=[g.select_nodes_nearest(xy) for xy in pnt_pair]
    node_path=metrics['node_path']=g.shortest_path(nodes[0],nodes[1])
    path_xy=metrics['path_xy']=g.nodes['x'][node_path]
    metrics['path_dist']=utils.dist_along(path_xy)
    metrics['path_z']=node_depths[node_path]

    resamp_xy=linestring_utils.resample_linearring(path_xy,
                                                   2.0,closed_ring=False)
    metrics['resamp_xy']=resamp_xy

    metrics['resamp_z']=resamp_z=dem(resamp_xy)
    metrics['resamp_dist']=resamp_dist=utils.dist_along(resamp_xy)

    # edge depths with bedlevtyp3, sampled to higher resolution
    seg_z_typ3=0.5*(path_z[:-1]+path_z[1:])
    resamp_z_typ3=seg_z_typ3[ np.searchsorted(path_dist,resamp_dist).clip(1,len(path_dist)-1)-1]
    metrics['resamp_z_typ3']=resamp_z_typ3

    metrics['ref_eta']=ref_eta=1.0

    A_dem =np.trapz( (ref_eta-resamp_z).clip(0,np.inf),
                    resamp_dist)
    L_dem=np.trapz( 1.*(ref_eta>resamp_z),
                    resamp_dist)
    A_typ3=np.trapz( (ref_eta-resamp_z_typ3).clip(0,np.inf),
                     resamp_dist)
    L_typ3=np.trapz( 1.*(ref_eta>resamp_z_typ3),
                     resamp_dist)

    metrics['A_dem']=A_dem
    metrics['L_dem']=L_dem
    metrics['A_typ3']=A_typ3
    metrics['L_typ3']=L_typ3

    metrics['A_err']=A_typ3-A_dem
    metrics['L_err']=L_typ3-L_dem
    metrics['z_err']=metrics['A_err']/(L_dem+0.1)

    return metrics


def plot_pnt_pair(pnt_pair):
    m=eval_pnt_pair(pnt_pair)

    ctr=np.array(pnt_pair).mean(axis=0)
    dist=utils.dist(pnt_pair[0],pnt_pair[1])

    zoom=[ctr[0]-2*dist,ctr[0]+2*dist,
          ctr[1]-2*dist,ctr[1]+2*dist]

    # analyze a pair of points
    plt.figure(11).clf()
    fig,axs=plt.subplots(1,2,num=11)

    ax_g,ax_z=axs

    g.plot_edges(ax=ax_g,clip=zoom)

    ax_g.plot( m['path_xy'][:,0],m['path_xy'][:,1],'g-o')
    ax_g.axis('equal')
    ax_g.axis(zoom)

    ax_z.plot(m['path_dist'],m['path_z'],
              'g-o',label='node elevations')

    ax_z.plot(m['resamp_dist'],m['resamp_z'],'k-',label='DEM')
    ax_z.plot(m['resamp_dist'],m['resamp_z_typ3'],'b-',label='typ 3')

    ax_z.axhline(m['ref_eta'],color='gray',lw=0.5,
                 label='Ref eta=%.1f'%m['ref_eta'])

    ax_z.legend()

    lines=["XS area error\n%.1f m2\n%.1f%%"%( m['A_err'], 100*m['A_err']/m['A_dem']),
           "Length\n%.1f m\n%.1f%%"%( m['L_err'],100*m['L_err']/m['L_dem'])]
    ax_z.text(0.01,0.15,
              "\n".join(lines),
              transform=ax_z.transAxes)
    return m,ax_g,ax_z

pnt_pair=pnt_pairs[0]
metrics,ax_g,ax_z=plot_pnt_pair(pnt_pair)

##

# Automatically pull channel cross sections up to 7 segments wide
max_segs=7


e2c=g.edge_to_cells()
e_boundary=np.any( e2c<0, axis=1)
boundary_nodes=np.unique( g.edges['nodes'][e_boundary] )
boundary_mask=np.zeros(g.Nnodes(),np.bool)
boundary_mask[boundary_nodes]=True
##
node_marks=np.zeros(g.Nnodes(),np.bool8)
all_pairs=[]

def boundary_weight(j):
    if e_boundary[j]:
        return 1.0
    else:
        return np.nan
def internal_weight(j):
    if e_boundary[j]:
        return np.nan
    else:
        return 1.0

for n in boundary_nodes:
    if node_marks[n]:
        continue
    node_marks[n]=True

    # search only boundary edges to rule out along-boundary neighbors
    my_nbrs=g.shortest_path(n,n2=boundary_nodes,
                            edge_weight=boundary_weight,
                            max_return=2*max_segs)
    my_nbrs=[ mn[0] for mn in my_nbrs]
    other_n2=np.setdiff1d(boundary_nodes,my_nbrs)
    tran_path=g.shortest_path(n,n2=other_n2,
                              edge_weight=internal_weight,
                              max_return=1)
    if len(tran_path)==0:
        continue
    n2=tran_path[0][1][0]
    all_pairs.append( [n,n2] )
    # this isn't a strictly commutative property, but close enough
    # for our purposes
    node_marks[n2]=True

# g.plot_nodes(mask=my_nbrs,ax=ax_g,color='g')
#tran_xy= g.nodes['x'][tran_path[0][1]]
#ax_g.plot(tran_xy[:,0],tran_xy[:,1],'r-o')

##
# separate approach to get node elevations to best reflect edge averages.

edge_vals=np.nan*np.zeros(g.Nedges(),np.float64)

xy_pairs=g.nodes['x'][all_pairs]

for i,xy_pair in enumerate(xy_pairs):
    if i%100==0:
        print("%d/%d"%(i,len(xy_pairs)))
    m=eval_pnt_pair(xy_pair)
    nodes=m['node_path']
    for a,b in zip( nodes[:-1],nodes[1:] ):
        j=g.nodes_to_edge(a,b)
        assert j is not None
        edge_vals[j]=m['z_err']

##

plt.figure(12).clf()
fig,ax=plt.subplots(num=12)
g.plot_edges(ax=ax,color='k',lw=0.3,alpha=0.5)
ecoll=g.plot_edges(ax=ax,values=edge_vals,mask=np.isfinite(edge_vals))

ax.axis('equal')

plt.colorbar(ecoll)

##

plt.figure(13).clf()

plt.hist(edge_vals[np.isfinite(edge_vals)],bins=np.linspace(-8,8,100))

# seems that there is something about the channel shape on the mid sac
# that causes a lot of under representation of cross-section.
