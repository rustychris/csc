import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse

from stompy.spatial import field,linestring_utils
from stompy import utils
from stompy.grid import depth_connectivity, unstructured_grid
import stompy.plot.cmap as scmap

from stompy.model.delft import dfm_grid
cmap=scmap.load_gradient('hot_desaturated.cpt')

##

def find_all_pairs(g,max_segs=7):
    """
    Automatically pull channel cross sections up to 7 segments wide.
    This does not use any hydro data, so relies on grid boundaries
    alone. For the CSC grid, it's right about 95% of the time, so
    don't assume that every channel section will be found, or that
    every node pair returns is a channel section.
    """

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

    for n in utils.progress(boundary_nodes,msg="find_all_pairs %s boundary nodes"):
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
    return all_pairs

##

def edge_depths_mean(g,dem,samples_per_edge=7,eta_max=2.5):
    """
    g: unstructured-grid
    dem: elevation data as Field
    samples_per_edge: edge mean estimated from this many samples
    eta_max: clip samples to <= eta_max.  This is to keep a high levee
      from swamping intertidal elevations.
    """
    alpha=np.linspace(0,1,samples_per_edge)[:,None,None] # for broadcasting

    x1=g.nodes['x'][g.edges['nodes'][:,0]]
    x2=g.nodes['x'][g.edges['nodes'][:,1]]

    seg_samples=(1-alpha)*x1 + (alpha)*x2
    seg_values=dem(seg_samples)
    seg_values=seg_values.clip(-np.inf,eta_max)

    dem_edge_depths=np.nanmean(seg_values,axis=0)

    return dem_edge_depths

# Sparse edge-mean fitting:
def node_depths_edge_mean_opt(g,
                              target_edge_depths,
                              target_node_depths,
                              section_weight=1.0,
                              node_weight=0.1,
                              nonsection_weight=0.0,
                              max_segments_per_section=7):
    """
    return per-node elevations such that the average of
    the endpoints per edge are close to the target edge depth.
    Edges which form channel cross-sections are given precedence.

    nodes per-edge to match the mean edge elevation from the DEM.
    """
    # construct a linear system where we try to solve good elevations for
    # all of the nodes at once

    # First, do this but ignore boundary edges:
    e2c=g.edge_to_cells()
    boundary_edges=e2c.min(axis=1)<0
    boundary_nodes=np.unique( g.edges['nodes'][boundary_edges] )

    weight_by_pairs=True

    if section_weight==nonsection_weight:
        log.info("Cross-section and non-cross-section weights the same.  No tracing needed")
    else:
        all_pairs=find_all_pairs(g,max_segments_per_section)
        cross_edges=[]
        for n1,n2 in utils.progress(all_pairs,msg="Extracting section edges %s"):
            cross_edges.append( g.shortest_path(n1,n2,return_type='edges') )
        cross_edges=np.concatenate(cross_edges)
        edge_weights=nonsection_weight*np.ones(g.Nedges(),np.float64)
        edge_weights[cross_edges]=section_weight

    rows=[]
    cols=[]
    data=[]
    rhs=[]

    for j in utils.progress(range(g.Nedges())):
        row=len(rhs)
        if boundary_edges[j]:
            continue
        if edge_weights[j]==0:
            continue
        n1,n2=g.edges['nodes'][j]
        rows.append(row)
        cols.append(n1)
        data.append(0.5 * edge_weights[j])
        rows.append(row)
        cols.append(n2)
        data.append(0.5 * edge_weights[j])
        rhs.append(target_edge_depths[j] * edge_weights[j])

    if node_weight>0:
        node_weights=node_weight*np.ones(g.Nnodes(),np.float64)
        # boundary nodes we consider free
        node_weights[boundary_nodes]=0

        for n in range(g.Nnodes()):
            if node_weights[n]!=0:
                rows.append(len(rhs))
                cols.append(n)
                data.append(node_weights[n])
                rhs.append(node_weights[n]*target_node_depths[n])

    M=sparse.coo_matrix((data,(rows,cols)),
                         shape=(len(rhs),g.Nnodes()))

    node_depths_sparse,status,itn,r1norm,r2norm,anorm,acond,arnorm,xnorm,var_=sparse.linalg.lsqr(M,rhs)

    return node_depths_sparse

##

def eval_pnt_pair(pnt_pair,node_depths):
    metrics={}
    nodes=metrics['nodes']=[g.select_nodes_nearest(xy) for xy in pnt_pair]
    node_path=metrics['node_path']=g.shortest_path(nodes[0],nodes[1])
    path_xy=metrics['path_xy']=g.nodes['x'][node_path]
    path_dist=metrics['path_dist']=utils.dist_along(path_xy)
    path_z=metrics['path_z']=node_depths[node_path]

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


def plot_pnt_pair(pnt_pair,node_depths):
    m=eval_pnt_pair(pnt_pair,node_depths=node_depths)

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


##

# DO THE WORK

# Load inputs:
g=unstructured_grid.UnstructuredGrid.from_ugrid('../grid/CacheSloughComplex_v100-edit06.nc')
dem=field.GdalGrid("../bathy/merged_2m-20181005.tif")

##

node_depths_dem=dem(g.nodes['x'])

##

if 1:
    target_edge_depths=edge_depths_mean(g,dem)
if 0:
    target_edge_depths=depth_connectivity.edge_connection_depth(g,dem,edge_mask=None,centers='lowest')

##
# this is kind of reasonable for BedLevType=4
# node_depths=depth_connectivity.greedy_edgemin_to_node(g,z_dem,edge_depths)

# This is trying to be reasonable for BedLevType=3
node_depths=node_depths_edge_mean_opt(g,target_edge_depths,node_depths_dem)

##


# some spot checks to see how that's doing.
if 0:
    plt.figure(10)
    ax=plt.gca()
    g.plot_edges(ax=ax)
    ax.axis('equal')

##
if 0:
    plt.figure(10)
    pnt_pair=plt.ginput(2)
    metrics,ax_g,ax_z=plot_pnt_pair(pnt_pair,node_depths)
##
if 0:
    # vertex locations defining ends of a transect:
    pnt_pairs=[
        [ [629710, 4250426], [629927, 4250402] ],
        [ [626423, 4238373], [626524, 4238469] ]
        ]


    pnt_pair=pnt_pairs[0]
    metrics,ax_g,ax_z=plot_pnt_pair(pnt_pair,node_depths_sparse)

    # an under-served edge
    metrics,ax_g,ax_z=plot_pnt_pair([(626690.0695316412, 4238152.963962159),
                                     (626766.1176475666, 4238269.110175572)],
                                    node_depths_sparse)


##
if 0: 
    all_pairs=find_all_pairs(g,7)
    #
    # Evaluate integral metrics on each cross-section
    edge_metrics=[None]*g.Nedges()

    xy_pairs=g.nodes['x'][all_pairs]

    for i,xy_pair in enumerate(xy_pairs):
        if i%100==0:
            print("%d/%d"%(i,len(xy_pairs)))
        m=eval_pnt_pair(xy_pair,node_depths=node_depths)
        nodes=m['node_path']
        for a,b in zip( nodes[:-1],nodes[1:] ):
            j=g.nodes_to_edge(a,b)
            assert j is not None
            edge_metrics[j]=m
##

if 0:
    edge_vals=np.nan*np.ones(g.Nedges())
    for j in range(g.Nedges()):
        m=edge_metrics[j]
        if m:
            edge_vals[j]=m['L_err']

    plt.figure(22).clf()
    fig,ax=plt.subplots(num=22)
    g.plot_edges(ax=ax,color='k',lw=0.3,alpha=0.5)
    ecoll=g.plot_edges(ax=ax,values=edge_vals,mask=np.isfinite(edge_vals))
    #ecoll.set_clim([-2,2])
    ecoll.set_clim([-100,100])
    ecoll.set_cmap('jet')
    ecoll.set_lw(4)
    ax.axis('equal')
    plt.colorbar(ecoll)

##

g.add_node_field('depth',node_depths,on_exists='overwrite')

# bathy2 => with dredged out CSC region
dfm_grid.write_dfm(g,'CacheSloughComplex_v100_bathy2_sparse_net.nc',overwrite=True)

##

g.add_cell_field('depth',g.interp_node_to_cell(node_depths),on_exists='overwrite')

g.write_cells_shp('derived/grid-cells.shp',overwrite=True)
