"""
Another shot at dispersion, but using a distance metric restricted to be within the
grid.

Looking for ways to refine grid-based approach.
 - ignore particles that leave the domain.

 - would it make sense to reverse the process, instead looking at the collection
   of the particles at a given location at t_n, and where they came from at t_m<t_n.
   also look at distance traveled by centroid in addition to increase in variance.
"""
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
import heapq

from shapely.ops import cascaded_union

from stompy import utils
from stompy.grid import unstructured_grid
from stompy.model.fish_ptm import ptm_tools
from stompy.plot import plot_wkb
from stompy.spatial import proj_utils

import logging
utils.log.setLevel(logging.INFO)

##

g=unstructured_grid.UnstructuredGrid.from_ugrid('../../dflowfm/runs/20180807_grid98_17/ptm_hydro.nc')
grid_poly=g.boundary_polygon_by_edges()
ll2utm=proj_utils.mapper('WGS84','EPSG:26910')

##

ptm_run_dir="../rate_estimator/run_10days"
ptm_groups=["INIT","SAC","SRV"]

ptm_data=[ ptm_tools.PtmBin(os.path.join(ptm_run_dir,grp+"_bin.out"))
           for grp in ptm_groups ]

ntimes=ptm_data[0].count_timesteps()

run_start=ptm_data[0].read_timestep(0)[0]
run_stop =ptm_data[0].read_timestep(ntimes-1)[0]

##

ptm_group=ptm_data[0] # initial release

##

# build up just the 2D coordinates for all of the initial release particles
part_locs=[]
part_ts=[]
for ti in utils.progress(range(0,1000,4)):
    t,parts=ptm_group.read_timestep(ti)
    part_locs.append( parts['x'][:,:2].copy() )
    part_ts.append(t)

# This will fail is there are additional releases in this group
part_locs=np.array(part_locs)
part_ts=np.array([utils.to_dt64(t) for t in part_ts])

##

from scipy import sparse
from scipy.sparse import csgraph

def grid_to_graph(g):
    # use cell-to-cell connections to make it easier to avoid
    # hopping over land
    e2c=g.edge_to_cells()
    internal=np.all(e2c>=0,axis=1)
    c1=e2c[internal,0]
    c2=e2c[internal,1]
    cc=g.cells_center()
    lengths=utils.dist(cc[c1],cc[c2])
    bi_lengths=np.concatenate([lengths,lengths])
    bi_lengths=bi_lengths.astype(np.float32) # more than enough precision.
    A=np.concatenate((c1,c2))
    B=np.concatenate((c2,c2))
    graph=sparse.csr_matrix( (bi_lengths,(A,B)),
                             shape=(g.Ncells(),g.Ncells()) )
    # use scipy graph algorithms to find the connections
    # dists=csgraph.shortest_path(graph,directed=False,indices=cloud_nodes)
    # that ends up being (32,Nnodes), where 32 is the number of unique
    # nodes we started with
    return graph

graph=grid_to_graph(g)
##

particles=part_locs

# one-time preprocessing, find the cell corresponding to each particle
particle0_cells=[g.select_cells_nearest(particles[0,i,:])
                for i in utils.progress(range(len(particles[0,:,:])))]
particle0_cells=np.array(particle0_cells)

##

def particle_cloud_for_point(x0,count=10):
    c0=g.select_cells_nearest(x0)

    # get ordering of all other cells relative to c0
    c0_dists=csgraph.shortest_path(graph,directed=False,indices=[c0])[0,:]

    # in-grid distance for all particles at time 0
    part_dists=c0_dists[particle0_cells]
    cloud_particles=np.argsort(part_dists)[:count]
    return cloud_particles

##

# a list of hashes for pair-wise distance was pretty slow
# full calculation also too slow/large
pair_map=[dict() for i in range(g.Ncells())]

def pairwise_grid_distance(cells):
    n_extra=10*len(cells)
    missing=[]
    for ai,a in enumerate(cells):
        a_dists=pair_map[a]
        for bi,b in enumerate(cells[ai+1:],start=ai+1):
            if b not in a_dists:
                missing.append(b)

    if missing:
        dists=csgraph.shortest_path(graph,directed=False,indices=missing)
        for mi,m in enumerate(missing):
            for n in cells:
                pair_map[n][m]=pair_map[m][n]=dists[mi,n]
            # opportunistically grab more?
            extras=np.argsort(dists[mi,:])[:n_extra]
            for n in extras:
                pair_map[n][m]=pair_map[m][n]=dists[mi,n]
    else:
        print("No need to call shortest path")

    result=np.zeros( (len(cells),len(cells)), np.float64)
    for ai,a in enumerate(cells):
        a_dists=pair_map[a]
        for bi,b in enumerate(cells[ai+1:],start=ai+1):
            assert b in a_dists
            result[ai,bi]=result[bi,ai]=a_dists[b]
    return result


##
from stompy import memoize
@memoize.memoize(lru=25000)
def one_cell_dists(c):
    return csgraph.shortest_path(graph,directed=False,indices=[c])[0,:]

def pairwise_grid_distance(cells):
    result=np.zeros( (len(cells),len(cells)), np.float64)
    for ai,a in enumerate(cells):
        a_dists=one_cell_dists(a)
        for bi,b in enumerate(cells[ai+1:],start=ai+1):
            result[ai,bi]=result[bi,ai]=a_dists[b]
    return result


## get polygons for particle exit:

if 0:
    t,parts=ptm_group.read_timestep(1000)

    fig=plt.figure(1)
    fig.clf()
    ax=fig.add_subplot(1,1,1)
    g.plot_edges(ax=ax,lw=0.5,color='k')
    ax.plot(parts['x'][:,0], parts['x'][:,1], 'g.')
    ax.axis('equal')

    from stompy.plot import plot_utils
    res=plot_utils.draw_polyline()


poly_geo=np.array([[ 629664., 4233211.],
                   [ 629823., 4233215.],
                   [ 629821., 4233119.],
                   [ 629661., 4233114.]])
poly_rio=np.array([[ 615061., 4224762.],
                   [ 616046., 4224084.],
                   [ 615810., 4223645.],
                   [ 614768., 4224230.]])
poly_bsp=np.array([[ 605189., 4237133.],
                   [ 605242., 4237153.],
                   [ 605257., 4237118.],
                   [ 605204., 4237098.]])

##

# fast lookups via matplotlib:
from matplotlib import path
ctrs=g.cells_centroid()
dead_cells=np.zeros(g.Ncells(),np.bool8)
for poly in [poly_geo,poly_rio]:
    dead_cells |= path.Path(poly).contains_points(ctrs)
g.plot_cells(mask=dead_cells,color='m',ax=ax)

##
#time_steps=range(len(part_ts))
time_steps=range(0,150,2)
# timeline same for all
track_times=part_ts[time_steps]
track_time_s=(track_times-track_times[0])/np.timedelta64(1,'s')

# truncate time series when particles leave domain

def track_cloud(cloud_particles):
    # track cloud.  this ti is an index into the time steps already extracted
    track_vars=[]

    for ti in utils.progress(time_steps,1.0,"processing %s timesteps"):
        cloud_xy=particles[ti,cloud_particles,:]
        cloud_cells=np.array( [g.select_cells_nearest(xy) for xy in cloud_xy] )
        if np.any( dead_cells[cloud_cells] ):
            # print("Particles hit boundary - truncate time")
            break
        pair_dists=pairwise_grid_distance(cloud_cells)
        variance_per_cell=(pair_dists**2).mean(axis=1)
        best_variance=variance_per_cell.min()
        track_vars.append(best_variance)

    track_vars=np.array(track_vars)
    t_s=track_time_s[:len(track_vars)]

    # give up if there is less than 12h or 5 data points.
    if len(t_s>5) and (t_s[-1]-t_s[0])>43200:
        mb=np.polyfit(t_s,track_vars,1)
        return mb[0]
    else:
        return np.nan

def dispersion_from_x(x0):
    cloud_particles=particle_cloud_for_point(x0)
    return track_cloud(cloud_particles)


def dispersion_for_cell(c):
    # wrapper which checks to make sure a particle actually started in the requested
    # cell
    cloud_particles=particle_cloud_for_point(cc[c])
    ti=0
    cloud_xy=particles[ti,cloud_particles,:]
    for xy in cloud_xy:
        if c==g.select_cells_nearest(xy):
            break
    else:
        # none of the cloud actually started in c, so this dispersion
        # "belong" to a different cell
        return np.nan
    return track_cloud(cloud_particles)

##

per_cell_K=-np.ones(g.Ncells())

##

retry=np.isnan(per_cell_K)
per_cell_K[retry]=-1

##

cc=g.cells_center()

missing=np.nonzero(per_cell_K<0)[0]

# random, to evenly fill out the map...
# missing=missing[ np.argsort( np.random.random(len(missing)) ) ]
# ordered, to capitalize on locality
missing=missing[ np.argsort( cc[missing,1] ) ]

for idx in utils.progress(missing):
    #per_cell_K[idx]=dispersion_from_x(cc[idx])
    per_cell_K[idx]=dispersion_for_cell(idx)

##

cell_nbrs=[ [c]+g.cell_to_cells(c) for c in range(g.Ncells())]

# drop negative boundary values
cell_nbrs=[ np.array(nbrs) for nbrs in cell_nbrs]
cell_nbrs=[ nbrs[nbrs>=0] for nbrs in cell_nbrs]

##
# v02b: deal with BSPP exits, too.
#ds=xr.open_dataset('dispersion_K_v02b.nc')
#per_cell_K=ds.K.values.copy()

smooth_K=per_cell_K.clip(0,np.inf)

for it in range(3):
    print(it)
    new_K=smooth_K.copy()
    for c,nbrs in enumerate(cell_nbrs):
        new_K[c]=np.nanmean( smooth_K[nbrs] )
    smooth_K=new_K

##

if 1: # save to disk along with geometry
    ds=g.write_to_xarray()
    ds['K']=('face',),per_cell_K
    ds['Ksmooth']=('face',),smooth_K
    ds.to_netcdf('dispersion_K_v02b_wsmooth.nc')
    ds.close()

##
#zoom=(601804.8951374379, 633292.0594476878, 4222037.326571205, 4252542.492583467)
zoom=(603276.0170315901, 623336.1938701726, 4230666.761713925, 4248493.5760155255)
plt.figure(3).clf()
fig,(ax1,ax2)=plt.subplots(1,2,sharex=True,sharey=True,num=3)

valid1=np.isfinite(per_cell_K) & (per_cell_K>=0)

coll1=g.plot_cells(values=per_cell_K.clip(1,np.inf),
                   norm=colors.LogNorm(),
                   mask=valid1,ax=ax1)

valid2=np.isfinite(smooth_K)&(smooth_K>=0)
coll2=g.plot_cells(values=smooth_K.clip(1,np.inf),
                   norm=colors.LogNorm(),
                   mask=valid2,ax=ax2)
for coll in [coll1,coll2]:
    coll.set_clim([3,300])
    coll.set_cmap('inferno_r')

for coll,ax in zip([coll1,coll2],[ax1,ax2]):
    plot_wkb.plot_polygon(grid_poly,ax=ax,fc='none',ec='k',lw=0.5)
    plt.colorbar(coll, label="log10 K",ax=ax)
    ax.axis('equal')
    ax.axis(zoom)
    ax.xaxis.set_visible(0)
    ax.yaxis.set_visible(0)

ax1.set_title('Output')
ax2.set_title('Smoothed')
fig.tight_layout()

fig.savefig('third-cutb-dispersion.png',dpi=200)

##

ax.plot(particles[0,:,0],
        particles[0,:,1],
        'm.',ms=4)

##
#zoom=(601804.8951374379, 633292.0594476878, 4222037.326571205, 4252542.492583467)
zoom=(603276.0170315901, 623336.1938701726, 4230666.761713925, 4248493.5760155255)
plt.ioff()

fig=plt.figure(3)
fig.clf()
fig.set_size_inches((7,10),forward=True)
ax=fig.add_axes([0,0,1,1])

valid=np.isfinite(smooth_K)&(smooth_K>=0)
coll=g.plot_cells(values=smooth_K.clip(1,np.inf),
                  norm=colors.LogNorm(),
                  mask=valid,ax=ax)
coll.set_clim([3,300])
coll.set_cmap('inferno_r')

plot_wkb.plot_polygon(grid_poly,ax=ax,fc='none',ec='k',lw=0.5)
cax=fig.add_axes([0.05,0.65,0.02,0.3])
plt.colorbar(coll, label="log10 K (m$^2$ s$^{-1}$)",cax=cax)
ax.axis('equal')
ax.axis(zoom)
ax.xaxis.set_visible(0)
ax.yaxis.set_visible(0)

fig.savefig('dispersion-map-2014-04.png',dpi=250)
plt.ion()


# depth is an issue, as particles which move into deeper water
# converge in the horizontal, which is advection in reality, but
# anti-dispersion in the math.


