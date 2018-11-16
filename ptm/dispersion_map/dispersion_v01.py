"""
Another shot at dispersion, but using a distance metric restricted to be within the
grid.
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
@memoize.memoize(lru=15000)
def one_cell_dists(c):
    return csgraph.shortest_path(graph,directed=False,indices=[c])[0,:]


def pairwise_grid_distance(cells):
    result=np.zeros( (len(cells),len(cells)), np.float64)
    for ai,a in enumerate(cells):
        a_dists=one_cell_dists(a)
        for bi,b in enumerate(cells[ai+1:],start=ai+1):
            result[ai,bi]=result[bi,ai]=a_dists[b]
    return result


##
#time_steps=range(len(part_ts))
time_steps=range(0,150,2)
# timeline same for all
track_times=part_ts[time_steps]
track_time_s=(track_times-track_times[0])/np.timedelta64(1,'s')


def track_cloud(cloud_particles):
    # track cloud.  this ti is an index into the time steps already extracted
    track_vars=[]

    for ti in utils.progress(time_steps,1.0,"processing %s timesteps"):
        cloud_xy=particles[ti,cloud_particles,:]
        cloud_cells=np.array( [g.select_cells_nearest(xy) for xy in cloud_xy] )
        pair_dists=pairwise_grid_distance(cloud_cells)
        variance_per_cell=(pair_dists**2).mean(axis=1)
        best_variance=variance_per_cell.min()
        track_vars.append(best_variance)

    track_vars=np.array(track_vars)

    mb=np.polyfit(track_time_s,track_vars,1)
    return mb[0]

def dispersion_from_x(x0):
    cloud_particles=particle_cloud_for_point(x0)
    return track_cloud(cloud_particles)

##

# fig=plt.figure(1)
# fig.clf()
# ax=fig.add_subplot(1,1,1)
# plot_wkb.plot_wkb(grid_poly,zorder=-1,fc='0.9',ec='k',lw=0.4)
# ax.axis('equal')
# 
# ti=np.arange(0,len(particles),10)
# cloud_t_xy=particles[ti][:,cloud_particles,:]
# cloud_times=ti[:,None]*np.ones_like(cloud_t_xy[...,0])
# ax.scatter(cloud_t_xy[...,0], cloud_t_xy[...,1],
#            20,cloud_times,cmap='jet')
# 
# ## and a time series
# fig2=plt.figure(2)
# fig2.clf()
# ax2=fig2.add_subplot(1,1,1)
# track_time_s=(track_times-track_times[0])/np.timedelta64(1,'s')
# mb=np.polyfit(track_time_s,track_vars,1)
# 
# ax2.plot(track_times,track_vars,'g-',label='Particles')
# ax2.plot(track_times,np.polyval(mb,track_time_s),'k-',lw=0.5,label="K=%.2f m2 s-1"%(mb[0]))
# ax2.legend()

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
    per_cell_K[idx]=dispersion_from_x(cc[idx])

##
if 0: # save to disk along with geometry
    ds=g.write_to_xarray()
    ds['K']=('face',),per_cell_K
    ds.to_netcdf('dispersion_K_v01.nc')
    ds.close()

##
#zoom=(601804.8951374379, 633292.0594476878, 4222037.326571205, 4252542.492583467)
zoom=(603276.0170315901, 623336.1938701726, 4230666.761713925, 4248493.5760155255)
plt.figure(3).clf()
fig,ax=plt.subplots(num=3)

valid=np.isfinite(per_cell_K) & (per_cell_K>=0)
coll=g.plot_cells(values=per_cell_K.clip(1,np.inf),
                  norm=colors.LogNorm(),
                  mask=valid)

coll.set_clim([1,300])
coll.set_cmap('jet')
plot_wkb.plot_polygon(grid_poly,ax=ax,fc='none',ec='k',lw=0.5)

plt.colorbar(coll, label="log10 K")
ax.axis('equal')
ax.axis(zoom)

fig.savefig('second-cut-dispersion.png')

##

ax.plot(particles[0,:,0],
        particles[0,:,1],
        'm.',ms=4)

##

# need to detect when particles exit!
# depth is also an issue, as particles which move into deeper water
# converge in the horizontal, which is advection in reality, but
# anti-dispersion in the math.


