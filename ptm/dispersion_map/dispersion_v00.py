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
grid_poly=g.boundary_polygon()

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
t,parts=ptm_group.read_timestep(0)

##

# Select a small number of particles and watch their spread.

def particles_within_r(parts,x0,r):
    deltas=parts['x'][:,:2] - x0
    dist=(deltas**2).sum(axis=1)
    sel=dist<r**2
    return np.nonzero(sel)[0]

##
x0=np.array([612129., 4238225.])
cloud_idxs=particles_within_r(parts,x0,r=300)

##

#zoom=(597913.7274775933, 648118.8262812896, 4217179.54644355, 4301202.344200377)
zoom=(610351.7663721164, 614143.730357768, 4236937.159272692, 4240093.599179518)
plt.figure(1).clf()
fig,ax=plt.subplots(num=1)

ax.plot(parts['x'][:,0],
        parts['x'][:,1],'g.')
ax.plot(parts['x'][cloud_idxs,0],
        parts['x'][cloud_idxs,1],'m.',ms=10)

plot_wkb.plot_polygon(grid_poly,ax=ax,fc='0.8',ec='k',lw=0.5)

ax.axis('equal')
ax.axis(zoom)

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

cloud_t_xy = part_locs[:,cloud_idxs,:]

cloud_cc = cloud_t_xy.mean(axis=1)

# In this case it remains more or less in the channel.
ax.plot(cloud_cc[:,0],cloud_cc[:,1],'k-')

spread=cloud_t_xy - cloud_cc[:,None,:] # [time,particle,xy]
cloud_var=(spread**2).sum(axis=2).mean(axis=1)

##

part_seconds=(part_ts-part_ts[0])/np.timedelta64(1,'s')
mb=np.polyfit(part_seconds,cloud_var,1)

plt.figure(2).clf()
plt.plot(part_ts,cloud_var)
plt.plot(part_ts,np.polyval(mb,part_seconds),'k-',lw=0.5,
         label='%.2f m$^2$/s'%mb[0])
plt.legend()


##

def within_r(part_locs,x0,r,at_least=0):
    deltas=part_locs[:,:] - x0
    dist=(deltas**2).sum(axis=1)
    sel=dist<r**2
    idxs=np.nonzero(sel)[0]
    if len(idxs)<at_least:
        ordering=np.argsort(dist)
        idxs=ordering[:at_least]
    return idxs

def dispersion_from_x(x0,particles,r=300,min_particles=10):
    cloud_idxs=within_r(particles[0,:,:],x0,r=r,at_least=min_particles)

    if len(cloud_idxs)<min_particles:
        return np.nan

    cloud_t_xy = particles[:,cloud_idxs,:]
    cloud_cc = cloud_t_xy.mean(axis=1)
    spread=cloud_t_xy - cloud_cc[:,None,:] # [time,particle,xy]
    cloud_var=(spread**2).sum(axis=2).mean(axis=1)
    mb=np.polyfit(part_seconds,cloud_var,1)
    return mb[0]


# 3ms per call?
# d=dispersion_from_x(x0,part_locs)
##

per_cell_K=-np.ones(g.Ncells())

##

retry=np.isnan(per_cell_K)
per_cell_K[retry]=-1

##
cc=g.cells_center()

missing=np.nonzero(per_cell_K<0)[0]
missing=missing[ np.argsort( np.random.random(len(missing)) ) ]

for idx in utils.progress(missing):
    per_cell_K[idx]=dispersion_from_x(cc[idx],part_locs,min_particles=5)

##

#zoom=(610351.7663721164, 614143.730357768, 4236937.159272692, 4240093.599179518)
zoom=(601804.8951374379, 633292.0594476878, 4222037.326571205, 4252542.492583467)
plt.figure(3).clf()
fig,ax=plt.subplots(num=3)

valid=np.isfinite(per_cell_K) & (per_cell_K>=0)
coll=g.plot_cells(values=per_cell_K.clip(1,np.inf),
                  norm=colors.LogNorm(),
                  mask=valid)

coll.set_clim([1,100])
coll.set_cmap('jet')
plot_wkb.plot_polygon(grid_poly,ax=ax,fc='none',ec='k',lw=0.5)

plt.colorbar(coll, label="log10 K")
ax.axis('equal')
ax.axis(zoom)
fig.savefig('first-cut-dispersion.png')
##

# need to detect when particles exit!

