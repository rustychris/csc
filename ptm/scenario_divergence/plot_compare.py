"""
Compare tracks from two PTM runs
"""
from stompy.grid import unstructured_grid
from stompy.model.fish_ptm import ptm_tools
from matplotlib import collections
from stompy import utils
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
##

g=unstructured_grid.UnstructuredGrid.from_ugrid('../../dflowfm/runs/20180807_grid98_16_single/ptm_hydro.nc')

##

ds_map=xr.open_dataset('../../dflowfm/runs/20180807_grid98_16_single/DFM_OUTPUT_flowfm/flowfm_map.nc')
def particle_water_depth(x,time):
    tidx=np.searchsorted(ds_map.time,time)
    cell_depths=ds_map.mesh2d_waterdepth.isel(time=tidx).values
    x_cells=[g.select_cells_nearest(xi) for xi in x]
    return cell_depths[x_cells]

##
dist_bspp=ptm_tools.PtmBin('nobspp/DIST_bin.out')
dist_nobspp=ptm_tools.PtmBin('bspp/DIST_bin.out')

##
ntimes=min(dist_bspp.count_timesteps(),dist_nobspp.count_timesteps())

##
zoom_lindsey=(603641.9059474986, 610614.5678962235, 4233101.860312233, 4237901.354960644)
zoom=zoom_lindsey

plt.figure(1).clf()
fig,ax=plt.subplots(num=1)

g.plot_edges(color='k',lw=0.4,ax=ax,clip=zoom)

ax.axis(zoom)
cax=fig.add_axes([0.93,0.2,0.02,0.6])

##
ntimes=100

def never_stuck(ptm_out,ntimes):
    t_a,state = ptm_out.read_timestep(ts=0)
    last_x=state['x'][:,:2]
    stuck_count=np.zeros(len(last_x),np.int32)

    for ti in utils.progress(range(1,ntimes)):
        t_a,state = ptm_out.read_timestep(ts=ti)
        this_x=state['x'][:,:2]
        stuck=np.all(last_x==this_x,axis=1)
        stuck_count[ stuck ] += 1
        last_x=this_x
    return stuck_count==0

if 1:
    valid_a=never_stuck(dist_bspp,ntimes)
    valid_b=never_stuck(dist_nobspp,ntimes)
    valid_never_stuck=valid_a&valid_b # a bit less than half.

def update(ti):
    del ax.collections[1:]
    del ax.lines[:]
    t_a,parts_a = dist_bspp.read_timestep(ts=ti)
    t_b,parts_b = dist_nobspp.read_timestep(ts=ti)

    segs= np.concatenate( [ parts_a['x'][:,None,:2],
                            parts_b['x'][:,None,:2] ],
                          axis=1)

    mode='dist_at_a'
    validate='never_stuck'

    valid=np.ones(len(parts_a),np.bool8)

    if validate=='depth':
        time=utils.to_dt64(dist_bspp.time[ti])
        depths_a=particle_water_depth(parts_a['x'][:,:2],time)
        depths_b=particle_water_depth(parts_b['x'][:,:2],time)
        valid=valid&(depths_a>0.05)&(depths_b>=0.05)
    if validate=='never_stuck':
        valid=valid&valid_never_stuck

    del ax.lines[:]

    segs=segs[valid,:,:]

    items=[]
    if mode=='segs':
        lcoll=collections.LineCollection(segs,color='m',lw=1.0)
        ax.add_collection(lcoll)
        items.append(lcoll)
    elif mode=='dist_at_a':
        separations=utils.dist( segs[:,0,:], segs[:,1,:] )
        order=np.argsort(separations)
        scat=ax.scatter( segs[order,0,0], segs[order,0,1], 20, separations[order],
                         cmap='jet')
        items.append(scat)
    elif mode=='dist_at_b':
        separations=utils.dist( segs[:,0,:], segs[:,1,:] )
        order=np.argsort(separations)
        scat=ax.scatter( segs[order,1,0], segs[order,1,1], 20, separations[order],
                         cmap='jet')
        items.append(scat)

    if mode!='segs':
        cax.cla()
        plt.colorbar(items[0],cax=cax)
        items[0].set_clim([0,500])

    return items

##
ti=500
update(ti)

##
for ti in range(0,500,5):
    update(ti)
    plt.draw()
    plt.pause(0.001)

##

# Getting a lot of beached particles.
# filtering out each time step based on local depth not that useful.

##


