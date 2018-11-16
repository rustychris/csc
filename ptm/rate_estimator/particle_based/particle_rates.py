import matplotlib.pyplot as plt
import numpy as np
import heapq

from shapely.ops import cascaded_union

from stompy import utils
from stompy.grid import unstructured_grid
from stompy.model.fish_ptm import ptm_tools
from stompy.plot import plot_wkb
from stompy.io.local import usgs_nwis
from stompy.spatial import proj_utils

##
cache_dir="cache"
os.path.exists(cache_dir) or os.mkdir(cache_dir)

##
g=unstructured_grid.UnstructuredGrid.from_ugrid('../../../dflowfm/runs/20180807_grid98_17/ptm_hydro.nc')
grid_poly=g.boundary_polygon()

ll2utm=proj_utils.mapper('WGS84','EPSG:26910')

##

ptm_run_dir="../run_10days"
ptm_groups=["INIT","SAC","SRV"]

ptm_data=[ ptm_tools.PtmBin(os.path.join(ptm_run_dir,grp+"_bin.out"))
           for grp in ptm_groups ]

ntimes=ptm_data[0].count_timesteps()

run_start=ptm_data[0].read_timestep(0)[0]
run_stop =ptm_data[0].read_timestep(ntimes-1)[0]

##

colors=['g','b','m','r','k']

zoom=(597913.7274775933, 648118.8262812896, 4217179.54644355, 4301202.344200377)
plt.figure(1).clf()
fig,ax=plt.subplots(num=1)

for d,col in zip(ptm_data,colors):
    d.plot(ntimes-1,ax=ax,zoom=zoom,update=False,ms=1,color=col)

plot_wkb.plot_polygon(grid_poly,ax=ax,fc='0.8',ec='k',lw=0.5)

ax.axis(zoom)

##

# Define the stations time series
# USGS stations:
usgs_station_codes=dict(
    fpt=11447650,
    sdc=11447890,
    ges=11447905,
    srv=11455420, # will have to get nudged a bit
    # ryi=11455350, # no data in this particular period
    hwb=11455165, # miner slough at hwy 84
    dws=11455335, # deep water ship channel, near RV
    ame=11446980, # American River, may have to get nudged
    ver=11425500,
    cli=11455315, # Cache Slough at S. Liberty Island
    cah=11455280, # Cache Slough near Hastings Tract
    ula=11455268, # Ulatis near Elmira
)

stations=[]
for stn_name in usgs_station_codes:
    stn_code=usgs_station_codes[stn_name]

    ds=usgs_nwis.nwis_dataset(stn_code,run_start,run_stop,products=[10],
                              cache_dir=cache_dir,cache_only=True)
    if ds is None:
        print("No data for %s -- skip."%stn_code)
        continue
    meta=usgs_nwis.station_metadata(stn_code,cache_dir=cache_dir)
    ds['latitude']=(),meta['lat']
    ds['longitude']=(),meta['lon']
    ll=np.array( [ds.longitude.values,ds.latitude.values])
    ds['ll']=('xy',),ll
    xy=ll2utm(ll)
    ds['xy']=('xy',),xy
    ds.attrs['name']=stn_name
    stations.append(ds)

# Based on the locations, the American River Station is maybe too far
# away, and since the model doesn't have American River flows, it
# doesn't make sense.
print("Omitting American River")
stations=[s for s in stations if s.site_no!=str(usgs_station_codes['ame'])]

station_xys=np.array([ station.xy.values
                       for station in stations ])

##

# Plot timeseries for each of those
fig=plt.figure(2)
fig.clf()
ax=fig.add_subplot(1,1,1)

for stn in stations:
    ax.plot(stn.time,stn.temperature_water,label=str(stn.attrs['name']).upper(),
            lw=1.0)

ax.axis(xmin=stn.time.values[0],xmax=stn.time.values[-1])

ax.legend(ncol=4,loc='lower right')
ax.set_ylabel('Water Temp $^\circ$C')
fig.autofmt_xdate()
fig.savefig('station-temp-timeseries.png',dpi=150)

##

if 0:
    colors=['g','b','m','r','k']

    zoom=(597913.7274775933, 648118.8262812896, 4217179.54644355, 4301202.344200377)
    fig=plt.figure(3)
    fig.clf()
    ax=fig.add_subplot(1,1,1)

    for d,col in zip(ptm_data,colors):
        d.plot(ntimes-1,ax=ax,zoom=zoom,update=False,ms=1,color=col)

    plot_wkb.plot_polygon(poly,ax=ax,fc='0.8',ec='k',lw=0.5)

    ax.plot(station_xys[:,0],station_xys[:,1],'ko',ms=5)

    ax.axis(zoom)

##

# Define footprints for each station
Ac=g.cells_area()
cc=g.cells_center()

for station in stations:
    xy=station.xy.values
    c0=g.select_cells_nearest(xy)

    station['cell0']=(),c0

    # arbitrary for now
    target_A=400*400

    def cost(c):
        return utils.dist(xy,cc[c])
    heap=[(cost(c0),c0)]
    nearby=[]

    # Note that this does not consider the case where regions
    # overlap between stations
    sum_area=0
    while heap and sum_area<target_A:
        _,c=heapq.heappop(heap)
        sum_area+=Ac[c]
        nearby.append(c)
        for nbr in g.cell_to_cells(c):
            if nbr in nearby: continue
            heapq.heappush(heap, (cost(nbr),nbr))

    #g.plot_cells(mask=[c0],color='r',ax=ax,zorder=10)
    #g.plot_cells(mask=nearby,ax=ax,zorder=9,color='b')

    station['cells']=('cells',),nearby

    # make that into a polygon for maybe quicker tests
    poly=cascaded_union( [ g.cell_polygon(c) for c in nearby] )
    station['footprint']=(), poly

##
import time

# This is the super slow step
# but it probably shouldn't be slow.
def group_trajectories():
    trajectories=dict() # (grp,idx in group) => [ (t,xy), ..]

    times=range(0,ntimes)

    start_time=time.time()

    # by avoiding keeping references to particle arrays, and
    # sorting stuck particles on the fly, the memory usage is
    # kept in check
    for part_grp in range(len(ptm_groups)):
        print("Particle Group %d"%part_grp)

        # read the last time step to get the max particle
        # count, and final resting place.
        step_t,parts=ptm_data[part_grp].read_timestep(ntimes-1)
        end_xy=parts['x']
        # we output particles once when stuck, but then record
        # them here to avoid any more
        dead=np.zeros(len(end_xy),np.bool8)
        max_particles=len(parts)

        # use list to speed up references
        grp_trajectories=[ [] for i in range(max_particles) ]

        last_xy=np.nan * np.ones((max_particles,2))
        lastlast_xy=np.nan * np.ones((max_particles,2))

        for ti,t in enumerate(utils.progress(times)):
            step_t,parts=ptm_data[part_grp].read_timestep(t)
            step_t=utils.to_dt64(step_t)
            Nstep=len(parts)
            # move as much referencing outside the loop
            part_xy=parts['x'][:,:]

            # particles which do not move between this step and
            # the end of the simulation are presumed dead and
            # not procesed anymore.
            # it might be worth keeping the first dead location.
            # this is probably stripping out some otherwise useful
            # points.
            for part_idx in np.nonzero(~dead[:Nstep])[0]: # range(Nstep):
                traj=grp_trajectories[part_idx]
                # avoid any references back to parts
                # assumes that traj[-1] is same location as traj[-2]
                # probably safe.
                rec=[step_t,
                     part_xy[part_idx,0],
                     part_xy[part_idx,1]]
                if len(traj)>=2 and (traj[-2][1]==rec[1]) and (traj[-2][2]==rec[2]):
                    # stuck particles are just updated by the latest time/position
                    traj[-1][0]=step_t
                else:
                    traj.append(rec)
            # if a particle is stuck from here on out, remove it from play
            dead[:Nstep]=np.all(part_xy==end_xy[:Nstep,:],axis=1)

        for part_idx,traj in enumerate(grp_trajectories):
            trajectories[ (part_grp,part_idx) ] = traj

    return trajectories

trajectories=group_trajectories()

##

trimmed=dict()

# Trim each trajectory, and convert to numpy array
for k in utils.progress(trajectories.keys()):
    traj=trajectories[k]
    if len(traj)>2:
        t,x,y = zip(*traj)
        traj_np=np.zeros(len(traj),dtype=[('t','M8[us]'),('x',np.float64),('y',np.float64)])
        traj_np['t']=t
        traj_np['x']=x
        traj_np['y']=y
        trimmed[k]=traj_np

trimmed_l=list(trimmed.values())

##

# Detect when those particles pass through a footprint
import matplotlib.path as mpltPath

paths=[mpltPath.Path( np.array(station['footprint'].values.item().exterior) )
       for station in stations]

trimmed_with_hits=[] # [ (trajectory, hit_indexes), ...]
# make a single no-hit record, which will be copied for anybody
# who actually gets hits

hit_recs0=np.zeros(ntimes,np.int16)-1

for traj in utils.progress(trimmed_l,msg="%s of trimmed trajectories"):
    pnts=np.c_[traj['x'], traj['y']]
    my_hits=hit_recs0[:len(pnts)] # copy before changing!
    for path_i,path in enumerate(paths):
        hits=path.contains_points(pnts)
        if hits.sum():
            if my_hits.base is hit_recs0:
                my_hits=my_hits.copy()
            my_hits[hits]=path_i
    trimmed_with_hits.append( [traj,my_hits] )
##

# How many trajectories intersect a station?
# 96%!  not bad. maybe a lot of those are exiting
# via SRV, tho.
count=0
for traj,hits in trimmed_with_hits:
    if hits.base is not hit_recs0:
        count+=1

##

station_data=[]

# preprocess that data a bit
for station in stations:
    # Apply the UTC => PST conversion here
    dnums=utils.to_dnum(station.time.values) - 8/24.
    temps=station.temperature_water.values
    station_data.append( np.c_[dnums,temps] )

##
# Fill in temperature for when particles encounter a station
def observed(t,station_i):
    return np.interp( utils.to_dnum(t),
                      station_data[station_i][:,0],
                      station_data[station_i][:,1] )

for traj_i in utils.progress(range(len(trimmed_with_hits))):
    rec=trimmed_with_hits[traj_i]
    traj=rec[0]
    hits=rec[1]
    rec=rec[:2] # idempotency!
    T=np.zeros( len(traj), np.float32 )
    T[:]=np.nan
    good = np.nonzero(hits>=0)[0]
    if len(good):
        for idx in good:
            t=traj['t'][idx]
            T[idx]=observed(t,hits[idx])
    rec.append(T)
    trimmed_with_hits[traj_i]=rec

##

# linear interpolate in time along each trajectory
#
traj_data=[]

missing=np.nan*np.ones(ntimes)

for traj_i in utils.progress(range(len(trimmed_with_hits))):
    traj,hits,T = trimmed_with_hits[traj_i]
    rec=dict(traj=traj,
             hits=hits,
             T=T)
    # Additional data:

    # some trajectories have no hits,
    # some with only a single hit or cluster of hits
    measured=np.nonzero(hits>=0)[0]

    # This will remain all nan if there is no valid data in T
    rec['Tfill']=utils.fill_invalid(T)
    if len(measured)==0:
        # Save on memory when there is no data.
        rec['fill_dist']=missing[:len(T)]
        rec['dTdt']=missing[:len(T)]
    else:
        # trajectory time in seconds
        t_s=(traj['t']-traj['t'][0])/np.timedelta64(1,'s')
        # just the times that were observed
        t_measured=t_s[measured]
        # second between each point and the nearest observation
        fill_dist=t_s - utils.nearest_val(t_measured,t_s)
        # that's signed seconds offset from the nearest in time measurement
        rec['fill_dist']=fill_dist

        if len(measured)>1:
            rates=np.nan*T
            for a,b in zip(measured[:-1],
                           measured[1:]):
                dt=(t_s[b]-t_s[a])
                # inclusive index, to avoid nan at the stations
                rates[a:b+1]=(T[b]-T[a])/dt
            rec['dTdt']=rates
        else:
            rec['dTdt']=missing[:len(T)]

    traj_data.append(rec)

##

def extract_particle_snapshot(t):
    # extract particle field at a given time
    snapshot=[]
    for traj_i,traj in enumerate(utils.progress(traj_data)):
        traj_t=traj['traj']['t']
        if traj_t[0]>t:
            continue
        if traj_t[-1]<t:
            continue
        # time index along this trajectory matching up with t
        idx=np.searchsorted(traj_t,t)
        # x,y,Tfill,dTdt,fill_dist
        rec=[traj['traj']['x'][idx],
             traj['traj']['y'][idx],
             traj['Tfill'][idx],
             traj['dTdt'][idx],
             traj['fill_dist'][idx],
             traj_i]
        snapshot.append(rec)

    snap_struct=np.array(snapshot)
    snap_data=dict( x=snap_struct[:,0],
                    y=snap_struct[:,1],
                    Tfill=snap_struct[:,2],
                    dTdt=snap_struct[:,3],
                    fill_dist=snap_struct[:,4],
                    traj_i=snap_struct[:,5].astype(np.int32) )
    return snap_data

##

class ParticlePlot(object):
    zoom=(601230, 634344, 4223225, 4246485)
    ax=None
    # if particle has not been observed within this time, do not display
    # (seconds)
    stale_thresh_s=np.inf

    def __init__(self,fig,t,**kw):
        self.__dict__.update(kw)

        snap_data=extract_particle_snapshot(t)
        if self.ax is None:
            fig.clf()
            ax=fig.add_subplot(1,1,1)
            ax.set_position([0,0,1,1])
        else:
            ax=self.ax

        plot_wkb.plot_polygon(grid_poly,ax=ax,fc='0.8',ec='none',lw=0.5,zorder=-2)
        plot_wkb.plot_polygon(grid_poly,ax=ax,fc='none',ec='k',lw=0.5,zorder=2)

        station_xys=np.array([ station.xy.values
                               for station in stations ])
        ax.plot(station_xys[:,0],station_xys[:,1],'ko',ms=5)

        values=self.snap_to_value(snap_data)

        sel=np.isfinite(values)
        stale_s=np.abs(snap_data['fill_dist'])
        sel=sel&(stale_s<self.stale_thresh_s)
        scat=ax.scatter(snap_data['x'][sel],
                        snap_data['y'][sel],
                        40*( 1/(1+stale_s[sel]/1800.) ),
                        values[sel])
        scat.set_cmap(self.cmap)
        scat.set_clim(self.clim)
        pos=ax.get_position()
        cax=fig.add_axes([pos.xmin+0.08,
                          pos.ymin+0.93*pos.height,
                          pos.width*0.25,
                          0.02])
        ax.text(0.08,0.96,utils.to_datetime(t).strftime("%Y-%m-%d %H:%M PST"),
                transform=ax.transAxes)

        plt.colorbar(scat,cax=cax,orientation='horizontal',
                     label=self.units_label)

        ax.xaxis.set_visible(0)
        ax.yaxis.set_visible(0)

        ax.axis('equal')
        ax.axis(self.zoom)

class ParticlePlotTemperature(ParticlePlot):
    clim=[12,20]
    cmap='jet'
    units_label=r'$^{\circ}$C'
    def snap_to_value(self,snap_data):
        return snap_data['Tfill']

class ParticlePlotdTdt(ParticlePlot):
    clim=[-1,1]
    cmap='seismic'
    units_label=r'$^{\circ}$C/day'
    def snap_to_value(self,snap_data):
        return snap_data['dTdt'] * 86400

##


fig=plt.figure(4)
fig.set_size_inches([12.5,8.75])

times=np.arange(np.datetime64("2014-04-05 00:00"),
                np.datetime64("2014-04-28 00:00"),
                np.timedelta64(900,"s"))

if 0: # temperature values
    for ti,t in enumerate(times):
        img_fn='frames2/temperature-particles-frame%04d.png'%ti
        #if os.path.exists(img_fn):
        #    continue
        print(img_fn)

        ParticlePlotTemperature(fig,t)
        break
        fig.savefig(img_fn)
if 0:
    for ti,t in enumerate(times):
        img_fn='frames3/rates-particles-frame%04d.png'%ti
        #if os.path.exists(img_fn):
        #    continue
        print(img_fn)

        ParticlePlotdTdt(fig,t)
        fig.savefig(img_fn)

##
##
# Take a closer look at one particle, at
# numpy.datetime64('2014-04-09T08:45:00')
# tgt=np.array( [615285.1632649068, 4224287.880627088] )
tgt=np.array([617599.4202082362, 4228851.501966425])

dx=snap_data['x'] - tgt[0]
dy=snap_data['y'] - tgt[1]

dist=(dx**2+dy**2)
best=np.argmin(dist) # 14195

ax.plot( [snap_data['x'][best]],
         [snap_data['y'][best]],
         'go',ms=10)

##

# take a closer look at the time series for one
# particle: 14195 from this snapshot.
# traj_i for this particle is 24347
traj_i=int(snap_data['traj_i'][best])

traj=traj_data[traj_i]

##

ax.plot( traj['traj']['x'],
         traj['traj']['y'],
         'g-')

##

class PlotJourney(object):
    zoom=(601230, 634344, 4223225, 4246485)

    ax_map=None
    ax_ts=None
    ax_ts2=None
    clim=[12,17]
    cmap='jet'
    def __init__(self,fig,traj,**kw):
        self.__dict__.update(kw)
        if not fig.axes:
            ax_map=fig.add_subplot(1,2,1)
            ax_ts=fig.add_subplot(2,2,2)
            ax_ts2=fig.add_subplot(2,2,4)
        else:
            # assume that ax_map, ax_ts, and/or ax_ts2 were passed in
            ax_map=self.ax_map
            ax_ts=self.ax_ts
            ax_ts2=self.ax_ts2

        if ax_map:
            pos=ax_map.get_position()

            plot_wkb.plot_polygon(grid_poly,ax=ax_map,fc='0.8',ec='0.6',lw=0.5,zorder=-2)
            ax_map.axis('equal')
            ax_map.axis( (606812, 633363, 4222818, 4256651) )
            ax_map.xaxis.set_visible(0)
            ax_map.yaxis.set_visible(0)

            ax_map.plot( traj['traj']['x'],
                         traj['traj']['y'],
                         'k-',lw=0.5)
            sel=np.zeros(len(traj['Tfill']),np.bool8)
            sel[:]=True # thinned out
            obs=traj['fill_dist']==0
            sel[obs]=True
            sizes=5*np.ones_like(sel)
            sizes[obs]=60
            scat=ax_map.scatter( traj['traj']['x'][sel],
                                 traj['traj']['y'][sel],
                                 sizes[sel],traj['Tfill'][sel],
                                 cmap=self.cmap)
            scat.set_clim(self.clim)

        if ax_ts:
            ax_ts.plot( traj['traj']['t'], traj['Tfill'], label='Tfill')
            ax_ts.set_ylabel('Lagrangian temp $^{\circ}$C')
        if ax_ts2:
            ax_ts2.plot(traj['traj']['t'], np.abs(traj['fill_dist'])/3600., label='fill_dist')
            ax_ts2.set_ylabel('Hours from observation')
            if ax_ts:
                ax_ts2.set_xlim( ax_ts.get_xlim() )

        fig.autofmt_xdate()
        if ax_map:
            cax=fig.add_axes( [pos.xmin,pos.ymin-0.08,pos.width,0.03] )
            plt.colorbar(scat,cax=cax,orientation='horizontal',label='$^{\circ}$C')

if 0:
    fig=plt.figure(5)
    fig.clf()

    fig.savefig('sample-trajectory.png',dpi=150)


##

fig=plt.figure(6)
fig.clf()
fig.set_size_inches([10.5,7.25])
t=times[140]
ax1=fig.add_subplot(1,2,2)
ax2=fig.add_subplot(1,2,1)
fig.subplots_adjust(bottom=0.15,left=0.05,right=0.95,top=0.97)

ppt=ParticlePlotTemperature(fig,t,ax=ax1,stale_thresh_s=2*86400)
pj=PlotJourney(fig,traj,ax_map=ax2,clim=ppt.clim)
ax1.axis(ax2.axis())

fig.axes[2].set_visible(0) # hide second colorbar

fig.savefig('trajectory_and_field.png',dpi=125)
