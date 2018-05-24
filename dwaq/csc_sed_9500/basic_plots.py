import six

import matplotlib.pyplot as plt
import xarray as xr 
import stompy.model.delft.io as dio
import stompy.model.delft.waq_scenario as waq
from stompy.model.delft import dfm_grid
from stompy.plot import plot_wkb

from stompy.grid import unstructured_grid

##

six.moves.reload_module(unstructured_grid)
six.moves.reload_module(dfm_grid)
six.moves.reload_module(waq)
six.moves.reload_module(dio)


##

map_ds,g=dio.read_map("csc_sed_9500/tbd.map",return_grid=True)

## 

ssc=map_ds.IM1.isel(time=-1,layer=0).values.copy()

ssc[ssc==-999]=np.nan


## 
fig=plt.figure(1)
fig.clf()
ax=fig.add_axes([0,0,1,1])
# g.plot_edges(lw=0.5,color='k',ax=ax)
ccoll2=g.plot_cells(lw=0.5,color='0.9',ax=ax,mask=np.isnan(ssc))
ccoll=g.plot_cells(lw=0.5,values=ssc,ax=ax,mask=np.isfinite(ssc))
ccoll.set_clim([0,50])
ccoll.set_edgecolor('face')
ax.axis('equal')

##

fig_dir='frames_csc_9500'
os.path.exists(fig_dir) or os.mkdir(fig_dir)

##
boundary=g.boundary_polygon()

##

# Animation of im1 and im2
fig=plt.figure(2)
fig.clf()
ax1=fig.add_axes([0,0,0.5,1])
cax1=fig.add_axes([0.05,0.5,0.03,0.3])
ax2=fig.add_axes([0.5,0,0.5,1])
cax2=fig.add_axes([0.55,0.5,0.03,0.3])

ax2.yaxis.set_visible(0)
ax1.text(0.05,0.95,'IM1',transform=ax1.transAxes)
ax2.text(0.05,0.95,'IM2',transform=ax2.transAxes)

ccolls=[]
plot_vars=['IM1','IM2']
for v,ax,cax in zip( plot_vars,
                     [ax1,ax2],
                     [cax1,cax2] ):
    ssc=map_ds[v].isel(time=-1,layer=0).values
    ssc=np.ma.array(ssc,mask=(ssc==-999)) # map_ds[v].isel(time=-1,layer=0).values.copy()

    ax.set_facecolor('0.6')
    plot_wkb.plot_polygon(boundary,facecolor='0.8',edgecolor='none',ax=ax,lw=0.25,zorder=-1)
    ccoll=g.plot_cells(lw=0.5,values=ssc,ax=ax,cmap='CMRmap_r') 
    ccoll.set_clim([0,50])
    ccoll.set_edgecolor('face')
    ccolls.append(ccoll)
    plt.colorbar(ccoll,cax=cax)
    ax.axis('equal')

##

for frame,ti in enumerate(range(0,800,2)):
    print "[%04d] %d/%d"%(frame,ti,len(map_ds.time))
    for v,ccoll in zip(plot_vars,ccolls):
        ssc=map_ds[v].isel(time=ti,layer=0).values
        ssc=np.ma.array(ssc,mask=(ssc==-999))
        ccoll.set_array(ssc)
    #fig.canvas.draw()
    #plt.pause(0.05)
    fig.set_size_inches( [6.4,4.8],forward=True )
    fig.savefig( os.path.join(fig_dir,'im1_im2_%04d.png'%frame),dpi=140 )

