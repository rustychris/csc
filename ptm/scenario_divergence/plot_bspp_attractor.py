from stompy.grid import unstructured_grid
from stompy.model.fish_ptm import ptm_tools

##

g=unstructured_grid.UnstructuredGrid.from_ugrid('../../dflowfm/runs/20180807_grid98_16_single/ptm_hydro.nc')

##

bspp=ptm_tools.PtmBin('bspp/DIST_bin.out')
bspp_dupe=ptm_tools.PtmBin('bspp/DIST_bin.out')

ntimes=ges.count_timesteps()

##

##

for ti in range(0,ntimes,2):
    bspp.plot(ti,ax=ax,zoom=zoom)
    plt.draw()
    plt.pause(0.001)

##
bspp_intake=np.array([(605206, 4237111.),
                      (605246, 4237130.),
                      (605238, 4237147.),
                      (605195, 4237129.)])

t_a,parts = bspp.read_timestep(ts=ntimes)

final_cells=[g.select_cells_nearest(x) for x in parts['x'][:,:2]]

from shapely import geometry
intake_poly=geometry.Polygon(bspp_intake)

intake_cells=g.select_cells_intersecting(intake_poly)

captured=intake_cells[final_cells]

##

#zoom_lindsey=(603641.9059474986, 610614.5678962235, 4233101.860312233, 4237901.354960644)
zoom_lindsey_out=(604352.5301599499, 611980.3180886208, 4233788.445937092, 4237602.339901428)
zoom=zoom_lindsey_out

fig=plt.figure(1)
fig.clf()
fig.set_size_inches((10,5),forward=True)
ax=fig.add_axes([0,0,1,1])
ax.xaxis.set_visible(0)
ax.yaxis.set_visible(0)

g.plot_edges(color='k',lw=0.4,ax=ax,clip=zoom)

ax.axis(zoom)

frame_dir='attractor'
os.path.exists(frame_dir) or os.makedirs(frame_dir)

first=True
for frame,ti in enumerate(range(0,ntimes,2)):
    bspp_dupe.plot(ti,ax=ax,zoom=zoom,mask=~captured,color='cyan',alpha=0.4,ms=3,update=not first)
    bspp.plot(ti,ax=ax,zoom=zoom,mask=captured,color='m',alpha=0.4,ms=4,update=not first)
    #plt.draw()
    #plt.pause(0.001)
    first=False
    fig.savefig('%s/frame%04d.png'%(frame_dir,frame))
