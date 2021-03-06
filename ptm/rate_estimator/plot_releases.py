from stompy.grid import unstructured_grid
from stompy.model.fish_ptm import ptm_tools
from stompy.plot import plot_wkb

##

g=unstructured_grid.UnstructuredGrid.from_ugrid('../../dflowfm/runs/20180807_grid98_17/ptm_hydro.nc')

##
init=ptm_tools.PtmBin('run_10days/INIT_bin.out')
sac=ptm_tools.PtmBin('run_10days/SAC_bin.out')
srv=ptm_tools.PtmBin('run_10days/SRV_bin.out')
ntimes=init.count_timesteps()
##

poly=g.boundary_polygon()

##

# zoom=(605889.6569457075, 638002.2586920519, 4217801.158715993, 4241730.226468915)
zoom=(597913.7274775933, 648118.8262812896, 4217179.54644355, 4301202.344200377)
# zoom=(611280.377359663, 632614.9072567355, 4222938.787804629, 4248182.140275016)
plt.figure(1).clf()
fig,ax=plt.subplots(num=1)

init.plot(ntimes-1,ax=ax,zoom=zoom,update=False,ms=1)
sac.plot(ntimes-1,ax=ax,zoom=zoom,update=False,color='cyan',ms=1)
srv.plot(ntimes-1,ax=ax,zoom=zoom,update=False,color='g',ms=1)

# g.plot_edges(color='k',lw=0.4,ax=ax,clip=zoom)
plot_wkb.plot_polygon(poly,ax=ax,fc='0.8',ec='k',lw=0.5)

ax.axis(zoom)

##

for ti in range(0,ntimes,2):
    init.plot(ti,ax=ax,zoom=zoom)
    sac.plot(ti,ax=ax,zoom=zoom)
    srv.plot(ti,ax=ax,zoom=zoom)
    plt.draw()
    plt.pause(0.001)

##


