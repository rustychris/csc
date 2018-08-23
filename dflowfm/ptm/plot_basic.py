from stompy.grid import unstructured_grid
from stompy.model.fish_ptm import ptm_tools

##

g=unstructured_grid.UnstructuredGrid.from_ugrid('dfm_ptm_hydro.nc')

##
ges=ptm_tools.PtmBin('GES_bin.out')
ntimes=ges.count_timesteps()
##
# zoom=(629582, 629832, 4233134, 4233560)
# zoom=(629375.7037314154, 630090.3968338037, 4232811.347364194, 4234029.184410663)
# zoom=(628897.95707908, 630568.1434861392, 4232811.347364194, 4234029.184410663)
# zoom= (620820.082274368, 633427.7560025024, 4225915.179544365, 4235108.220555253)
zoom=(629674.547743109, 629827.6731611794, 4233121.869190959, 4233235.97232507)

plt.figure(1).clf()
fig,ax=plt.subplots(num=1)

ges.plot(ntimes-1,ax=ax,zoom=zoom,update=False)

g.plot_edges(color='k',lw=0.4,ax=ax,clip=zoom)

ax.axis(zoom)

##

for ti in range(0,2000,2):
    ges.plot(ti,ax=ax,zoom=zoom)
    plt.draw()
    plt.pause(0.001)
