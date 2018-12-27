import stompy.model.delft.dflow_model as dfm
import stompy.grid.unstructured_grid as ugrid
from stompy.plot import plot_utils,plot_wkb
from stompy.spatial import field

from shapely import geometry
import os

#dfm_bin_dir=os.path.join(os.environ['HOME'],"src/dfm/r53925-opt/bin")
dfm_bin_dir=os.path.join(os.environ['HOME'],"src/dfm/1.5.0/bin")
os.environ['LD_LIBRARY_PATH']=dfm_bin_dir.replace('bin','lib')


dfm.DFlowModel.dfm_bin_dir=dfm_bin_dir
dfm.DFlowModel.mpi_bin_dir=dfm_bin_dir

##
W=100
L=2000
nx=21
ny=2

def base_dem(res=2.0):
    F=np.zeros( (int(W/res),int(L/res)), np.float64)
    f=field.SimpleGrid(extents=[0,L,0,W],F=F)
    x,y=f.xy()
    f.F[:,:] = (-y/10)[:,None]
    return f

dem=base_dem()

##

def base_grid():
    g=ugrid.UnstructuredGrid(max_sides=4)

    g.add_rectilinear([0,0],[L,W],nx,ny)

    # z=0 along lower edge, and 10 on the upper edge.
    # node_depths=-g.nodes['x'][:,1]/10.
    node_depths=dem(g.nodes['x'])

    g.add_node_field('depth',node_depths)
    return g

##
def base_model(run_dir):
    model=dfm.DFlowModel()
    model.load_mdu('template.mdu')
    model.set_run_dir(run_dir,mode='pristine')
    model.run_start=np.datetime64("2000-01-01 00:00")
    model.run_stop =np.datetime64("2000-01-03 12:00")
    model.mdu['geometry','Conveyance2D'] = 2
    model.mdu['geometry','Nonlin2D']=1
    model.mdu['geometry','WaterLevIni'] = -9
    g=base_grid()
    model.set_grid(g)

    model.add_bcs( dfm.FlowBC(name='inflow',
                              geom=np.array( [ [0,0], [0,W]]) ,
                              Q=10.0))

    model.add_monitor_sections([dict(name='middle',geom=geometry.LineString( [ [L/2,0], [L/2,W]]))])
    model.add_monitor_points([dict(name='midpoint',geom=geometry.Point( [L/2,W/2] ))])
    return model

##

# run conveyance2d=-1, conveyance2d=2, and conveyance2d=2/Nonlin2d=1
mod_conv0=base_model('run-no_conv')
mod_conv0.mdu['geometry','Conveyance2D']=-1
mod_conv0.mdu['geometry','Nonlin2D']=0

mod_conv2=base_model('run-conv2')
mod_conv2.mdu['geometry','Conveyance2D']=2
mod_conv2.mdu['geometry','Nonlin2D']=0

mod_conv2nonlin=base_model('run-conv2-nonlin')
mod_conv2nonlin.mdu['geometry','Conveyance2D']=2
mod_conv2nonlin.mdu['geometry','Nonlin2D']=1

for model in [mod_conv0, mod_conv2, mod_conv2nonlin]:
    model.write()
    model.partition()
    model.run_model()

##

g=model.grid

# in the real deal, this gets more complicated. for now, hardcode
region_cells=g.cells_center()[:,0]>L/2

# And extract the "correct" answer from the DEM
section_ls=np.array( model.mon_sections[0]['geom'].coords )
region_poly=g.boundary_polygon_by_union(np.nonzero(region_cells)[0])

dem_in_poly = dem.polygon_mask(region_poly)
pixA=dem.dx*dem.dy
pix_z=dem.F[dem_in_poly]

def eta_to_V(eta):
    return (pixA*(eta-pix_z).clip(0,np.inf)).sum()

etas=np.linspace(-10,5,200)
dem_volumes=np.array([eta_to_V(eta) for eta in etas])


## 
from stompy.spatial import linestring_utils
section_segs=linestring_utils.resample_linearring(section_ls,0.5*dem.dx,closed_ring=0)
section_z=dem(section_segs)
section_s=utils.dist_along(section_segs)

def eta_to_xA(eta):
    return np.trapz( (eta-section_z).clip(0,np.inf), section_s)

dem_xas=np.array([eta_to_xA(eta) for eta in etas])

##

plt.figure(2).clf()
fig,axs=plt.subplots(3,2,sharex='col',sharey='row',num=2)

for model in [
        mod_conv0,
        mod_conv2,
        mod_conv2nonlin
]:
    his=xr.open_dataset(model.his_output())
    ds=xr.open_dataset(model.map_outputs()[0])
    
    mod_label=model.run_dir.replace('run-','')

    # Look at history to see if it's possible to extract hypsometry
    xareas=his.cross_section_area.isel(cross_section=0).values
    eta=   his.waterlevel.isel(stations=0).values

    axs[0,0].plot(his.time,eta,label=mod_label)
    axs[1,0].plot(his.time,xareas,label=mod_label)

    # volumes direct from the map
    # How comparable these are depends on the settings:
    #  conv2d=-1, bedlevtyp=3 => good, not exact, v_and_q constant offset below volumes.
    #  conv2d=2, nonlin2d=0 => as long as the freesurface is fairly flat, this is very close.
    #  conv2d=2, nonlin2d=1 => persistent and large difference, with an error in volumes
    #     due to the pseudo-subgrid treatment of the bed.
    volumes=((ds.mesh2d_flowelem_ba.isel(nmesh2d_face=region_cells) * ds.mesh2d_waterdepth.isel(nmesh2d_face=region_cells))).sum(dim='nmesh2d_face')
    # in the future will have to figure out sign of each section
    vol_and_Q=volumes[0].values + his.cross_section_cumulative_discharge.isel(cross_section=0)

    # axs[2,0].plot(volumes.time,volumes,label='%s - vols'%mod_label)
    axs[2,0].plot(vol_and_Q.time,vol_and_Q,label='%s'%mod_label)

    # And as a function of eta
    order=np.argsort(eta)
    axs[1,1].plot(eta[order], xareas[order],label=mod_label)
    axs[2,1].plot(eta[order], vol_and_Q[order], label="%s"%mod_label)
    # including map output - not reliable for nonlin2d=1
    #axs[2,1].plot( ds.mesh2d_s1.isel(nmesh2d_face=region_cells).mean(dim='nmesh2d_face'),
    #               volumes,label='%s vols'%mod_label )
    
    his.close()
    ds.close()

axs[2,1].plot(etas,dem_volumes,
              label='DEM',color='k',lw=0.6)
axs[1,1].plot(etas,dem_xas,label='DEM',color='k',lw=0.6)

axs[0,0].set_ylabel('eta')
axs[1,0].set_ylabel('xarea')
axs[2,0].set_ylabel('vol')
axs[2,1].set_xlabel('eta')

axs[0,0].legend()
axs[1,0].legend()
axs[2,0].legend()
axs[1,1].legend()
axs[2,1].legend()

##

# Ultimately I want plots that
# (1) show the volume ~ eta relationship for regions in the model
# (2) show the xs-area ~ eta relationship for cross-sections in the model
# and show the "correct" functions for both of those based on DEM.

