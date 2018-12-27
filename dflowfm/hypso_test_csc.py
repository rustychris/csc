"""
Adapt the test-case code in hypso_test.py to the CSC domain.
"""
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

from shapely import geometry
import os

import stompy.model.delft.dflow_model as dfm
import stompy.grid.unstructured_grid as ugrid
from stompy import utils
from stompy.plot import plot_utils,plot_wkb
from stompy.spatial import linestring_utils
from stompy.spatial import field

##
class DfmHypso(object):
    def __init__(self,model):
        self.model=model
        if model.num_procs==1:
            self.ds=xr.open_dataset(self.model.map_outputs()[0])
            self.g=ugrid.UnstructuredGrid.read_dfm(self.ds)
        else:
            self.ds=None
            self.g=ugrid.UnstructuredGrid.read_dfm(model.mdu.filepath(['geometry','NetFile']))
            
        self.his=xr.open_dataset(self.model.his_output())
        self.g_poly=self.g.boundary_polygon()
        self.get_cross_sections()
        self.get_xs_edge_mapping()
        self.label_cells()
        
    def get_cross_sections(self):
        """
        extract coordinates of cross sections from the history output
        populates self.xs_xys as [ Nx2 array, Nx2 array, ... ]
        """
        self.xs_xys=[]

        for xs_idx in self.his.cross_section:
            pnts_x=self.his.cross_section_x_coordinate.isel(cross_section=xs_idx)
            pnts_y=self.his.cross_section_y_coordinate.isel(cross_section=xs_idx)
            valid=pnts_x<1e30
            xs_xy=np.c_[pnts_x[valid], pnts_y[valid]]
            self.xs_xys.append(xs_xy)
    def plot_cross_sections(self,ax=None):
        ax=ax or plt.gca()
        for xs_xy in self.xs_xys:
            ax.plot(xs_xy[:,0],xs_xy[:,1],color='r',lw=0.8)

    def get_xs_edge_mapping(self):
        """
        Annotate grid edges with membership in a cross-section, saving
        the mapping to self.edge_to_xs
        """
        self.edge_to_xs=edge_to_xs=np.zeros( self.g.Nedges(), dtype=[('xs',np.int32),
                                                                     ('sign',np.int32)])
        edge_to_xs['xs']=-1

        for xs_idx,xs_xy in enumerate(self.xs_xys):
            node_string=[]
            for a,b in zip( xs_xy[:-1:,:], xs_xy[1:,:]):
                na=self.g.select_nodes_nearest(a)
                nb=self.g.select_nodes_nearest(b)
                b_to_a=self.g.shortest_path(na,nb)
                a_to_b=b_to_a[::-1]
                node_string.extend(a_to_b)
            # each path is inclusive of endpoints, so expect those to get
            # repeated
            node_string=utils.remove_repeated(node_string)
            # but two segments could trace to the same node, and that's bad
            assert len(node_string)==len(np.unique(node_string))

            for na,nb in zip(node_string[:-1],node_string[1:]):
                j=self.g.nodes_to_edge([na,nb])
                edge_to_xs['xs'][j]=xs_idx
                if self.g.edges['nodes'][j][0]==na:
                    s=1
                elif self.g.edges['nodes'][j][0]==nb:
                    s=-1
                else: raise Exception("what?")
                edge_to_xs['sign'][j]=s
    def plot_xs_edges(self,ax=None):
        ax=ax or plt.gca()
        xs=self.edge_to_xs['xs']
        return self.g.plot_edges(mask=(xs>=0),
                                 values=xs,lw=4,cmap='jet')
        
    def label_cells(self):
        # Label cells separated by sections
        cell_labels=self.g.cells_connected_components(edge_mask=(self.edge_to_xs['xs']<0),
                                                      randomize=False)
        self.cell_labels=cell_labels.data.astype(np.int32)

    def plot_cell_labels(self):
        self.g.plot_cells(values=self.cell_labels,cmap='jet')

    def nearest_station(self,origin):
        g=self.g
        his=self.his.isel(time=0)
        station_nodes=[]
        for station in his.stations.values:
            xy=np.array( [his.station_x_coordinate.isel(stations=station),
                          his.station_y_coordinate.isel(stations=station) ] )
            n=g.select_nodes_nearest(xy)
            station_nodes.append(n)
        o_node=g.select_nodes_nearest( origin )

        results=g.shortest_path(o_node,station_nodes,max_return=1,return_type='cost')
        node,dist = results[0]
        station=station_nodes.index(node)

        station_name=self.his.station_name.values[station]
        try:
            station_name=station_name.decode()
        except AttributeError:
            pass

        return station,dist,station_name
        
class Region(object):
    def __init__(self,hyp,region_i,bathy_fn):
        self.hyp=hyp
        self.region_i=region_i

        self.region_cells=(self.hyp.cell_labels==self.region_i)

        self.dem_analysis(bathy_fn)
        self.vol_and_Q_analysis()
        
    def plot_region_cells(self,**kw):
        defaults=dict(color='orange',alpha=0.5)
        defaults.update(kw)
        return self.hyp.g.plot_cells(mask=self.region_cells,**defaults)

    def get_region_dem(self,bathy_fn):
        # And extract the "correct" answer from the DEM
        region_poly=self.hyp.g.boundary_polygon_by_union(np.nonzero(self.region_cells)[0])
        self.region_poly=region_poly

        xyxy=region_poly.bounds

        region_dem=field.GdalGrid(bathy_fn,geo_bounds=[xyxy[0],xyxy[2],xyxy[1],xyxy[3]])
        self.dem=region_dem

        dem_in_poly = region_dem.polygon_mask(region_poly)

        dem_bad=np.isnan(region_dem.F)
        if np.any(dem_bad & dem_in_poly):
            print("%d of %d dem pixels are missing"%( (dem_bad & dem_in_poly).sum(),
                                                      dem_in_poly.sum() ))
            dem_in_poly=dem_in_poly & (~dem_bad)

        pix_A=region_dem.dx*region_dem.dy
        pix_z=region_dem.F[dem_in_poly]
        self.pix_A=pix_A
        self.pix_z=pix_z
        return pix_A,pix_z
    
    def dem_analysis(self,bathy_fn):
        pix_A,pix_z = self.get_region_dem(bathy_fn)
        
        def eta_to_V(eta):
            return (pix_A*(eta-pix_z).clip(0,np.inf)).sum()
        self.dem_eta_to_V = eta_to_V
        
        # faster approach:
        order=np.argsort(pix_z)
        
        dem_etas=pix_z[order]
        dem_volumes=np.zeros_like(dem_etas)
        dem_volumes[0]=0.0
        for i in range(1,len(pix_z)):
            A=i*pix_A # wet area is all previous pixels
            dV=A*(dem_etas[i]-dem_etas[i-1])
            dem_volumes[i]=dem_volumes[i-1]+dV

        # for plotting, interpolate to coarser
        self.dem_etas=np.arange(dem_etas[0],dem_etas[-1],0.01)
        self.dem_volumes=np.interp(self.dem_etas,dem_etas,dem_volumes)
        
    def map_volume0(self):
        """
        Calculated volume of region at start of run from the map file.
        This can be pretty wrong, though.
        """
        ds=self.hyp.ds
        volume0=(ds.mesh2d_flowelem_ba * ds.mesh2d_waterdepth.isel(time=0)).isel(nmesh2d_face=self.region_cells).sum()
        return volume0

    def bounding_cross_sections(self):
        """
        returns dictionary of cross-section id => sign for flow into this region
        """
        g=self.hyp.g
        
        g.edge_to_cells()

        # edges that straddle the boundary of this region
        bound_edges=(self.region_cells[ g.edges['cells'][:,0] ]) != ( self.region_cells[ g.edges['cells'][:,1] ] )
        # and are not closed boundaries
        bound_edges= bound_edges & (g.edges['cells'].min(axis=1)>=0)

        bound_xs={} # map xs id to sign
        for j in np.nonzero(bound_edges)[0]:
            xs=self.hyp.edge_to_xs['xs'][j]
            if xs<0:
                print("Got some bad edges")
                continue
            # Note that edge_to_xs gives sign relative to the *edge*
            j_sgn=self.hyp.edge_to_xs['sign'][j]
            c1,c2=g.edges['cells'][j]

            # sign convention by trial and error with DWS
            if self.region_cells[c1]:
                sgn=-j_sgn
            elif self.region_cells[c2]:
                sgn=j_sgn
            else: raise Exception("bad")

            if xs in bound_xs:
                assert bound_xs[xs]==sgn
            else:
                bound_xs[xs]=sgn
        self.bound_xs=bound_xs
        return bound_xs

    def get_region_station(self):
        # and need an eta value for the region
        region_station=None
        his=self.hyp.his
        his=his.isel(time=0) # in case station coords have time dim
        for station in his.stations.values:
            xy=np.array( [his.station_x_coordinate.isel(stations=station),
                          his.station_y_coordinate.isel(stations=station) ] )
            c=self.hyp.g.select_cells_nearest(xy)
            if self.region_cells[c]:
                region_station=station
                break
        else:
            c=np.nonzero(self.region_cells)[0][0]
            region_station,dist,station_name=self.hyp.nearest_station(self.hyp.g.cells_center()[c])
        self.region_station=region_station
    
    def vol_and_Q_analysis(self):
        # Which sections are adjacent to this region, and what is the sign for cumulative discharge
        # *into* the volume?

        self.get_region_station()
        
        eta=self.hyp.his.waterlevel.isel(stations=self.region_station).values
        
        bound_xs=self.bounding_cross_sections()

        vol_and_Q=0.0
        for xs in list(bound_xs.keys()):
            sgn=bound_xs[xs]
            vol_and_Q = vol_and_Q + sgn*self.hyp.his.cross_section_cumulative_discharge.isel(cross_section=xs)

        if 0:
            # volume0 isn't very reliable.
            vol_offset=self.map_volume0()
        else:
            # it's really the prism we care about, not the subtidal volume,
            # so pull volume0 from the DEM
            # further, there are sometimes transient effects at the start, so line this
            # up at the end
            vol_offset=self.dem_eta_to_V(eta[-1]) - vol_and_Q.values[-1]
        
        vol_and_Q=vol_and_Q+vol_offset
        self.vol_and_Q=vol_and_Q
        self.eta=eta

    def calculate_cell_volumes(self):
        """
        Calculate the eta/volume relationship based on cell depths
        for the full range of elevations present in the dem.
        return tuple (etas,volumes)
        """
        g=self.hyp.g
        if self.hyp.ds is not None:
            reg_cell_z=self.hyp.ds.mesh2d_flowelem_bl.values[self.region_cells]
        else:
            # have to go to nodes
            print("Getting cell depth from average of nodes - assumes bedlevtyp=6")
            node_depth=g.nodes['depth']
            reg_cell_z=np.array( [ node_depth[g.cell_to_nodes(c)].mean()
                                   for c in np.nonzero(self.region_cells)[0] ] )

        cell_vols=[]
        areas=g.cells_area()[self.region_cells]
        for eta in self.dem_etas:
            depths=(eta-reg_cell_z).clip(0,np.inf)
            cell_vols.append( (depths*areas).sum() )
        cell_vols=np.array(cell_vols)
        return self.dem_etas,cell_vols
        
    def summary_figure(self,num=3):
        fig=plt.figure(num)
        fig.clf()

        ax_map=fig.add_subplot(1,2,1)
        ax_V=fig.add_subplot(1,2,2)

        xyxy=self.region_poly.bounds
        clip=[xyxy[0],xyxy[2],xyxy[1],xyxy[3]]
        hyp=self.hyp
        
        # hyp.g.plot_edges(clip=clip,ax=ax_map,color='k',lw=0.4,alpha=0.4)
        plot_wkb.plot_wkb( hyp.g_poly, fc='0.8', ec='none',ax=ax_map,zorder=-4)
        hyp.g.plot_cells(mask=self.region_cells,ax=ax_map,color='orange')

        # Show the location of the reference station
        ax_map.plot([self.hyp.his.station_x_coordinate.isel(stations=self.region_station)],
                    [self.hyp.his.station_y_coordinate.isel(stations=self.region_station)],
                    'go')
        # bound_xs=self.bounding_cross_sections()
        xys=[]
        uvs=[]
        for k in self.bound_xs.keys():
            sgn=self.bound_xs[k]
            sec_x=self.hyp.his.cross_section_x_coordinate.isel(cross_section=k).values
            sec_y=self.hyp.his.cross_section_y_coordinate.isel(cross_section=k).values
            sec_xy=np.c_[sec_x,sec_y]
            sec_xy=sec_xy[sec_xy[:,0]<1e10]
            ax_map.plot(sec_xy[:,0],sec_xy[:,1],'r-')
            uvs.append(sgn*(sec_xy[-1] - sec_xy[0]))
            xys.append(sec_xy.mean(axis=0))
        xys=np.array(xys)
        uvs=np.array(uvs)
        uvs=utils.rot(-np.pi/2,utils.to_unit(uvs))
        ax_map.quiver(xys[:,0],xys[:,1],uvs[:,0],uvs[:,1])

        ax_map.axis('equal')

        order=np.argsort(self.eta)
        ax_V.plot(self.eta[order], self.vol_and_Q[order], label="model")
        ax_V.plot(self.dem_etas,self.dem_volumes,
                  label='DEM',color='k',lw=0.6)

        # and what we expect from the model, assuming nonlin2d=0
        cell_etas,cell_vols = self.calculate_cell_volumes()
        ax_V.plot(cell_etas,cell_vols,
                  label='Cell depth',color='r',lw=0.6)

        ax_V.set_ylabel('Vol')
        ax_V.set_xlabel('eta')

        ax_V.legend()
        return fig

# similar for cross-sections.
class Section(object):
    def __init__(self,hyp,section_i,bathy_fn):
        self.hyp=hyp
        self.section_i=section_i
        self.section_name=self.hyp.his.cross_section_name.values[self.section_i]
        try:
            self.section_name=self.section_name.decode()
        except AttributeErro:
            pass

        self.section_ls=self.hyp.xs_xys[self.section_i]
        self.get_section_station()
        self.dem_analysis(bathy_fn)
        self.A_analysis()

    def dem_analysis(self,bathy_fn):
        # 1.0 is from half the current DEM resolution
        pad=10.0
        bounds=[self.section_ls[:,0].min()-pad,
                self.section_ls[:,0].max()+pad,
                self.section_ls[:,1].min()-pad,
                self.section_ls[:,1].max()+pad]
        self.dem=field.GdalGrid(bathy_fn,geo_bounds=bounds)
        res=0.5*self.dem.dx
        self.section_segs=linestring_utils.resample_linearring(self.section_ls,res,closed_ring=0)
        self.section_z=self.dem(self.section_segs)
        self.section_s=utils.dist_along(self.section_segs)

    
    def get_section_station(self):
        """choose a point station that is representative of this section
        """
        self.section_station,dist,station_name=self.hyp.nearest_station( self.section_ls.mean(axis=0) )
        
        print("Match section %s => station %s, distance=%.1f"%(self.section_name,
                                                               station_name,
                                                               dist))
        
    def dem_eta_to_xA(self,eta):
        return np.trapz( (eta-self.section_z).clip(0,np.inf), self.section_s)

    def A_analysis(self):
        """ Extract model cross-section area and eta.
        """
        self.areas=self.hyp.his.cross_section_area.isel(cross_section=self.section_i)
        self.eta=self.hyp.his.waterlevel.isel(stations=self.section_station).values

    def summary_figure(self,num=3):
        fig=plt.figure(num)
        fig.clf()

        ax_map=fig.add_subplot(1,2,1)
        ax_A=fig.add_subplot(1,2,2)

        #
        pnts=self.section_ls
        L=5*self.section_s.max()
        clip=[ pnts[:,0].min() - L,
               pnts[:,0].max() + L,
               pnts[:,1].min() - L,
               pnts[:,1].max() ]

        plot_wkb.plot_wkb( hyp.g_poly, fc='0.8', ec='none',ax=ax_map,zorder=-4)
        # hyp.g.plot_cells(mask=self.region_cells,ax=ax_map,color='orange')
        ax_map.plot( self.section_ls[:,0],self.section_ls[:,1],'b-')

        # Show the location of the reference station
        ax_map.plot([self.hyp.his.station_x_coordinate.isel(stations=self.section_station)],
                    [self.hyp.his.station_y_coordinate.isel(stations=self.section_station)],
                    'go')

        ax_map.text( self.section_ls[0,0], self.section_ls[0,1], self.section_name)
        ax_map.axis('equal')
        ax_map.axis(clip)

        order=np.argsort(self.eta)
        ax_A.plot(self.eta[order], self.areas[order], label="model")


        dem_xas=np.array( [self.dem_eta_to_xA(eta) for eta in self.eta] )
        ax_A.plot(self.eta[order],dem_xas[order],label='DEM',color='k',lw=0.6)

        ax_A.set_ylabel('area')
        ax_A.set_xlabel('eta')

        ax_A.legend()

        return fig

## 

bathy_fn="../bathy/merged_2m-20181113.tif"

model=dfm.DFlowModel.load('runs/grid100_20')

hyp=DfmHypso(model)

##

fig_dir=os.path.join(model.run_dir,"figs-conveyance")

os.path.exists(fig_dir) or os.makedirs(fig_dir)

plt.ioff()
region_ids=np.unique(hyp.cell_labels)
force=False
for reg_i in region_ids:
    img_fn= os.path.join(fig_dir,"region-%03d.png"%reg_i)
    if not force and os.path.exists(img_fn): continue
    print(reg_i)
    reg=Region(hyp,reg_i,bathy_fn)
    reg.summary_figure(10+reg_i)
    plt.savefig(img_fn,dpi=150)
    plt.close('all')
plt.ion()

## 
# summary figure for a section
plt.ioff()
for section_i in hyp.his.cross_section.values:
    section=Section(hyp,section_i,bathy_fn)
    fig=section.summary_figure(30+section_i)
    fig.savefig(os.path.join(fig_dir,"section-%03d.png"%section_i),dpi=150)

##

reg_i=5
region=Region(hyp,reg_i,bathy_fn)
## 
region.summary_figure(10+reg_i)
fig=plt.gcf()

##

# can I recreate the same curve directly?
# region.eta

# confirmed that these are the same:
areas= hyp.g.cells_area()[region.region_cells]
ba=hyp.ds.mesh2d_flowelem_ba.isel(nmesh2d_face=region.region_cells).values

bl=hyp.ds.mesh2d_flowelem_bl.isel(nmesh2d_face=region.region_cells).values

##

def eta_to_cellV(eta):
    depths=(eta-bl).clip(0,np.inf)
    return (depths*areas).sum()

cell_volumes=np.array( [eta_to_cellV(e) for e in region.dem_etas] )

ax=fig.axes[1]
ax.plot(region.dem_etas,cell_volumes,color='r',label='Cells')

##

# calculate hypsometry-preserving cell elevations
# for the grid, one region at a time.
g=hyp.g
hyp_cell_z=np.zeros(g.Ncells(),np.float64)
orig_node_z=hyp.ds.mesh2d_node_z.values # could also come from point-sampling DEM
orig_cell_z=g.interp_node_to_cell(orig_node_z)
Ac=g.cells_area()

region_ids=np.unique(hyp.cell_labels)
for reg_i in region_ids:
    print("Processing cells in region %d"%reg_i)
    region=Region(hyp,reg_i,bathy_fn)
    cell_idxs=np.nonzero(region.region_cells)[0]

    region_orig_cell_z=orig_cell_z[cell_idxs]

    # first, need to get the full range of etas from the DEM.
    pix_A=region.pix_A
    pix_z=region.pix_z
    pix_z_sort=np.sort(pix_z) # in ascending order

    # sort the cells by original elevation
    region_cell_order=np.argsort(region_orig_cell_z)
    region_Ac=Ac[cell_idxs]

    area_sum=0.0
    for ci in utils.progress(region_cell_order):
        pix_idx=int(round( (area_sum+0.5*region_Ac[ci])/pix_A ))
        if pix_idx<0:
            print("Whoa - negative?")
            pix_idx=0
        elif pix_idx>=len(pix_z_sort):
            print("Whoa - off the top end")
            pix_idx=len(pix_z_sort)-1

        c=cell_idxs[ci]
        hyp_cell_z[c]=pix_z_sort[pix_idx]
        area_sum+=region_Ac[ci]
    
##

# Plot the results:

orig_cell_z=hyp.ds.mesh2d_flowelem_bl.values

fig=plt.figure(10)
fig.clf()
fig,axs=plt.subplots(1,3,sharex=True,sharey=True,num=10)

g.plot_cells(values=orig_cell_z,ax=axs[0],cmap='jet',clim=[-10,5])
g.plot_cells(values=hyp_cell_z ,ax=axs[1],cmap='jet',clim=[-10,5])
g.plot_cells(values=(hyp_cell_z-orig_cell_z),ax=axs[2],cmap='seismic',clim=[-2,2])

axs[0].axis('equal')
fig.tight_layout()

##

# Still need to translate that to node elevations.

# want a matrix that takes node elevations, and averages to cell elevations.
# then it's a short trip to get to volume

#node_z_to_cell_z=np.zeros( (g.Ncells(),g.Nnodes()), np.float64)
from scipy import sparse
g=hyp.g
node_z_to_cell_z=sparse.dok_matrix( (g.Ncells(),g.Nnodes()), np.float64 )
for c in utils.progress(range(g.Ncells())):
    nodes=g.cell_to_nodes(c)
    node_z_to_cell_z[c,nodes]=1./len(nodes)

# test that -- yep, works.
# mat_cell_z=np.dot(node_z_to_cell_z,region_node_z)

##

res=sparse.linalg.lsqr(node_z_to_cell_z.tocsr(),
                       hyp_cell_z)
# with no damping, node elevations are -38 to 36.
# with damp=0.3, the range of node values is -26 -- 14,
# seems pretty decent. but this has a distinct bias
# since the expected node elevations are not centered on zero.
# 
hyp_node_z, istop, itn, r1norm  = res[:4]

res_z=node_z_to_cell_z.dot(hyp_node_z)-hyp_cell_z

##

g.add_node_field('depth',hyp_node_z)

g.write_ugrid('CacheSloughComplex_v111-edit19-hyp_bathy.nc')

##

plt.figure(11).clf()
g.plot_cells(values=res_z,cmap='seismic',clim=[-2,2])
plt.axis('equal')


##
def new_eta_to_cellV(eta):
    depths=(eta-hyp_cell_z[region.region_cells]).clip(0,np.inf)
    return (depths*areas).sum()

new_cell_volumes=np.array( [new_eta_to_cellV(e) for e in region.dem_etas] )

fig=plt.figure(15)
ax=fig.axes[1]
ax.plot(region.dem_etas,new_cell_volumes,'m-')

##

cell_etas,cell_vols = calculate_cell_volumes(reg)

##

# Is it close enough to calculate cell means directly?
region_dem=region.dem
cell_mean_from_dem=np.zeros(g.Ncells(),np.float64)

for c in utils.progress(np.nonzero(region.region_cells)[0]):
    msk=region_dem.polygon_mask(g.cell_polygon(c))
    cell_mean_from_dem[c]=np.nanmean(region_dem.F[msk])

cell_mean_from_dem[np.isnan(cell_mean_from_dem)]=0.0    
##

def celldem_eta_to_cellV(eta):
    depths=(eta-cell_mean_from_dem[region.region_cells]).clip(0,np.inf)
    return (depths*areas).sum()

celldem_cell_volumes=np.array( [celldem_eta_to_cellV(e) for e in region.dem_etas] )

fig=plt.figure(15)
ax=fig.axes[1]
ax.plot(region.dem_etas,celldem_cell_volumes,'g-',label='Cell-DEM')


##

# recalculate DEM values
dem_Vs=[]
for e in region.eta:
    idx=np.searchsorted(pix_z_sort,e)
    dem_Vs.append( (e-pix_z_sort[:idx]).sum() * pix_A )
    
ax.plot(region.eta,dem_Vs,color='r')
