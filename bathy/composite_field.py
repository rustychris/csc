"""
Combine other sources based on csc_bathy_sources, trim to domain,
and output 2m DEM.

"""
import numpy as np
import six
import subprocess

from stompy.spatial import field
import matplotlib.pyplot as plt

six.moves.reload_module(field)

import os
opj=os.path.join

tif_paths=[".",
           "/home/rusty/data/bathy_dwr/gtiff",
           ("/home/rusty/mirrors/ucd-X/Arc_Hydro/CSC_Project/MODELING/"
            "1_Hydro_Model_Files/Geometry/Bathymetry_Tiles/NDelta_hydro_2m_v4")]
dwr_dem_path="/home/rusty/data/bathy_dwr/gtiff"

##

# Load grid to find minimal area to output
from stompy.grid import unstructured_grid
from stompy.model.delft import dfm_grid

g=unstructured_grid.UnstructuredGrid.from_ugrid('../grid/CacheSloughComplex_v111-edit19fix.nc')

poly=g.boundary_polygon()

poly_buff=poly.buffer(100.0)

##

res=dx=2
xyxy=np.array(poly_buff.bounds)
xyxy=dx*np.round(xyxy/dx)
bleed=50
total_bounds=[xyxy[0],xyxy[2],xyxy[1],xyxy[3]]
xxyy_pad=[total_bounds[0]-bleed,total_bounds[1]+bleed,
          total_bounds[2]-bleed,total_bounds[3]+bleed]

##

def factory(attrs):
    geo_bounds=attrs['geom'].bounds

    if attrs['src_name']=='dwr_2m_dems':
        mrf=field.MultiRasterField([os.path.join(dwr_dem_path,"*_2m*tif")])
        return mrf
    elif attrs['src_name']=='dwr_10m_dems':
        #mrf=field.MultiRasterField([os.path.join(dwr_dem_path,"dem_bay_delta_10m_v3_20121109_*.tif")])
        # In this case there is only one that matters -- tile 4.
        gg=field.GdalGrid(os.path.join(dwr_dem_path,"dem_bay_delta_10m_v3_20121109_4.tif"),
                          geo_bounds=geo_bounds)
        gg.default_interpolation='linear'
        return gg
    elif attrs['src_name'].endswith('.tif'):
        for p in tif_paths:
            fn=opj(p,attrs['src_name'])
            if os.path.exists(fn):
                gg=field.GdalGrid(fn,geo_bounds=geo_bounds)
                gg.default_interpolation='linear'
                return gg
        raise Exception("Could not find tif %s"%attrs['src_name'])
    elif attrs['src_name'].startswith('py:'):
        expr=attrs['src_name'][3:]
        # something like 'ConstantField(-1.0)'
        # a little sneaky... make it look like it's running
        # after a "from stompy.spatial.field import *"
        # and also it gets fields of the shapefile
        field_hash=dict(field.__dict__)
        # convert the attrs into a dict suitable for passing to eval
        attrs_dict={}
        for name in attrs.dtype.names:
            attrs_dict[name]=attrs[name]
        return eval(expr,field_hash,attrs_dict)

    assert False

#src_shp='csc_bathy_sources.shp'
src_shp='csc_bathy_sources_post2014.shp'

mbf=field.CompositeField(shp_fn=src_shp,
                         factory=factory,
                         priority_field='priority',
                         data_mode='data_mode',
                         alpha_mode='alpha_mode')
##

def f(args):
    fn,xxyy,res = args

    bleed=150 # pad out the tile by this much to avoid edge effects

    if not os.path.exists(fn):
        #try:
        xxyy_pad=[ xxyy[0]-bleed,
                   xxyy[1]+bleed,
                   xxyy[2]-bleed,
                   xxyy[3]+bleed ]
        dem=mbf.to_grid(dx=res,dy=res,bounds=xxyy_pad,
                        mask_poly=poly_buff)
        # not great, since we're not padding the borders, but
        # there is very little filling now that the 2m dataset
        # is so pervasive.
        # dem.fill_by_convolution(iterations='adaptive',smoothing=2,kernel_size=7)

        if bleed!=0:
            dem=dem.crop(xxyy)
        dem.write_gdal(fn)
        #except Exception as exc:
        #    print("Failed with exception")
        #    print(repr(exc))

##
if 1: # __name__ == '__main__':
    dem_dir="tiles_2m_20190122"
    os.path.exists(dem_dir) or os.mkdir(dem_dir)

    res=2.0
    tile_x=tile_y=2000.0

    total_tile_bounds= [np.floor(total_bounds[0]/tile_x) * tile_x,
                        np.ceil(total_bounds[1]/tile_x) * tile_x,
                        np.floor(total_bounds[2]/tile_y) * tile_y,
                        np.ceil(total_bounds[3]/tile_y) * tile_y ]

    calls=[]
    for x0 in np.arange(total_tile_bounds[0],total_tile_bounds[1],tile_x):
        for y0 in np.arange(total_tile_bounds[2],total_tile_bounds[3],tile_y):
            xxyy=(x0,x0+tile_x,y0,y0+tile_y)
            fn=os.path.join(dem_dir,"%.0f_%.0f.tif"%(x0,y0))
            calls.append( [fn,xxyy,res] )

    if 0:
        p = Pool(4)
        p.map(f, calls )
    else:
        for i,call in enumerate(calls):
            print("Call %d/%d"%(i,len(calls)))
            f(call)


            # and then merge them with something like:
            # if the file exists, its extents will not be updated.
    output_fn='merged_2m-20190122.tif'
    os.path.exists(output_fn) and os.unlink(output_fn)
    subprocess.call("gdal_merge.py -init nan -a_nodata nan -o %s %s/*.tif"%(output_fn,dem_dir),
                    shell=True)


