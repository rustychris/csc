# Try to convert the relevant bits of DSS data into
# a netcdf with enough spatial information to put things onto
# a grid
# There is not a good way around depending on HEC-DSSVUE for this.
# First have to get DSSVUE to run dss_export_all.py, which
# pulls each timeseries out to an excel file
from stompy import utils
import xarray as xr
import pandas as pd
import numpy as np
import os, glob

import six
from stompy.model import dsm2
from stompy.grid import unstructured_grid
import matplotlib.pyplot as plt
from matplotlib import collections

##

dss=[]
for f in glob.glob('csv/*.csv'):
    print(f)
    # First read the header info
    meta={}
    with open(f,'rt') as fp:
        # key=value pairs in first line after # comment
        for kv in fp.readline().strip('#').strip().split():
            k,v=kv.split('=')
            meta[k]=v

    p=meta.pop('path')
    _,A,B,C,D,E,F,_ = p.split('/')
    try:
        B=int(B)
    except ValueError:
        print("%s not a valid node id. Skip"%B)
        continue
            
    df=pd.read_csv(f,skiprows=1,parse_dates=['time'])
    ds=xr.Dataset.from_dataframe(df.set_index('time'))

    ds['A']=(),A
    ds['B']=(),B
    ds['C']=(),C
    ds['E']=(),E
    ds['F']=(),F
    ds['value'].attrs.update(meta)

    dss.append(ds)


## 

# Bring in the geometry
# Do this before creating netcdf, b/c this will tell us the number
# of nodes required.
grid_inp='dsm2_gis/channel_std_delta_grid_NAVD_20150129.inp'
node_shp='dsm2_gis/dsm2_8_2_0_calibrated_model_grid/dsm2_8_2_0_calibrated_nodes.shp'
channel_shp='dsm2_gis/dsm2_8_2_0_calibrated_model_grid/dsm2_8_2_0_calibrated_channels_centerlines.shp'
dsm_grid=dsm2.DSM2Grid(grid_inp,node_shp,channel_shp=channel_shp)

## 

# Force same time axis, have the 3 types of flow (C) as separate
# variables

A_values=np.unique([ds['A'].values for ds in dss]) # just 1
valid_nodes=np.unique( [ds['B'].values for ds in dss])
Nnodes=1+np.max([ds['B'].values for ds in dss])
assert dsm_grid.Nnodes() >= max_node
Nnodes=dsm_grid.Nnodes()

C_values=np.unique([ds['C'].values for ds in dss]) # 'DIV-FLOW', 'DRAIN-FLOW', 'SEEP-FLOW'

# First, create a concatenated dataset for each of
# the unique C_values

ds=xr.Dataset()

# Just go for dense
ds['node']=('node',),np.arange(Nnodes)

date_min=np.min([ds['time'].values.min() for ds in dss])
date_max=np.max([ds['time'].values.max() for ds in dss])

dt=np.timedelta64(1,'D')

ds['time']=('time',),np.arange(date_min,date_max+dt,dt)

def C_to_var(C):
    return C.lower().replace('-','_')

for C_val in C_values:
    varname=C_to_var(C_val)
    fill=np.nan*np.ones( (ds.dims['time'],ds.dims['node']))
    ds[varname]=('time','node'),fill

for d in dss:
    # Where to shove it?
    varname=C_to_var(d['C'].item())
    ds[varname].loc[dict(time=d.time,node=d.B.values)]=d.value

##

# read in grid topology from inp, geometry from shps
# then maybe ready to see how it lines up with the
# CSC grid....

dsm_grid.write_to_xarray(ds)

## 

compress=dict(zlib=True,complevel=4)
ds.to_netcdf('dcd-1922_2016.nc',encoding=dict(seep_flow=compress,
                                              div_flow=compress,
                                              drain_flow=compress))

##

# Load the CSC grid to see how well things line up, and whether nodes
# need to be jittered at this stage.
# Get a grid with bathy to test out the node matching
g_csc=unstructured_grid.UnstructuredGrid.read_ugrid("../../dflowfm/CacheSloughComplex_v111-edit21-bathy.nc")
g_csc.edge_to_cells(recalc=True)



##

# Probably do need to look at "reservoirs" since I think that covers Liberty Island,
# so DCD terms there are very relevant. But the reservoirs are just given string labels,
# e.g. "Liberty Is."
# Does DSM2 include any sort of ET for the reservoirs?
# reservoir_std_delta_grid_NAVD_20121214.inp:
#   provides area (226.904) and BOT_ELEV (-7.652) for Liberty Island
#   And that it connects to node 322.

# Looks like DCD might supply some data related to clifton court forebay
# check for /DICU-HIST+RSVR in the dss files.
# Just the entries which appear to be for CCF (called BBID for reasons I don't know)
# Unclear if/how ET would be included in DSM2 for Liberty Island.
# For now, have to press on.

# Node around the split of prospect sl and liberty cut falls in between, while
# the channel centerlines go up liberty cut
# there's nothing in the toe drain.

# last node in DWSC is not in the DFM domain.
# A node on Miner Slough is not in the domain, but it's channel centerlines are in the domain


# The key question here is whether it is sufficient to nudge DCD nodes to nearest
# grid nodes, or if I need something manual.


## 
valid_dcd_nodes=( (np.isfinite(ds.seep_flow).any(dim='time')).values |
                  (np.isfinite(ds.drain_flow).any(dim='time')).values |
                  (np.isfinite(ds.div_flow).any(dim='time')).values )

valid_dcd_nodes=np.nonzero(valid_dcd_nodes)[0]


##
pairs=[]
bad_pairs=[]

z_cell=g_csc.interp_node_to_cell(g_csc.nodes['depth'])


cc=g_csc.cells_center()
thresh=400 # 500 includes one bad node at the DCC.  400 is good for the current grid.

for n in valid_dcd_nodes:
    dsm_x=dsm_grid.nodes['x'][n]
    c_near=g_csc.select_cells_nearest(dsm_x)

    # And then find the deepest cell nearby, subject
    # to distance from DSM point
    c_nbrs=[c_near]
    for _ in range(5):
        c_nbrs=np.unique([ nbr
                           for c in c_nbrs
                           for nbr in g_csc.cell_to_cells(c)
                           if nbr>=0 and utils.dist(dsm_x,cc[nbr])<thresh])
    if len(c_nbrs)==0:
        bad_pairs.append( [dsm_x, cc[c_near]] )
    else:
        match=c_nbrs[ np.argmin(z_cell[c_nbrs]) ]
        pairs.append( [dsm_x, cc[match]] )

if 1:

    plt.figure(1).clf()
    dsm_grid.plot_nodes()
    ecoll=dsm_grid.plot_edges(centerlines=True)

    g_csc.plot_cells(alpha=0.5,zorder=-1,color='0.7')

    plt.axis('tight')
    plt.axis('equal')

    lcoll1=collections.LineCollection(pairs,color='k',lw=3)
    lcoll2=collections.LineCollection(bad_pairs,color='m',lw=1.5)

    ax=plt.gca()
    ax.add_collection(lcoll1)
    ax.add_collection(lcoll2)


##


