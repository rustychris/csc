"""
"""
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from stompy import utils
from stompy.grid import unstructured_grid
from stompy.model.fish_ptm import ptm_tools
from stompy.plot import plot_wkb
from stompy.spatial import proj_utils

import logging
utils.log.setLevel(logging.INFO)

##
#g=unstructured_grid.UnstructuredGrid.from_ugrid('../../dflowfm/runs/20180807_grid98_17/ptm_hydro.nc')
#grid_poly=g.boundary_polygon_by_edges()
#ll2utm=proj_utils.mapper('WGS84','EPSG:26910')

##
ptm_run_dir="../rate_estimator/run_10days"
ptm_groups=["INIT","SAC","SRV"]

ptm_data=[ ptm_tools.PtmBin(os.path.join(ptm_run_dir,grp+"_bin.out"))
           for grp in ptm_groups]

##

ntimes=ptm_data.count_timesteps()

##
# filter out exited particles
poly_geo=np.array([[ 629664., 4233211.],
                   [ 629823., 4233215.],
                   [ 629821., 4233119.],
                   [ 629661., 4233114.]])
poly_rio=np.array([[ 615061., 4224762.],
                   [ 616046., 4224084.],
                   [ 615810., 4223645.],
                   [ 614768., 4224230.]])

# fast lookups via matplotlib:
dead_polys=[poly_geo,poly_rio]
from matplotlib import path
dead_paths=[path.Path(poly) for poly in dead_polys]

# dead_cells |= .contains_points(ctrs)

##


data_dir="csvs"
os.path.exists(data_dir) or os.mkdir(data_dir)

for step in utils.progress(range(1000)):
    fn=os.path.join(data_dir,"combined_%04d.csv"%step)
    if os.path.exists(fn): continue

    dfs=[]
    for src_i,src in enumerate(ptm_data):
        t,parts=src.read_timestep(step)
        sel=np.ones( len(parts), np.bool8 )
        for dp in dead_paths:
            sel &= ~dp.contains_points(parts['x'][:,:2])

        df=pd.DataFrame()
        df['x']=parts['x'][sel,0]
        df['y']=parts['x'][sel,1]
        df['group']=src_i
        dfs.append(df)

    comb=pd.concat(dfs)
    comb.to_csv(fn,index=False,float_format="%.1f")

