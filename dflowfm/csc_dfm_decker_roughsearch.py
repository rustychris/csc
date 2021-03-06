#!/usr/bin/env python
"""
Script to configure and generate mdu ready for running dflow fm
for Cache Slough Complex
 - *_decker version switches to a grid that extends down to decker island,
   necessitating different forcing config.

 - sensitivity test for roughness by region
"""

import os
import glob
import shutil
import numpy as np
import logging
import pandas as pd
import xarray as xr
import six
from collections import OrderedDict

from stompy import utils, filters
import stompy.model.delft.io as dio
from stompy.io.local import usgs_nwis
from stompy.model.delft import dfm_grid
from stompy.grid import unstructured_grid
from stompy.spatial import wkb2shp, field
from stompy.io.local import cdec

if __name__=='__main__':
    logging.basicConfig(level=logging.INFO)

log=logging.getLogger('csc_dfm')

import barker_data
import stompy.model.delft.dflow_model as dfm
import local_config
import stompy.model.hydro_model as hm
cache_dir='cache'

## --------------------------------------------------

# for DCC, when the gates are closed they sometimes do not report
# flow, so here we will those gaps with 0.
class FillGaps(hm.BCFilter):
    max_gap_interp_s=2*60*60
    large_gap_value=0.0

    def transform_output(self,da):
        # have self.bc, self.bc.model
        # self.bc.data_start, self.bc.data_stop
        if len(da)==0:
            log.warning("FillGaps called with no input data")
            da=xr.DataArray(self.large_gap_value)
            return da
        log.warning("FillGaps code incomplete")
        return da


def base_config(model):
    """
    model: Model instance, should already have run_dir set.

    """
    model.num_procs=4
    model.z_datum='NAVD88'
    model.projection='EPSG:26910'
    model.utc_offset=np.timedelta64(-8,'h') # PST

    model.run_start=np.datetime64('2017-12-01')
    # very short!
    model.run_stop=np.datetime64('2017-12-05')

    model.load_mdu('template.mdu')

    src_grid='../grid/CacheSloughComplex_v111-edit19fix.nc'
    dst_grid=os.path.basename(src_grid).replace('.nc','-bathy.nc')
    bathy_fn="../bathy/merged_2m-20181113.tif"
    if utils.is_stale(dst_grid,[src_grid,bathy_fn]):
        g=unstructured_grid.UnstructuredGrid.from_ugrid(src_grid)
        dem=field.GdalGrid(bathy_fn)
        # node_depths=dem(g.nodes['x'])
        import dem_cell_node_bathy
        node_depths=dem_cell_node_bathy.dem_to_cell_node_bathy(dem,g)
        g.add_node_field('depth',node_depths,on_exists='overwrite')
        g.write_ugrid(dst_grid,overwrite=True)
    else:
        g=unstructured_grid.UnstructuredGrid.from_ugrid(dst_grid)
    model.set_grid(g)

    # I think this doesn't matter with conveyance2d>0
    model.mdu['geometry','BedlevType']=6
    # ostensibly analytic 2D conveyance with 3.
    # 2 is faster, and appears less twitchy while showing very similar
    #  calibration results.
    # 0 blocks some flow, 1 was a little twitchy.
    model.mdu['geometry','Conveyance2D']=-1
    # enabling this seems to cause a lot of oscillation.
    model.mdu['geometry','Nonlin2D']=0

    model.mdu['output','MapInterval']=7200
    #model.mdu['output','WaqInterval']=1800
    model.mdu['physics','UnifFrictCoef']= 0.025

    model.set_cache_dir('cache')

    model.add_gazetteer('gis/model-features.shp')
    model.add_gazetteer('gis/point-features.shp')

    hm.SourceSinkBC.dredge_depth=-1
    hm.FlowBC.dredge_depth=-1

    # check_bspp.py has the code that converted original tim to csv.
    # awkward reloading of that, but at least it's independent of the
    # simulation period.
    model.add_bcs(barker_data.BarkerPumpsBC(name='Barker_Pumping_Plant'))

    # Decker only exists post-2015
    if model.run_start>np.datetime64("2015-11-16"):
        model.add_bcs(hm.NwisStageBC(name='decker',station=11455478,cache_dir='cache'))
    else:
        # maybe fall back to Rio Vista, or some adjustment thereof
        raise Exception("Decker tidal data starts 2015-11-16, too late for this simulation period")
    model.add_bcs(hm.NwisFlowBC(name='threemile',station=11337080,cache_dir='cache'))
    # GSS: from compare_flows and lit, must be flipped.
    model.add_bcs(hm.NwisFlowBC(name='Georgiana',station=11447903,cache_dir='cache',
                                filters=[hm.Transform(lambda x: -x)] ))

    model.add_bcs(hm.NwisFlowBC(name='dcc',station=11336600,cache_dir='cache',
                                     filters=[FillGaps(),
                                              hm.Transform(lambda x: -x)] ) )

    sac=hm.NwisFlowBC(name="SacramentoRiver",station=11447650,
                           pad=np.timedelta64(5,'D'),cache_dir='cache',
                           filters=[hm.LowpassGodin(),
                                    hm.Lag(np.timedelta64(-2*3600,'s'))])
    model.add_bcs(sac)

    ##
    pad=np.timedelta64(5,'D')
    lisbon_ds=cdec.cdec_dataset(station='LIS',
                                start_date=model.run_start-pad, end_date=model.run_stop+pad,
                                sensor=20, cache_dir=cache_dir)
    # to m3/s
    lisbon_ds['Q']=lisbon_ds['sensor0020'] * 0.028316847
    # to hourly average
    lisbon_ds.groupby( lisbon_ds.time.astype('M8[h]')).mean()
    lisbon_bc=hm.FlowBC(name='lis',Q=lisbon_ds.Q)

    # Roughness - handled by caller.
    # model.add_RoughnessBC(shapefile='forcing-data/manning_n.shp')

    # Culvert at CCS
    if 1:
        # rough accounting of restoration
        if model.run_start < np.datetime64("2014-11-01"):
            # name => id
            # polylinefile is filled in from shapefile and name
            model.add_Structure(name='ccs_breach',
                                type='gate',
                                door_height=15, # no overtopping?
                                lower_edge_level=0.89,
                                opening_width=0.0, # pretty sure this is ignored.
                                sill_level=0.85, # gives us a 0.05m opening?
                                horizontal_opening_direction = 'symmetric')
            # these are really just closed
            for levee_name in ['ccs_west','ccs_east']:
                model.add_Structure(name=levee_name,
                                    type='gate',
                                    door_height=15,
                                    lower_edge_level=3.5,
                                    opening_width=0.0,
                                    sill_level=3.5,
                                    horizontal_opening_direction = 'symmetric')

    mon_sections=model.match_gazetteer(monitor=1,geom_type='LineString')
    mon_points  =model.match_gazetteer(geom_type='Point')
    model.add_monitor_sections(mon_sections)
    model.add_monitor_points(mon_points)

    return model
##


run_coll_dir=os.path.join(local_config.run_dir_root,"roughsearch")
run_log_file=os.path.join(run_coll_dir,"runlog")

if os.path.exists(run_log_file):
    existing_runs=pd.read_csv(run_log_file)

os.path.exists(run_coll_dir) or os.makedirs(run_coll_dir)


## 11 regions
regions=wkb2shp.shp2geom('gis/roughness_regions.shp')

##

# parameters for one run associate each region name with a roughness

def baseline_settings():
    settings=OrderedDict()
    for region in regions:
        settings[region['name']]=region['n_nominal']
    return settings


def settings_to_roughness_xyz(model,settings):
    xy=model.grid.nodes['x']
    z=np.zeros(len(xy),np.float64)
    z[:]=np.nan

    for region in regions:
        sel=model.grid.select_nodes_intersecting(geom=region['geom'])
        z[sel]=settings[region['name']]
    missing=np.sum(np.isnan(z))
    if missing:
        print("%d nodes will get default of 0.02"%missing)
        z[np.isnan(z)]=0.02
    return np.c_[xy,z]


from stompy import memoize

def settings_to_txt(settings):
    names=list(settings.keys())
    names.sort()

    lines=[]
    for n in names:
        lines.append("settings['%s']=%r"%(n,settings[n]))
    return "\n".join(lines)

def check_and_run_settings(settings,retries=1):
    # don't pass as keyword argument, in case order is
    # scrambled
    key=memoize.memoize_key(settings=settings)
    run_dir=os.path.join(run_coll_dir,key)

    if dfm.DFlowModel.run_completed(run_dir):
        # print("%s exists -- will skip"%run_dir)
        return dfm.DFlowModel.load(run_dir)

    for retry in range(retries):
        print("running in %s"%run_dir)

        model=dfm.DFlowModel()
        model.set_run_dir(run_dir, mode='pristine')
        base_config(model)

        txt_fn=os.path.join(run_dir,"settings.txt")
        with open(txt_fn,'wt') as fp:
            fp.write( settings_to_txt(settings) )
        with open(txt_fn,'rt') as fp:
            print(fp.read())
        print("-"*20)

        xyz=settings_to_roughness_xyz(model,settings)
        # Turn that into a DataArray
        da=xr.DataArray( xyz[:,2],dims=['location'],name='n' )
        da=da.assign_coords(x=xr.DataArray(xyz[:,0],dims='location'),
                            y=xr.DataArray(xyz[:,1],dims='location'))
        da.attrs['long_name']='Manning n'

        model.add_RoughnessBC(data_array=da)

        model.write()
        print("---- About to partition----")
        model.partition()
        print("---- About to run -----")
        model.run_model()
        if model.is_completed():
            break
    else:
        print("looks like that didn't finish")
    return model

##---

# define a cost function for a run

def calc_model_metrics(model):
    import run_cal_plotter
    his_file=model.his_output()
    mr = run_cal_plotter.hcp.DflowfmModelResults([his_file],trim_spinup_days=1.0)

    defaults=run_cal_plotter.defaults

    recs=[]

    for station in run_cal_plotter.stations:
        args = station[1]
        args['ID'] = station[0]
        args['plot_dir']=os.path.dirname(his_file)
        for key in defaults:
            args.setdefault(key, defaults[key])  # append default parameters if missing

        p = run_cal_plotter.MyPlotter(mr, args)
        if not p.valid: continue

        if args['ID']=='TSL' and args['var']=='flow' and 'grid100_13' in his_file:
            print("Monkey patch sign on TSL flow")
            p.pred *= -1
        metrics=p.calculate_metrics()
        if metrics['ratio']<-1000: continue # probably no data

        rec={}
        rec['id']=args['ID']
        rec['var']=args['var']
        rec.update(metrics)
        recs.append(rec)
    return recs

def get_model_metrics(model):
    metrics_fn=os.path.join( os.path.dirname(model.his_output()),
                             "metrics.csv")
    if not os.path.exists(metrics_fn):
        recs=calc_model_metrics(model)
        df=pd.DataFrame(recs)
        df.to_csv(metrics_fn)
        if len(df)==0:
            import pdb
            pdb.set_trace()
    else:
        df=pd.read_csv(metrics_fn)
    return df

def get_model_score(model):
    metrics=get_model_metrics(model)
    return metrics_to_score(metrics)

def metrics_to_score(metrics):
    # Easiest is to use skill, but it gets less sensitive as
    # you get close, and doesn't really penalize a 20 minute lag
    # return np.mean((1-metrics['skill'].values)**2)

    # Try a more explicit amplitude and lag metric
    # mean square minutes error
    lag_cost=((metrics['lag'].values*24*60)**2).mean()
    # mean square pct error, normalized to 2%
    amp_cost=(((metrics['ratio'].values-1.0)/0.02)**2).mean()
    return lag_cost+amp_cost

##

def print_settings_and_score(settings,model):
    print()
    print(settings_to_txt(settings))
    print("SCORE: %.5f"%get_model_score(model))

def optimize():
    # the starting point
    best_settings=baseline_settings()

    for abs_delta in 0.005,0.0025:
        while 1:
            accepted_new=False
            for name in regions['name']:
                print("Optimizing %s, accepted_new=%s"%(name,accepted_new))
                # baseline score before changing this parameter
                model=check_and_run_settings(best_settings)
                # print_settings_and_score(best_settings,model)

                # best of this iteration
                best_score=init_score=get_model_score(model)

                # try increasing, then decreasing
                for delta_n in [abs_delta,-abs_delta]:
                    # try increasing this region's roughness
                    while 1:
                        settings=dict(best_settings) # copy
                        settings[name]+=delta_n # step
                        # these should be pulled from the shapefile, but for now
                        # just hardcode broad limits
                        if (settings[name]>0.040) or (settings[name]<0.015):
                            print(" [hit limit] ")
                            break

                        model=check_and_run_settings(settings)
                        # print_settings_and_score(settings,model)
                        this_score=get_model_score(model)
                        if this_score<best_score:
                            print('----New best----  %.5f => %.5f'%(best_score,this_score))
                            best_settings=settings
                            best_score=this_score
                            accepted_new=True
                            continue # try again
                        else:
                            break
            if not accepted_new:
                print("Converged at abs_delta=%.5f"%abs_delta)
                break


def print_best():
    metrics_fns=glob.glob("runs/roughsearch/*/DFM_OUTPUT_flowfm/metrics.csv")
    scores=[]
    for metrics_fn in metrics_fns:
        metrics=pd.read_csv(metrics_fn)
        score=metrics_to_score(metrics)
        scores.append(score)

    order=np.argsort(scores)
    for i in order:
        print("%s: %.5f"%(metrics_fns[i],scores[i]))


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Optimize roughness.')

    parser.add_argument('--optimize','-o',help='run optimization',action='store_true')
    parser.add_argument('--print','-p',help='print status',action='store_true')

    args = parser.parse_args()

    if args.optimize:
        optimize()
    elif args.print:
        print_best()
    else:
        parser.print_help()
