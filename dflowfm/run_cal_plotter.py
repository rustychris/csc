import os
import pandas as pd
import xarray as xr
import six
import logging
import numpy as np
import stompy.model.delft.dflow_model as dfm
from stompy.plot import plot_wkb, plot_utils
from stompy import utils
import stompy.grid.unstructured_grid as ugrid

# kludge because I don't keep this in my PYTHONPATH
utils.path("/home/rusty/src/hydro")

import matplotlib.pyplot as plt
import hydro.plotting.hydro_cal_plotter as hcp

plot_dir="NOT_SET_YET" # gets set at the end.

try:
    base_dir=os.path.dirname(__file__)
except NameError:
    logging.info("Assuming current directory is the dflowfm base")
    base_dir="."
    
# Moving cal data into the repo, now reference relative directory
csv_dir=os.path.join(base_dir,"calibration_data")
def csv(fn): return os.path.join(csv_dir,fn)

defaults = {
    'corr_datum': True,
    'wl_range': [0.0, 2.5],
    'plot_type': 'plot_3panels', # 'plot_basic'
    'station_dir': csv_dir,
    # 'plot_dir': plot_dir, # set dynamically below
    'point_limit': 500,  # maximum number of points to show in scatter plots, or None for all points
}

stations = [
    ('CourtlandToe',dict(var='flow',units='cfs',filename=csv('CourtToe-flow.csv'))),
    ('CourtlandToe',dict(var='stage',units='ft',filename=csv('CourtToe-stage.csv'))),

    ('HollandNBreach',dict(var='flow',units='cfs',filename=csv("HNB-flow.csv"))),
    ('HollandNBreach',dict(var='stage',units='ft',filename=csv("HNB-stage.csv"))),
    
    ("ToeAtLiberty", dict(var='flow',units='cfs',filename=csv("LibertyToe-flow.csv"))),
    ("ToeAtLiberty", dict(var='stage',units='ft',filename=csv("LibertyToe-stage.csv"))),

    ("LibertyIslandCut", dict(var='flow',units='cfs',filename=csv("LIC-flow.csv"))),
    ("LibertyIslandCut", dict(var='stage',units='ft',filename=csv("LIC-stage.csv"))),

    # is this the correct one, rather than LiberyIslandCut?
    ("Stairstep", dict(var='flow',units='cfs',filename=csv("LIC-flow.csv"))),

    ("LibCutHolland", dict(var='flow',units='cfs',filename=csv("LIH-flow.csv"))),

    ("WildlandsUpMarsh", dict(var='flow',units='cfs',filename=csv("WildlandsUpMarsh-flow.csv"))),

    ("LIY", dict(var='stage',units='ft',filename=csv('LIY-stage.csv'))),

    ("LIB", dict(var='flow',units='cfs',filename=csv("LIB-flow.csv"))),
    ("LIB", dict(var='stage',units='ft',filename=csv("LIB-stage.csv"))),

    ("MIR", dict(var='flow',units='cfs',filename=csv("MIR-flow.csv"))),
    ("MIR", dict(var='stage',units='ft',filename=csv("MIR-stage.csv"))),
    
    ('TSL',{'var':'flow','units':'cfs',
            'filename':csv('TSL-flow.csv')}),
    ('TSL',{'var':'stage','units':'ft',
            'filename':csv('TSL-stage.csv')}),
    ('SDI',{'var':'stage','units':'ft',
            'filename':csv('SDI-stage.csv')}),
    ('DLC',{'var':'flow','units':'cfs','pred_xs_name':'DLC',
            'filename':csv('DCC-flow.csv')}),
    ('GSS',{'var':'flow','units':'cfs','pred_xs_name':'Georgiana_xc',
            'filename':csv('GSS-flow.csv')}),
    ('SDC',{'var':'flow','units':'cfs','pred_xs_name':'SDC',
            'filename':csv('SDC-flow.csv')}),
    
    ('DOP',{'var':'stage', 'units':'m',
            'filename':csv('HecDssExcel9146174433491355875_DOP_wsel.csv')}),
    ('CCS',{'var':'stage', 'units':'m',
            # RH: this file has an hour error, it appears
            #'filename':'HecDssExcel1634846602065331917_CCS_wsel.csv'
            # RH: use original file from A Rabidoux
            'filename':csv('CCS_orig_stage_m.csv')
    }),
    ('CCSC',{'var':'stage', 'units':'m',
             # 'filename':'HecDssExcel1634846602065331917_CCS_wsel.csv'
             # RH see above
             'filename':csv('CCS_orig_stage_m.csv')
    }),
    ('LSHB',{'var':'stage', 'units':'m',
             'filename':csv('LSHB_stage_m.csv')}),
    ('CCEH',{'var':'stage', 'units':'m',
             'filename':csv('CCEH_syn_stage_m.csv')}),
    
    ('HS1',{'var':'stage', 'units':'m',
            'filename':csv('HecDssExcel41820602247344004_HS1_wsel.csv')}),
    ('HAAS',{'var':'stage', 'units':'m',
             'filename':csv('HecDssExcel8019730502736234932_HAAS_wsel.csv')}),
    ('GES',{'var':'stage', 'units':'ft',
            'filename':csv('GES-stage.csv')}),
    ('LN2',{'var':'stage', 'units':'m',
            'filename':csv('HecDssExcel6060592544570966047_LN2_wsel.csv')}),
    ('SG1',{'var':'stage', 'units':'m',
            'filename':csv('HecDssExcel7247236844972400088_SG1_wsel.csv')}),
    
    ('DWS',{'var':'stage', 'units':'ft', 'filename':csv('DWS-stage.csv')}),
    ('FPX',{'var':'stage', 'units':'ft', 'filename':csv('FPX-stage.csv')}),
    # transposed station name?
    ('HBW',{'var':'stage', 'units':'ft', 'filename':csv('HWB-stage.csv')}),

    ('LIS',{'var':'flow', 'units':'cfs', 'filename':csv('LIS-flow.csv')}),

    ('RYI',{'var':'stage', 'units':'ft', 'filename':csv('RYI-stage.csv')}),
    ('SRV',{'var':'stage', 'units':'ft', 'filename':csv('SRV-stage.csv')}),
    ('SSS',{'var':'stage', 'units':'ft', 'filename':csv('SSS-stage.csv')}),
    ('SUT',{'var':'stage', 'units':'ft', 'filename':csv('SUT-stage.csv')}),
    
    # No data until July 7, 2014 according to the csv in calibration_data
    # ('UL1',{'var':'stage', 'units':'m',
    #         'filename':'HecDssExcel7632503257539881164_UL1_wsel.csv'}),
    
    ('LSHB',{'var':'flow','units':'cfs','pred_xs_name':'LSHB',
             'filename':csv('LSHB_flow_cfs.csv')}),
    ('GES',{'var':'flow','units':'cfs','pred_xs_name':'GES',
            'filename':csv('GES-flow.csv')}),
    ('RYI',{'var':'flow','units':'cfs','pred_xs_name':'RYI',
            'filename':csv('RYI-flow.csv')}),
    ('RioVista',{'var':'flow','units':'cfs','pred_xs_name':'RioVista',
            'filename':csv('SRV-flow.csv')}),
    ('HWB',{'var':'flow','units':'cfs','pred_xs_name':'HWB',
            'filename':csv('HWB-flow.csv')}),
    ('SSS',{'var':'flow','units':'cfs','pred_xs_name':'SSS',
            'filename':csv('SSS-flow.csv')}),
    ('SUT',{'var':'flow','units':'cfs','pred_xs_name':'SUT',
            'filename':csv('SUT-flow.csv')}),
    ('Freeport',{'var':'flow','units':'cfs','pred_xs_name':'Freeport',
            'filename':csv('FPX-flow.csv')}),
    ('DWS',{'var':'flow','units':'cfs','pred_xs_name':'DWS',
            'filename':csv('DWS-flow.csv')}),

    # DOP stage is taken from Thomas's data, but that did not appear to
    # have flows -- flow is pulled from the txt in ARabidoux
    ('DOP',{'var':'flow', 'units':'cfs','pred_xs_name':'DOP',
            'filename':csv('DOP_flow_cfs.csv')}),
]



class MyPlotter(hcp.HydroCalPlot):
    ax_txt=None # where we actually put the texts
    lines=None
    ax_txt_text=None

    def add_text_to_scatter(self,ax,string,value=None):
        """
        Monkey patch text formatting that is a bit more flexible
        """ 
        if self.ax_txt is None:
            return super(MyPlotter,self).add_text_to_scatter(ax,string,value=value)
            
        if value is not None:
            string="%s: %.3f"%(string,value)
        if self.lines is None:
            self.lines=[]
        self.lines.append(string)
        txt="\n".join(self.lines)
        if self.ax_txt_text is None:
            self.ax_txt_text=self.ax_txt.text(0.05,0.98,txt,ha='left',va='top')
        else:
            self.ax_txt_text.set_text(txt)

    def plot_basic(self):
        fig=plt.figure(figsize=[12,8])
        ax0 = fig.add_subplot(2,1,1)
        ax1 = fig.add_subplot(2,2,3)
        self.ax_txt = ax2 = fig.add_subplot(2,2,4)
        ax=[ax0,ax1]

        self.plot_inst(ax[0])
        self.plot_scatter(ax[1], self.args['point_limit'])
        figname = self.plot_dir + '/%s_%s_inst_scat.png' % (self.ID, self.var)
        fig.suptitle(self.ID)
        ax[0].xaxis.set_minor_locator(hcp.mdates.DayLocator(interval=1))
        ax[0].grid(which='minor', axis='x', linestyle="--")
        self.ax_txt.xaxis.set_visible(0)
        self.ax_txt.yaxis.set_visible(0)
        plt.setp( list(self.ax_txt.spines.values()), visible=0)
        plt.savefig(figname, dpi=600)

##
# Develop energy flux plots

def trim_time(mr,da):
    # trim to match model results, and fill in gaps with tidal data.
    start_time=utils.to_dt64(mr.htimes[0])
    stop_time=utils.to_dt64(mr.htimes[-1])
    return da.isel(time=(da.time>=start_time) & (da.time<=stop_time))

def name_to_da(mr,name,var):
    for stn_name,details in stations:
        if name==stn_name and details['var']==var:
            break
    else:
        return None # wasn't found
    df=pd.read_csv(details['filename'],parse_dates=['Time'])
    df=df.set_index('Time')
    ds=xr.Dataset.from_dataframe(df)
    if details['units']=='cfs':
        print("Adjusting cfs to cms")
        ds.Flow.values *= 0.028316847
        ds.Flow.attrs['units']='m3 s-1'
    elif details['units']=='ft':
        ds.Stage.values *= 0.3048
        ds.Stage.attrs['units']='m'

    ds=ds.rename(Time='time')
    if var=='flow':
        da=ds['Flow']
    elif var=='stage':
        da=ds['Stage']
    else:
        raise Exception("Not sure what the data variable is")
    da=trim_time(mr,da)
    da=utils.fill_tidal_data(da)

    return da
## 
def name_to_model_da(mr,name,var):
    if 1:
        for stn_name,details in stations:
            if name==stn_name and details['var']==var:
                break
        else:
            #return None # wasn't found
            raise Exception() # DEV
    model_id=details.get('pred_xs_name',name)

    # array of values
    data=mr.get_pred_ts(model_id,var)
    ds=xr.Dataset()
    ds['time']=('time',),mr.htimes
    ds[var]=('time',),data

    return ds[var]

## get the info for each:
# names in the station list
# all are negative, consistent with the upstream flux of tidal energy.
# Most are quite close, with the exception of DWS, justifying the
# decrease in friction there.

def flux_figure(mr,flow_station_name,stage_station_name):
    fig=plt.figure(1)
    fig.clf()
    fig.set_size_inches((7,9), forward=True)
    fig,axs=plt.subplots(3,1,num=1,sharex=True)

    for src in ['obs','mod']:
        if src=='obs':
            flow_da=name_to_da(mr,flow_station_name,'flow')
            stage_da=name_to_da(mr,stage_station_name,'stage')
        else:
            flow_da=name_to_model_da(mr,flow_station_name,'flow')
            stage_da=name_to_model_da(mr,stage_station_name,'stage')

        # E=h*Q ish?

        # range of valid overlapping data:
        t0=max(flow_da.time.values.min(), stage_da.time.values.min() )
        tN=min(flow_da.time.values.max(), stage_da.time.values.max() )
        common_time=np.arange(t0,tN,np.timedelta64(15*60,'s'))

        Q_int = np.interp( utils.to_dnum(common_time),
                           utils.to_dnum(flow_da.time.values),flow_da.values)
        h_int = np.interp( utils.to_dnum(common_time),
                           utils.to_dnum(stage_da.time.values), stage_da.values)

        # just tides:
        from stompy import filters

        Q_int_bar=filters.lowpass( Q_int, cutoff=40, dt=0.25 )
        h_int_bar=filters.lowpass( h_int, cutoff=40, dt=0.25 )
        Q_tidal=Q_int - Q_int_bar
        h_tidal=h_int - h_int_bar

        Qh=filters.lowpass(Q_tidal*h_tidal, cutoff=40, dt=0.25)

        pad=np.timedelta64(40*3600,'s')
        sel=(common_time>common_time[0]+pad) & (common_time<common_time[-1]-pad)

        axs[0].plot(flow_da.time,flow_da)
        axs[0].plot(common_time[sel],Q_tidal[sel],label='Q_tidal %s'%src)
        axs[1].plot(stage_da.time,stage_da)
        axs[1].plot(common_time[sel],h_tidal[sel],label='h_tidal %s'%src)
        axs[2].plot(common_time[sel], Qh[sel],label='Qh %s'%src)

        for ax in axs:
            ax.legend()
    xxyy=axs[2].axis()
    axs[2].axis(ymax=max(xxyy[3],0))
        
    axs[0].set_title(flow_station_name)
    return fig

# GES: -15    -14
# RYI: -250  -200
# FPX: -7.5    -7
# SRV: -400  -400
# SSS: -6      -7
# SUT: -3      -4
# HWB: -6      -9
# DWS: -30    -50

flow_stage_stations=[ ['GES','GES'],
                      ['RYI','RYI'],
                      ['Freeport','FPX'],
                      ['RioVista','SRV'],
                      ['SSS','SSS'],
                      ['SUT','SUT'],
                      ['HWB','HBW'],
                      ['DWS','DWS'] ]

def make_flux_figures(model,mr):
    for flow_station_name,stage_station_name in flow_stage_stations:
        fig=flux_figure(mr,flow_station_name,stage_station_name)
        fig.savefig(os.path.join(plot_dir,'energy_flux-%s.png'%flow_station_name))

##

# Summary map figures
def make_summary_map(model,mr):
    poly=model.grid.boundary_polygon()
    his_file=model.his_output()
    fig=plt.figure(10)
    fig.set_size_inches((6,9),forward=True)
    fig.clf()
    ax=fig.add_axes([0,0,1,1])
    ax.xaxis.set_visible(0)
    ax.yaxis.set_visible(0)

    plot_wkb.plot_wkb(poly,ax=ax,fc='0.75',ec='0.75',lw=0.8)
    ax.axis('equal')

    for station in stations:
        args = station[1]
        args['ID'] = station[0]
        for key in defaults:
            args.setdefault(key, defaults[key])  # append default parameters if missing

        p = MyPlotter(mr, args)
        if not p.valid: continue

        if args['ID']=='TSL' and args['var']=='flow' and 'grid100_13' in his_file:
            print("Monkey patch sign on TSL flow")
            p.pred *= -1
        metrics=p.calculate_metrics()
        if metrics['ratio']<-1000: continue # probably no data

        xy=p.xy
        if xy.ndim==2:
            xy=xy.mean(axis=0)
        # a little tricky with flow+stage
        if args['var']=='flow':
            kw=dict(va='top')
        else:
            kw=dict(va='bottom')
        kw['clip_on']=True
        kw['fontsize']=8
        ax.text(xy[0],xy[1],"%s %s: %.2f/%.1fmin"%(args['ID'],args['var'],
                                                   metrics['ratio'],
                                                   24*60*metrics['lag']) ,
                **kw)

    ax.axis( (605686., 639321., 4211471, 4267529))
    plot_utils.reduce_text_overlap(ax)
    fig.savefig(os.path.join(plot_dir,'amp_phase_map.png'),dpi=120)
    return fig

##

def calc_phases(model,
                period_s=12.4206*3600,
                ref_phase_loc=[610031., 4216014.]):

    map_fns=model.map_outputs()

    grids=[]
    results={}
    ref_eta_phase=0.0
    ref_u_phase=0.0
    ref_v_phase=0.0
    ref_uf_phase=0.0

    for proc,map_fn in enumerate(map_fns):
        print("-----------------%3d------------------"%proc)
        pmap=xr.open_dataset(map_fn)

        # short spinup of 1 day
        trim=pmap.time.values>pmap.time.values[0]+np.timedelta64(86400,'s')

        pmap_sel=pmap.isel(time=trim)
        proc_data={}
        results[proc]=proc_data
        gsub=ugrid.UnstructuredGrid.read_dfm(pmap)
        proc_data['grid']=gsub

        eta=pmap_sel.mesh2d_s1.values   # time,cell
        u = pmap_sel.mesh2d_ucx.values  # time,cell
        v = pmap_sel.mesh2d_ucy.values  # time,cell

        # Combined phase - first get flood-directed velocity
        u_flood=np.zeros(u.shape, np.float64)
        for c in utils.progress(range(gsub.Ncells())):
            c_uv=np.c_[u[:,c],v[:,c]]
            c_eta=eta[:,c]
            c_rot=utils.rotate_to_principal(c_uv,c_eta,ambiguous='silent',detrend=True)
            u_flood[:,c]=c_rot[:,0]

        Ncells=gsub.Ncells()
        # fit M2 tide only
        t_s=(pmap_sel.time.values - pmap_sel.time.values[0])/np.timedelta64(1,'s')
        fit_cos=np.cos(t_s/period_s*2*np.pi)
        fit_sin=np.sin(t_s/period_s*2*np.pi)

        def phase_map(val,fit_cos=fit_cos,fit_sin=fit_sin,small=1e-6):
            phases=np.zeros(Ncells,np.float64)

            for c in utils.progress(range(Ncells)):
                c_val = val[:,c]
                C=np.cov( np.c_[fit_cos,fit_sin,c_val].T )

                if C[2,2] > small: # some signal at this cell
                    phase=np.arctan2(C[0,2],C[1,2]) * 180/np.pi
                else:
                    phase=np.nan
                phases[c]=phase
            return phases

        proc_data['eta_phases']=phase_map(eta)
        #proc_data['u_phases']  =phase_map(u)
        #proc_data['v_phases']  =phase_map(v)
        proc_data['uf_phases'] =phase_map(u_flood)

        c_ref=gsub.select_cells_nearest(ref_phase_loc,inside=True)
        if c_ref is not None:
            results['ref_eta_phase']=proc_data['eta_phases'][c_ref]
            #results['ref_u_phase']=u_phases[c_ref]
            #results['ref_v_phase']=v_phases[c_ref]
            results['ref_uf_phase']=proc_data['uf_phases'][c_ref]
    return results

##

def make_phase_figure(model,mr):
    poly=model.grid.boundary_polygon()
    
    period_s=12.4206*3600
    results=calc_phases(model,period_s=period_s)
    
    fig=plt.figure(11)
    fig.clf()
    fig.set_size_inches((10,10),forward=True)
    fig,axs=plt.subplots(1,2,sharex=True,sharey=True,num=11)
    axs[0].set_position([0,0,0.5,1.0])
    axs[1].set_position([0.5,0,0.5,1.0])

    caxs=[ fig.add_axes([0.05,0.5,0.03,0.45]),
           fig.add_axes([0.55,0.5,0.03,0.45]) ]

    ccolls=[]
    for proc in range(model.num_procs):
        proc_data=results[proc]
        for ax,cax,ph,name in zip(axs,caxs,
                                  [proc_data['eta_phases']-results['ref_eta_phase'],
                                   proc_data['uf_phases']-results['ref_uf_phase']],
                                  ['eta','u']):
            gsub=proc_data['grid']

            phase_relative=(ph + 180)%360 - 180
            phase_minutes=-phase_relative/360. * (period_s/60.)

            period_minutes=period_s/60.
            ccoll=gsub.plot_cells(values=phase_minutes%period_minutes,
                                  clim=[0,period_minutes],ax=ax)
            if 0: # contours are messy
                ecoll=gsub.scalar_contour(phase_minutes,V=np.arange(-30,500,15))
                ax.add_collection(ecoll)
            ccoll.set_cmap('hsv')

            if proc==0:
                ax.xaxis.set_visible(0)
                ax.yaxis.set_visible(0)
                plot_wkb.plot_wkb(poly,ax=ax,fc='none',ec='0.5',lw=0.8,zorder=-10)
                plt.colorbar(ccoll,orientation='vertical',cax=cax,
                             label='%s phase (minutes)'%name)

    ax.axis('equal')
    fig.savefig(os.path.join(plot_dir,"phases.png"))


def make_time_series_figures(model,mr):
    for station in stations:
        args = station[1]
        args['ID'] = station[0]
        for key in defaults:
            args.setdefault(key, defaults[key])  # append default parameters if missing
        p = MyPlotter(mr, args)
        if p.valid:
            p.plot()
    
##

fig_version="20190430"

if __name__=='__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Plot calibration figures for Delta DFM models.')

    parser.add_argument('mdupath',help='path to mdu file')

    args = parser.parse_args()
    print(args.mdupath)

    model=dfm.DFlowModel.load(args.mdupath)
    
    his_file=model.his_output()
    his_dir=os.path.dirname(his_file)
    plot_dir = os.path.join(his_dir,'figs-%s'%fig_version)
    os.path.exists(plot_dir) or os.makedirs(plot_dir)
    
    mr = hcp.DflowfmModelResults([his_file],trim_spinup_days=1.0)

    # A bit awkward to mutate stations...
    for station in stations:
        station[1]['plot_dir']=plot_dir

    plt.ioff()
    make_summary_map(model,mr)
    make_time_series_figures(model,mr)
    make_flux_figures(model,mr)
    make_phase_figure(model,mr)
    plt.ion()
    
