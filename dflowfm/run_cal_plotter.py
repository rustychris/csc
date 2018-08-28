import os
import six

from stompy import utils
utils.path("/home/rusty/src/hydro")

# from Plotting import filter_wrapper
# six.moves.reload_module(filter_wrapper)
import matplotlib.pyplot as plt
import Plotting.hydro_cal_plotter as hcp

six.moves.reload_module(hcp)

# globals
his_dir = 'runs/20180807_grid98_16_single/DFM_OUTPUT_flowfm'
his_file = os.path.join(his_dir,'flowfm_his.nc')

#his_dir = 'runs/base20180701/DFM_OUTPUT_FlowFM'
#his_file = os.path.join(his_dir,'FlowFM_his.nc')

plot_dir = os.path.join(his_dir,'figs-20180820/')
os.path.exists(plot_dir) or os.makedirs(plot_dir)

# Moving cal data into the repo, now reference relative directory
csv_dir="calibration_data"
defaults = {
    'corr_datum': True,
    'wl_range': [0.0, 2.5],
    # 'plot_type': 'plot_3panels',
    'plot_type': 'plot_basic',
    'station_dir': csv_dir,
    'plot_dir': plot_dir,
    'point_limit': 500,  # maximum number of points to show in scatter plots, or None for all points
}

stations = [
    ('DOP',{'var':'stage', 'units':'m',
            'obs_file':'HecDssExcel9146174433491355875_DOP_wsel.csv'}),
    ('CCS',{'var':'stage', 'units':'m',
            # RH: this file has an hour error, it appears
            #'obs_file':'HecDssExcel1634846602065331917_CCS_wsel.csv'
            # RH: use original file from A Rabidoux
            'obs_file':'CCS_orig_stage_m.csv'
    }),
    ('CCSC',{'var':'stage', 'units':'m',
             # 'obs_file':'HecDssExcel1634846602065331917_CCS_wsel.csv'
             # RH see above
             'obs_file':'CCS_orig_stage_m.csv'
    }),
    ('LSHB',{'var':'stage', 'units':'m',
             'obs_file':'LSHB_stage_m.csv'}),
    ('CCEH',{'var':'stage', 'units':'m',
             'obs_file':'CCEH_syn_stage_m.csv'}),
    
    ('HS1',{'var':'stage', 'units':'m',
            'obs_file':'HecDssExcel41820602247344004_HS1_wsel.csv'}),
    ('HAAS',{'var':'stage', 'units':'m',
             'obs_file':'HecDssExcel8019730502736234932_HAAS_wsel.csv'}),
    ('GES',{'var':'stage', 'units':'ft',
            'obs_file':'GES-2014-04-stage.csv'}),
    ('LN2',{'var':'stage', 'units':'m',
            'obs_file':'HecDssExcel6060592544570966047_LN2_wsel.csv'}),
    ('SG1',{'var':'stage', 'units':'m',
            'obs_file':'HecDssExcel7247236844972400088_SG1_wsel.csv'}),
    
    ('DWS',{'var':'stage', 'units':'ft', 'obs_file':'DWS-2014-04-stage.csv'}),
    ('FPX',{'var':'stage', 'units':'ft', 'obs_file':'FPX-2014-04-stage.csv'}),
    # transposed station name?
    ('HBW',{'var':'stage', 'units':'ft', 'obs_file':'HWB-2014-04-stage.csv'}),
    ('LIS',{'var':'stage', 'units':'ft', 'obs_file':'LIS-stage-WY2014.csv'}),
    ('RYI',{'var':'stage', 'units':'ft', 'obs_file':'RYI-2014-04-stage.csv'}),
    ('SRV',{'var':'stage', 'units':'ft', 'obs_file':'SRV-2014-04-stage.csv'}),
    ('SSS',{'var':'stage', 'units':'ft', 'obs_file':'SSS-2014-04-stage.csv'}),
    ('SUT',{'var':'stage', 'units':'ft', 'obs_file':'SUT-2014-04-stage.csv'}),
    
    # No data until July 7, 2014 according to the csv in calibration_data
    # ('UL1',{'var':'stage', 'units':'m',
    #         'obs_file':'HecDssExcel7632503257539881164_UL1_wsel.csv'}),
    
    ('LSHB',{'var':'flow','units':'cfs','pred_xs_name':'LSHB',
             'obs_file':'LSHB_flow_cfs.csv'}),
    
    ('GES',{'var':'flow','units':'cfs','pred_xs_name':'GES',
            'obs_file':'GES-2014-04-flow.csv'}),
    ('RYI',{'var':'flow','units':'cfs','pred_xs_name':'RYI',
            'obs_file':'RYI-2014-04-flow.csv'}),
    ('SRV',{'var':'flow','units':'cfs','pred_xs_name':'RioVista',
                 'obs_file':'SRV-2014-04-flow.csv'}),
    ('HWB',{'var':'flow','units':'cfs','pred_xs_name':'HWB',
            'obs_file':'HWB-2014-04-flow.csv'}),
    ('SSS',{'var':'flow','units':'cfs','pred_xs_name':'SSS',
            'obs_file':'SSS-2014-04-flow.csv'}),
    ('SUT',{'var':'flow','units':'cfs','pred_xs_name':'SUT',
            'obs_file':'SUT-2014-04-flow.csv'}),
    ('FPX',{'var':'flow','units':'cfs','pred_xs_name':'Freeport',
            'obs_file':'FPX-2014-04-flow.csv'}),

    ( 'DWS',{'var':'flow','units':'cfs','pred_xs_name':'DWS',
             'obs_file':'DWS-2014-04-flow.csv'}),

    # DOP stage is taken from Thomas's data, but that did not appear to
    # have flows -- flow is pulled from the txt in ARabidoux
    ('DOP',{'var':'flow', 'units':'cfs','pred_xs_name':'DOP',
            'obs_file':'DOP_flow_cfs.csv'}),

]

mr = hcp.ModelResults(his_file)
trim_spinup_days=1.0
dnums=utils.to_dnum(np.array(mr.htimes))
trim_count=np.searchsorted( dnums, dnums[0]+trim_spinup_days)
mr.htimes=mr.htimes[trim_count:]
mr.his=mr.his.isel(time=slice(trim_count,None))

class MyPlotter(hcp.HydroCalPlot):
    ax_txt=None # where we actually put the texts
    lines=None
    ax_txt_text=None

    def add_text_to_scatter(self,ax,string,value=None):
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

plt.ioff()
for station in stations:
    args = station[1]
    args['ID'] = station[0]
    for key in defaults:
        args.setdefault(key, defaults[key])  # append default parameters if missing
    p = MyPlotter(mr, args)
    if p.valid:
        p.plot()
plt.ion()
# plt.show()
