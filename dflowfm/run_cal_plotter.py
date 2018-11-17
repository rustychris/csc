import os
import six
import numpy as np

from stompy import utils
utils.path("/home/rusty/src/hydro")

import matplotlib.pyplot as plt
import Plotting.hydro_cal_plotter as hcp

six.moves.reload_module(hcp)

# globals
his_dir = 'runs/grid100_00/DFM_OUTPUT_flowfm'
his_file = os.path.join(his_dir,'flowfm_0000_his.nc')

#his_dir = 'runs/base20180701/DFM_OUTPUT_FlowFM'
#his_file = os.path.join(his_dir,'FlowFM_his.nc')

plot_dir = os.path.join(his_dir,'figs-20181116')
os.path.exists(plot_dir) or os.makedirs(plot_dir)

# Moving cal data into the repo, now reference relative directory
csv_dir="calibration_data"
def csv(fn): return os.path.join(csv_dir,fn)

defaults = {
    'corr_datum': True,
    'wl_range': [0.0, 2.5],
    'plot_type': 'plot_3panels',
    # 'plot_type': 'plot_basic',
    'station_dir': csv_dir,
    'plot_dir': plot_dir,
    'point_limit': 500,  # maximum number of points to show in scatter plots, or None for all points
}

stations = [
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
            'filename':csv('GES-2014-04-stage.csv')}),
    ('LN2',{'var':'stage', 'units':'m',
            'filename':csv('HecDssExcel6060592544570966047_LN2_wsel.csv')}),
    ('SG1',{'var':'stage', 'units':'m',
            'filename':csv('HecDssExcel7247236844972400088_SG1_wsel.csv')}),
    
    ('DWS',{'var':'stage', 'units':'ft', 'filename':csv('DWS-stage.csv')}),
    ('FPX',{'var':'stage', 'units':'ft', 'filename':csv('FPX-stage.csv')}),
    # transposed station name?
    ('HBW',{'var':'stage', 'units':'ft', 'filename':csv('HWB-stage.csv')}),
    ('LIS',{'var':'stage', 'units':'ft', 'filename':csv('LIS-stage-WY2014.csv')}),
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
    ('SRV',{'var':'flow','units':'cfs','pred_xs_name':'RioVista',
            'filename':csv('SRV-flow.csv')}),
    ('HWB',{'var':'flow','units':'cfs','pred_xs_name':'HWB',
            'filename':csv('HWB-flow.csv')}),
    ('SSS',{'var':'flow','units':'cfs','pred_xs_name':'SSS',
            'filename':csv('SSS-flow.csv')}),
    ('SUT',{'var':'flow','units':'cfs','pred_xs_name':'SUT',
            'filename':csv('SUT-flow.csv')}),
    ('FPX',{'var':'flow','units':'cfs','pred_xs_name':'Freeport',
            'filename':csv('FPX-flow.csv')}),

    ('DWS',{'var':'flow','units':'cfs','pred_xs_name':'DWS',
            'filename':csv('DWS-flow.csv')}),

    # DOP stage is taken from Thomas's data, but that did not appear to
    # have flows -- flow is pulled from the txt in ARabidoux
    ('DOP',{'var':'flow', 'units':'cfs','pred_xs_name':'DOP',
            'filename':csv('DOP_flow_cfs.csv')}),

]
mr = hcp.DflowfmModelResults([his_file])
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

plt.ioff()
for station in stations:
    if station[0]!='SSS': continue # DBG

    args = station[1]
    args['ID'] = station[0]
    for key in defaults:
        args.setdefault(key, defaults[key])  # append default parameters if missing
    p = MyPlotter(mr, args)
    if p.valid:
        p.plot()
plt.ion()
# plt.show()
