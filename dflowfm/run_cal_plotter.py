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
his_dir = 'runs/20180807_grid98_06/DFM_OUTPUT_flowfm'
his_file = os.path.join(his_dir,'flowfm_0000_his.nc')

#his_dir = 'runs/base20180701/DFM_OUTPUT_FlowFM'
#his_file = os.path.join(his_dir,'FlowFM_his.nc')

plot_dir = os.path.join(his_dir,'figs-20180820/')
os.path.exists(plot_dir) or os.makedirs(plot_dir)

# Moving cal data into the repo, now reference relative directory
csv_dir="calibration_data"
defaults = {
    'corr_datum': True,
    'wl_range': [0.0, 2.5],
    'plot_type': 'plot_3panels',
    'station_dir': csv_dir,
    'plot_dir': plot_dir,
    'point_limit': 500,  # maximum number of points to show in scatter plots, or None for all points
}

# additional stations:
# 'SOI', 'SXS', 'RYI', 'HBW', 'SSS', 'SUT', 'DWS', 'LIB', 'LIS',
# 'FPX', 'SRV'

# additional cross-sections:
#   RioVista Freeport RYI HWB SSS SUT GES

stations = [
    ('DOP',{'var':'stage', 'units':'m',
            'obs_file':'HecDssExcel9146174433491355875_DOP_wsel.csv'}),
    ('CCS',{'var':'stage', 'units':'m',
            'obs_file':'HecDssExcel1634846602065331917_CCS_wsel.csv'}),
    ('CCSC',{'var':'stage', 'units':'m',
             'obs_file':'HecDssExcel1634846602065331917_CCS_wsel.csv'}),
    ('LSHB',{'var':'stage', 'units':'m',
             'obs_file':'LSHB_stage_m.csv'}),

    # ('HS1',{'var':'stage', 'units':'m',
    #         'obs_file':'HecDssExcel41820602247344004_HS1_wsel.csv'}),
    # ('HAAS',{'var':'stage', 'units':'m',
    #          'obs_file':'HecDssExcel8019730502736234932_HAAS_wsel.csv'}),
    # ('GES',{'var':'stage', 'units':'ft',
    #         'obs_file':'GES-2014-04-stage.csv'}),
    ('LN2',{'var':'stage', 'units':'m',
            'obs_file':'HecDssExcel6060592544570966047_LN2_wsel.csv'}),
    # ('SG1',{'var':'stage', 'units':'m',
    #         'obs_file':'HecDssExcel7247236844972400088_SG1_wsel.csv'}),

    # # These seem to be okay re units and time zone
    # ('DWS',{'var':'stage', 'units':'ft', 'obs_file':'DWS-2014-04-stage.csv'}),
    # ('FPX',{'var':'stage', 'units':'ft', 'obs_file':'FPX-2014-04-stage.csv'}),
    # # transposed station name?
    # ('HBW',{'var':'stage', 'units':'ft', 'obs_file':'HWB-2014-04-stage.csv'}),
    # ('LIS',{'var':'stage', 'units':'ft', 'obs_file':'LIS-stage-WY2014.csv'}),
    # ('RYI',{'var':'stage', 'units':'ft', 'obs_file':'RYI-2014-04-stage.csv'}),
    # ('SRV',{'var':'stage', 'units':'ft', 'obs_file':'SRV-2014-04-stage.csv'}),
    # ('SSS',{'var':'stage', 'units':'ft', 'obs_file':'SSS-2014-04-stage.csv'}),
    # ('SUT',{'var':'stage', 'units':'ft', 'obs_file':'SUT-2014-04-stage.csv'}),

    # No data until July 7, 2014 according to the csv in calibration_data
    # ('UL1',{'var':'stage', 'units':'m',
    #         'obs_file':'HecDssExcel7632503257539881164_UL1_wsel.csv'}),

    # temporarily use -cfs to indicate flipped sign.
    ('LSHB',{'var':'flow','units':'cfs','pred_xs_name':'LSHB',
             'obs_file':'LSHB_flow_cfs.csv'}),

    # ('GES',{'var':'flow','units':'cfs','pred_xs_name':'GES',
    #         'obs_file':'GES-2014-04-flow.csv'}),
    # ('RYI',{'var':'flow','units':'cfs','pred_xs_name':'RYI',
    #         'obs_file':'RYI-2014-04-flow.csv'}),
    # ('SRV',{'var':'flow','units':'cfs','pred_xs_name':'RioVista',
    #              'obs_file':'SRV-2014-04-flow.csv'}),
    # ('HWB',{'var':'flow','units':'cfs','pred_xs_name':'HWB',
    #         'obs_file':'HWB-2014-04-flow.csv'}),
    # ('SSS',{'var':'flow','units':'cfs','pred_xs_name':'SSS',
    #         'obs_file':'SSS-2014-04-flow.csv'}),
    # ('SUT',{'var':'flow','units':'cfs','pred_xs_name':'SUT',
    #         'obs_file':'SUT-2014-04-flow.csv'}),
    # ('FPX',{'var':'flow','units':'cfs','pred_xs_name':'Freeport',
    #         'obs_file':'FPX-2014-04-flow.csv'}),

    #  DWS: not currently in the model as a flow section

    # DOP stage is taken from Thomas's data, but that did not appear to
    # have flows -- flow is pulled from the txt in ARabidoux
    ('DOP',{'var':'flow', 'units':'cfs','pred_xs_name':'DOP',
            'obs_file':'DOP_flow_cfs.csv'}),

]

mr = hcp.ModelResults(his_file)
plt.ioff()
for station in stations:
    args = station[1]
    args['ID'] = station[0]
    for key in defaults:
        args.setdefault(key, defaults[key])  # append default parameters if missing
    p = hcp.HydroCalPlot(mr, args)
    if p.valid:
        p.plot()

plt.ion()
