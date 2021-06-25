import six
import os
import numpy as np
import numpy.fft as fft
import pandas as pd
from pandas.plotting import register_matplotlib_converters;

register_matplotlib_converters()
import datetime as dt
import xarray as xr
from matplotlib import ticker

import seaborn as sb
from scipy.stats import pearsonr
import statsmodels.api as sm
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from itertools import combinations

from stompy import utils, filters
from stompy.model import data_comparison
import stompy.model.hydro_model as hm
from stompy.model.delft import dflow_model

# define relative path to model results here
model = dflow_model.DFlowModel.load('runs/sb_rv4_20190601')
# check against cal_figures to see if offset is necessary (not used for decker BC runs)
model.utc_offset = np.timedelta64(-8, 'h')  # model runs in PST
# time to cut from beginning of model
spin_up = 14  # [d]
# cutoff for lowpass filtering of data
cutoff = 3. / 24  # [d]

# this station only has NO3+NO2 starting from mid-June 2020
# dict(sources=[hm.NwisScalarBC(name='LIS', station=11455095, scalar='NO3+NO2',
#                              cache_dir='cache', model=model)],
#     station_name='Deep Water Ship Channel near Freeport')

# sorted from upstream to downstream (roughly), so then when pairing stations 1 is always upstream of 2
# format = {'model PLI codename': (USGS station no., 'verbose name')}
bgc_stations = {'FPX': (11447650, 'Sac R. at Freeport'),
                'SDC': (11447890, 'Sac R. at Delta Cross-Channel'),
                'SRV': (11455478, 'Sac R. at Decker Island (model age at Rio Vista)'),
                'RYI': (11455385, 'Cache at Ryer Island'),
                'LIB': (11455315, 'Cache at Liberty'),
                # 'ToeAtLiberty':(11455140, 'Toe Drain at Liberty')  # this one doesn't have nit?
                }

# BGC params to extract from each station
scalars = ['NO3+NO2',
           'temperature',
           'pH',
           'fDOM']

# names used for plot labeling
aliases = {'frac_1': 'Sac. fraction 1',
           'frac_2': 'Sac. fraction 2',
           'd_age': r'$\Delta$age [d]',
           'temperature_avg': r'$\langle$T$\rangle$ [$^{\circ}$C]',
           'pH_avg': r'$\langle$pH$\rangle$',
           'NO3+NO2_delta': r'$\Delta$(NO$_3$+NO$_2$) [mg/l]',
           'NO3+NO2_delta_frac': r'$\Delta$(NO$_3$+NO$_2$)/(NO$_3$+NO$_2$)',
           'fDOM_delta': r'$\Delta$fDOM [$\mu$g/l QSE]',
           'fDOM_delta_frac': r'$\Delta$fDOM/fDOM'}

# use BGC station info to create a set of NwisScalarBC objects for model/gage data retrieval
bgc_gages = [dict(sources=[hm.NwisScalarBC(name=st,
                                           station=bgc_stations[st][0],
                                           scalar=s,
                                           cache_dir='cache',
                                           model=model) for s in scalars],
                  station_name=bgc_stations[st]) for st in bgc_stations.keys()]

# copied this from cal_figures.py
defaults = dict(zoom_period=[model.run_stop - np.timedelta64(5, 'D'), model.run_stop],
                models=[model])

station_names = []
dats = []

for bgc_gage in bgc_gages:
    name = bgc_gage['station_name']
    print(f'getting {name }...')
    settings = dict(defaults)  # copy
    settings.update(bgc_gage)

    # get model/gage data at the station
    sources, combined = data_comparison.assemble_comparison_data(settings['models'],
                                                                 settings['sources'])

    # calculate age from model outputs
    conc = sources[0]  # ZNit
    age_conc = sources[1]  # NO3 (model age tracer, not actual nit)
    age = age_conc / conc
    valid = conc > 0.0001
    age[~valid] = np.nan

    # organize data for this station in a dict
    colnames = ['age'] + scalars
    data = [age] + sources[2:]
    dat = dict(zip(colnames, data))

    station_names.append(name)
    dats.append(dat)

# pair up stations for comparison
age_nit_tuples = zip(station_names, dats)
station_pairs = combinations(age_nit_tuples, 2)

# set output path (subdir of model run 'DFM_OUTPUT_flowfm' directory)
fig_dir = os.path.join(os.path.dirname(model.his_output()), "bgc_age_corrs")
os.path.exists(fig_dir) or os.mkdir(fig_dir)

# iterate over all station pairs
for station_pair in station_pairs:
    (name1, dat1), (name2, dat2) = station_pair
    print(name1[1], name2[1])

    # fraction of Sac River water at upstream and downstream station
    # frac1, frac2 = dat1['frac'], dat2['frac']

    # calculate age difference (d_age)
    age1, age2 = dat1['age'], dat2['age']
    d_age = age2 - age1
    d_age = d_age.dropna('time')
    # lowpass d_age
    d_age.values = filters.lowpass(d_age.values, utils.to_dnum(d_age.time), cutoff=cutoff)
    t0 = d_age.time[0]  # start time
    # remove starting spin-up time from d_age
    d_age = d_age[d_age.time > t0 + pd.to_timedelta(spin_up, 'd')]
    # calculate times to grab second station nitrate/BGC values (offset by d_age)
    t2s = [t + pd.to_timedelta(da, 'd') for da, t in zip(d_age.values, d_age.time.values)]

    # dataframe for joining all observations and making correlograms, etc.
    df = pd.DataFrame(data={# 'frac_1': frac1.interp(time=d_age.time).values,
                            # 'frac_2': frac2.interp(time=d_age.time).values,
                            'd_age': d_age.values},
                      index=d_age.time.values)

    # calculate averages and differences of scalars between stations
    for scalar in scalars:
        d1, d2 = dat1[scalar], dat2[scalar]
        # offsetting second station measurements by d_age
        avg = (d1.interp(time=d_age.time) + d2.interp(time=t2s).values) / 2
        delta = d2.interp(time=t2s).values - d1.interp(time=d_age.time)
        delta = delta.dropna('time')
        if len(delta) < 10:
            print(f'ERROR: missing d_{scalar} data.')
            continue
        # fractional difference
        delta_frac = delta / d1.interp(time=delta.time)

        # remove starting spin-up time
        avg = avg[avg.time > t0 + pd.to_timedelta(spin_up, 'd')]
        delta = delta[delta.time > t0 + pd.to_timedelta(spin_up, 'd')]

        # ensure there is enough data
        if len(d_age) < 10 or len(avg) < 10 or len(delta) < 10:
            print(f'Not enough valid data for pair {name1} + {name2}')
            continue

        # look at fractional deltas for nitrate, fDOM
        if scalar in ['NO3+NO2', 'fDOM']:
            temp_df = pd.DataFrame(data={f'{scalar}_delta_frac': delta_frac.values}, index=delta_frac.time.values)
        # look at averages for the other params (temperature, pH)
        else:
            temp_df = pd.DataFrame(data={f'{scalar}_avg': avg.values}, index=avg.time.values)
        df = df.join(temp_df)

    # only keep 15 min data (so we aren't interpolating 15 min USGS data to 5 min resolution)
    df = df[df.index.minute % 15 == 0]
    # interpolate gaps before lowpassing (linear backfilling), and then dropping NAs at end of data
    df = df.interpolate(method='linear', limit_direction='backward').dropna()

    # lowpass everything (except d_age which has already been lowpassed)
    for col in df.columns[1:]:
        df[col] = filters.lowpass(df[col], utils.to_dnum(df.index), cutoff=3./24)

    # identify dominant period of each signal
    periods_table = [[' ', '$T_1$ [hr]', '$|\^f|$', '$T_2$ [hr]', '$|\^f|$']]  # table header
    for col in df.columns:
        # check for gaps (must be regularly spaced data for FFT)
        i = 0
        for t1, t2 in zip(df.index, df.index[1:]):
            if (t2-t1).seconds != 900:
                print(t1, t2, (t2 - t1).seconds/900)
                i += 1
        print(f'data gaps = {i}')
        spectrum = fft.fft(df[col])
        freq = fft.fftfreq(len(spectrum))
        T = 1 / freq * 15 / 60  # [hrs]
        # get frequencies with greatest spectral magnitude corresponding to period < 100 hrs
        # (< 100 hrs to ignore more long-term/seasonal trends, which may also be spurious)
        peak_T = T[abs(spectrum) == max(abs(spectrum[abs(T) < 100]))][0]
        print(f' dominant period = {abs(peak_T):.2f} hrs')

        # get 6 periods with greatest spectral magnitudes (periods < 100 hrs)
        mask = abs(spectrum) > sorted(abs(spectrum[(T > 0) & (T < 100)]))[-6]
        peak_Ts = T[mask]
        peak_spectrum = spectrum[mask]
        # format periods, spectral magnitudes into table rows
        rows = [[round(t, 2), int(abs(spec))] for t, spec in zip(peak_Ts, peak_spectrum) if t < 100 and t > 0]
        rows = sorted(rows, key=lambda x: x[1], reverse=True)
        periods_table.append([aliases[col].split(' [')[0]] + rows[0] + rows[1])

    # rearrange columns so that dependent variable (nitrate difference) comes first
    # df = df[[df.columns[1]] + [df.columns[0]] + df.columns.tolist()[2:]]

    ##############################################################

    # multiple linear regression
    # fit nitrate difference as function of age difference, temp, pH
    x = df[['d_age', 'temperature_avg', 'pH_avg']]
    x = x.rename(columns=aliases)
    y = df['NO3+NO2_delta_frac']
    y = y.rename(columns=aliases)
    # add constant offset to fit
    x = sm.add_constant(x)
    est = sm.OLS(y, x).fit()
    # print(est.summary())

    ##############################################################
    # make plots

    def corrfunc(x, y, ax=None, **kws):
        """Plot the correlation coefficient in the top left hand corner of a plot."""
        r, _ = pearsonr(x, y)
        ax = ax or plt.gca()
        ax.annotate(f'ρ = {r:.2f}', xy=(.1, .9), xycoords=ax.transAxes)

    # make pairplot with overlaid regression lines, kde contours, and correlation coefficients
    df = df.rename(columns=aliases)
    df.index = df.index.to_pydatetime()
    g = sb.pairplot(df, corner=True, plot_kws={'s': 1, 'c': df.index, 'cmap': 'viridis'})
    # make colorbar (probably not the most elegant way to do this)
    cax = g.fig.add_axes([0.27, 0.9, 0.5, 0.025])
    cb = g.fig.colorbar(g.fig.axes[1].collections[0], cax=cax, cmap='viridis', orientation='horizontal')
    cb.ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    cb.ax.set_xticklabels([dt.datetime.fromtimestamp(t // 1000000000).strftime('%Y-%m-%d') for t in cb.get_ticks()])
    g.fig.suptitle(f'{name1[1]} $\\rightarrow$ {name2[1]}')
    g.map_lower(sb.regplot, scatter=False, line_kws={'color': 'red'})
    g.map_lower(sb.kdeplot, levels=5, color=".2")
    g.map_lower(corrfunc)

    # reformat table 0 from statsmodel summary
    d = est.summary().tables[0].data
    # rename dependent variable from nondescript 'y'
    d[0][1] = aliases['NO3+NO2_delta']
    # remove date and time from summary table
    for i in range(3, 7):
        d[i][:2] = d[i+2][:2]
    d = d[:7]

    # plot model fitness table
    ax = g.fig.add_subplot(2, 2, 2)
    ax.axis('off')
    t = ax.table(d, loc='center', colLoc='center', cellLoc='center')
    t.auto_set_font_size(False)
    t.set_fontsize(8)
    t.scale(1.3, 1.7)  # scale width, height
    t.auto_set_column_width(False)

    # plot regression coefficients table
    ax = g.fig.add_subplot(3, 3, 6)
    ax.axis('off')
    t = ax.table(est.summary().tables[1], loc='upper center', colLoc='center', cellLoc='center')
    t.auto_set_font_size(False)
    t.set_fontsize(8)
    t.scale(1.8, 1.7)  # scale width, height
    t.auto_set_column_width(False)

    # plot dominant periods table
    ax = g.fig.add_subplot(4, 4, 12)
    ax.axis('off')
    t = ax.table(periods_table, bbox=[0.2, 0.3, 1.3, 0.8], colLoc='center', cellLoc='center')
    t.auto_set_font_size(False)
    t.set_fontsize(8)
    t.scale(1.4, 1.9)  # scale width, height
    t.auto_set_column_width(False)

    # plt.show()

    # save plot
    img_fn = os.path.join(fig_dir, f'{name1[1]}_{name2[1]}.png')
    fig = g.fig
    fig.set_size_inches(10, 8)
    fig.savefig(img_fn, dpi=200, bbox_inches='tight')
    print(img_fn)
