# Reformat barker data as a time series with real time stamps
# this relied on loading csc_dfm.py enough to read in barker
# relative to that starting time stamp, then write it back out
# as CSV.

import stompy.model.delft.io as dio

#barker=model.read_tim('forcing-data/Barker_Pumping_Plant.tim',columns=['Q','s','T'],
#                      time_unit='S')

barker2=dio.read_dfm_tim('forcing-data/Barker_Pumping_Plant.tim',columns=['Q','s','T'],
                         time_unit='S',ref_time=np.datetime64("2014-04-01"))

##
import matplotlib.pyplot as plt
plt.figure(10).clf()
plt.plot(barker2.time,barker2.Q,label='reparsed')
plt.plot(barker.time,barker.Q,label='from model',lw=0.8)
plt.plot(barker3.time,barker3.Q,label='verify',lw=0.3,color='m')
plt.legend()

##

barker.to_dataframe().to_csv('forcing-data/Barker_Pumping_Plant.csv')

##
import pandas as pd

barker3=xr.Dataset.from_dataframe(pd.read_csv('forcing-data/Barker_Pumping_Plant.csv',parse_dates=['time']))


