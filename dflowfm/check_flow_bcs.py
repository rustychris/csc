import stompy.model.delft.io as dio

mdu=dio.MDUFile("runs/20180802_grid97/flowfm.mdu")

t_ref,t_start,t_stop = mdu.time_range()

t0=np.datetime64("2014-04-01 00:00:00")

##
min_Q = np.loadtxt("runs/20180802_grid97/SacramentoRiver_flow_0001.tim")
sac_t =t0+min_Q[:,0]*np.timedelta64(60,'s')
sac_Q=min_Q[:,1]

min_Q = np.loadtxt("runs/20180802_grid97/AmericanRiver_flow_0001.tim")
amer_t =t0+min_Q[:,0]*np.timedelta64(60,'s')
amer_Q=min_Q[:,1]

##
plt.figure(1).clf()
plt.plot(sac_t,min_Q[:,1],label='Sac flow from tim')

##

# download sac flow
from stompy.io.local import usgs_nwis

# this will get converted from PDT to UTC internally
ds_verona=usgs_nwis.nwis_dataset(station=11425500,
                                 start_date=np.datetime64("2014-04-01"),
                                 end_date=np.datetime64("2014-06-01"),
                                 products=[60,65],cache_dir="cache")

ds_fpx=usgs_nwis.nwis_dataset(station=11447650,
                              start_date=np.datetime64("2014-04-01"),
                              end_date=np.datetime64("2014-06-01"),
                              products=[60,65],cache_dir="cache")


##
plt.figure(2).clf()

plt.plot(sac_t,sac_Q,label='Sac flow from tim')
plt.plot(amer_t,amer_Q,label='American Q from tim')
# plt.plot(sac_t,sac_Q,label='Sac+American from tim')

plt.plot(ds_verona.time-np.timedelta64(8,'h'),
         ds_verona.stream_flow_mean_daily*0.028316847,label='USGS Verona')
plt.plot(ds_fpx.time-np.timedelta64(8,'h'),
         ds_fpx.stream_flow_mean_daily*0.028316847,label='USGS Freeport')

plt.plot(ds_fpx.time-np.timedelta64(8,'h')-np.timedelta64(int(14.5*3600),'s'),
         0.9*ds_fpx.stream_flow_mean_daily*0.028316847,label='0.9(USGS Freeport)-14.5h')
ds_adj=ds_fpx.resample(time='1h').interpolate('linear')
ds_adj.time.values -= np.timedelta64(int(1.5*3600),'s')
plt.plot(ds_adj.time-np.timedelta64(8,'h'),
         0.9*ds_adj.stream_flow_mean_daily*0.028316847,label='hourly(USGS Freeport)-2h')

plt.legend()

# that shows very little tidal action at Verona.
plt.axis((735360.909580711, 735364.940986722, -291.1237478431372, 592.8520805497794))
