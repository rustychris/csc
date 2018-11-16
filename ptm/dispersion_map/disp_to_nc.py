# doctor that up so visit will show it
ds=xr.open_dataset('dispersion_K_v02_wsmooth.nc')

ds['time']=('time',), [np.datetime64("2014-04-01 00:00")]

for v in ['K','Ksmooth']:
    K=ds[v].values[None,:].copy()
    K[np.isnan(K)] = -999

    ds[v]=('time','face'),K

ds.to_netcdf('K_v02_with_time.nc')
ds.close()
