import stompy.model.delft.dflow_model as dfm
from stompy.io.local import usgs_nwis

class NwisBC(object):
    cache_dir=None
    product_id="set_in_subclass"

    def __init__(self,station,**kw):
        """
        station: int or string station id, e.g. 11455478
        """
        self.station=str(station)
        super(NwisBC,self).__init__(**kw)

class NwisStageBC(NwisBC,dfm.StageBC):
    product_id=65 # gage height
    def src_data(self):
        ds=self.fetch_for_period(self.data_start,self.data_stop)
        return ds['z']
    def write_bokeh(self,**kw):
        defaults=dict(title="Stage: %s (%s)"%(self.name,self.station))
        defaults.update(kw)
        super(NwisStageBC,self).write_bokeh(**defaults)
    def fetch_for_period(self,period_start,period_stop):
        """
        Download or load from cache, take care of any filtering, unit conversion, etc.
        Returns a dataset with a 'z' variable, and with time as UTC
        """
        ds=usgs_nwis.nwis_dataset(station=self.station,start_date=period_start,
                                  end_date=period_stop,
                                  products=[self.product_id],
                                  cache_dir=self.cache_dir)
        ds['z']=('time',), 0.3048*ds['height_gage']
        ds['z'].attrs['units']='m'
        return ds

class NwisFlowBC(NwisBC,dfm.FlowBC):
    product_id=60 # discharge
    def src_data(self):
        ds=self.fetch_for_period(self.data_start,self.data_stop)
        return ds['Q']
    def write_bokeh(self,**kw):
        defaults=dict(title="Flow: %s (%s)"%(self.name,self.station))
        defaults.update(kw)
        super(NwisFlowBC,self).write_bokeh(**defaults)
    def fetch_for_period(self,period_start,period_stop):
        """
        Download or load from cache, take care of any filtering, unit conversion, etc.
        Returns a dataset with a 'z' variable, and with time as UTC
        """
        ds=usgs_nwis.nwis_dataset(station=self.station,start_date=period_start,
                                  end_date=period_stop,
                                  products=[self.product_id],
                                  cache_dir=self.cache_dir)
        ds['Q']=('time',), 0.028316847*ds['stream_flow_mean_daily']
        ds['Q'].attrs['units']='m3 s-1'
        return ds


if 0: # dev / testing
    decker=NwisStageBC(name='decker',station=11455478,cache_dir='cache')
    decker.data_start=np.datetime64("2018-08-01")
    decker.data_stop =np.datetime64("2018-09-01")
    decker.write_bokeh()

    ges=NwisFlowBC(name='Georgiana',station=11447903,cache_dir='cache')
    ges.data_start=np.datetime64("2018-08-01")
    ges.data_stop =np.datetime64("2018-09-01")
    ges.write_bokeh()


