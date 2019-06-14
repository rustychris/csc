"""
Semi-automatic generation of sediment BC info.
Using the DFlowModel machinery, but not the full
model class.
"""
import stompy.model.delft.dflow_model as dfm
import numpy as np
import os
from stompy.io.local import usgs_nwis
import matplotlib.pyplot as plt

## 

model=dfm.DFlowModel.load('flowfm.mdu')
model.utc_offset=np.timedelta64(-8,'h') # PST.  
model.z_datum='NAVD88' # not really used right now.
model.projection='EPSG:26910'

cache_dir="../../cache"
##

if 0:
    verona_bc=dfm.NwisScalarBC(11425500,name="VeronaSSC",
                               model=model,
                               product_id=63680,cache_dir=cache_dir,
                               filters=[dfm.Transform(lambda x: 2.2*x,units='mg/l'),
                                        # should signal issues down the line 
                                        dfm.FillGaps(max_gap_interp_s=20*86400,large_gap_value=np.nan)])
if 1:
    # Verona is missing data around the 2014 pulse I've been looking at.
    # Try putting FPT sediment at Verona, but adjust time by 20h (based
    # on late 2014 observed pulses).  Leave scale alone - different
    # events show different offsets.
    # => This puts too much tidal signal in, and the lags are hard to
    #    interpret as they change drastically with flow conditions.
    # not quite sure about that -- and maybe I got the lag sign wrong.
    verona_bc=dfm.NwisScalarBC(11447650,name="VeronaSSC",
                               model=model,
                               product_id=63680,cache_dir=cache_dir,
                               filters=[dfm.Transform(lambda x: 2.2*x,units='mg/l'),
                                        dfm.Lag(lag=np.timedelta64(22,'h')),
                                        dfm.Lowpass(cutoff_hours=36),
                                        # should signal issues down the line 
                                        dfm.FillGaps(max_gap_interp_s=20*86400,large_gap_value=np.nan)])
    
verona_bc.write_bokeh(filename="bc_verona_sed.html")
model.write_tim(verona_bc.data(),
                os.path.join(model.run_dir,"SacramentoRiver_IM1_0001.tim"))

##

riovista_bc=dfm.NwisScalarBC(11455420,name="RioVistaSSC",
                             model=model,
                             product_id=63680,cache_dir=cache_dir,
                             filters=[dfm.Transform(lambda x: 2.2*x,units='mg/l'),
                                      # should signal issues down the line 
                                      dfm.FillGaps(max_gap_interp_s=20*86400,large_gap_value=np.nan)])

riovista_bc.write_bokeh(filename="bc_riovista_sed.html")
model.write_tim(riovista_bc.data(),
                os.path.join(model.run_dir,"SRV_IM1_0001.tim"))

