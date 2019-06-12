import numpy as np
from stompy.io.local import usgs_nwis
import stompy.model.delft.dflow_model as dfm
import matplotlib.pyplot as plt

##

six.moves.reload_module(usgs_nwis)
six.moves.reload_module(dfm)

verona_bc=dfm.NwisScalarBC(11425500,name="VeronaSSC",
                           product_id=63680,cache_dir='cache',
                           data_start=np.datetime64("2014-03-01"),
                           data_stop =np.datetime64("2014-06-01"),
                           filters=[dfm.Transform(lambda x: 2.2*x,units='mg/l'),
                                    # should signal issues down the line 
                                    dfm.FillGaps(max_gap_interp_s=20*86400,large_gap_value=np.nan)])

verona_bc.write_bokeh(filename="verona_sed.html")

##

rio_bc=dfm.NwisScalarBC(11455420,name="RioVistaSSC",
                        product_id=63680,cache_dir='cache',
                        data_start=np.datetime64("2014-03-01"),
                        data_stop =np.datetime64("2014-06-01"),
                        filters=[dfm.Transform(lambda x: 2.2*x,units='mg/l'),
                                 # should signal issues down the line 
                                 dfm.FillGaps(max_gap_interp_s=20*86400,large_gap_value=np.nan)])

# Rio Vista is weird because there are two turbidity sensors. but this appears to work okay.
rio_bc.write_bokeh(filename="riovista_sed.html")
# da=rio_bc.data()

