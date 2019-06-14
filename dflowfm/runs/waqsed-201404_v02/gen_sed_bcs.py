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

cache_dir="../../cache"
##

verona_bc=dfm.NwisScalarBC(11425500,name="VeronaSSC",
                           model=model,
                           product_id=63680,cache_dir=cache_dir,
                           filters=[dfm.Transform(lambda x: 2.2*x,units='mg/l'),
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

