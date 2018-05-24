import six
import os
import shutil

from stompy.model.delft import waq_scenario as waq
from stompy.model.delft import dfm_grid

##

six.moves.reload_module(waq)

## 

hydro=waq.HydroFiles("../../dflowfm/csc_95/DFM_DELWAQ_FlowFM/FlowFM.hyd")

Sub=waq.Substance

class Scen(waq.Scenario):
    # csc_sed_9500: map output interval too small
    # csc_sed_9501: output every 30 minutes
    base_path="csc_sed_9501"
    map_formats=['binary']

    # With integration_option="15.60"
    # 3000 was unstable.
    # 1000 was unstable.
    #  500 stable at least to 5%, running ~ 600x realtime
    # Also tried 1000 with integration_option="21.70" and "16.70", no better.
    # 21.70 appears to be about 15% faster than 15.60, a little over 700x.
    #  that was outputting a lot.  
    integration_option="21.70" 
    time_step=500 
    map_time_step=3000

    def init_parameters(self):
        params=super(Scen,self).init_parameters()


        params['TauShields']= 5.00000e-001 #    1 uitwisseling met bufferlaag in Scheldemodel circa 2x sneller, verblijftijd 4x korter
        params['GRAIN50']= 3.00000e-004 # 
        params['GRAV']= 9.80000e+000 # 
        params['KinViscos']= 1.00000e-006 # 
        params['RHOSAND']= 2.60000e+006 # 
        params['RhoWater']= 1.02000e+003 # 
        params['PORS2']= 4.00000e-001 # 
        params['ThickS2']= 5.00000e-002 #       0.1
        params['MinDepth']= 5.00000e-002 #      0.01
        params['MaxResPup']= 1.00000e+020 # 
        params['FactResPup']= 1.75000e-007 #    3.5e-07
        params['V0SedIM1']= 1.10000e+001 #      129.6 m/d
        params['CrSS']= 2.50000e+001 # 
        params['nIM1']= 1.00000e-001 # 
        params['TaucSIM1']= 1.00000e+003 # 
        params['TcSED']= 1.01000e+002 #  
        params['FrIM1SedS2']= 1.00000e-001 #    0.05
        params['FrTIMS2Max']= 1.0 # 
        params['SWResIM1']= 1.0 # 
        params['SWResusp']= 1.0 # 
        params['ZResIM1']= 8.64000e+005 #       8.64e3  dikte flufflaag in Scheldemodel circa 10x kleiner, mag wel 5x tot 10x groter zijn
        params['VResIM1']= 1.00000e+001 #       4e-2
        params['TaucRS1IM1']= 2.00000e-001 # 
        params['TaucRS2IM1']= 1.00000e+002 # 
        params['V0SedIM2']= 2.20000e+001 # 
        params['nIM2']= 1.00000e-001 # 
        params['TaucSIM2']= 1.00000e+003 # 
        params['FrIM2SedS2']= 1.00000e-001 #         0.05
        params['SWResIM2']= 1.0 # 
        params['ZResIM2']= 4.23000e+005 #            8.64e3
        params['VResIM2']= 2.00000e-001 #            4e-2
        params['TaucRS1IM2']= 1.50000e-001 # 
        params['TaucRS2IM2']= 1.00000e+002 # 

        params['Manncoef']= 2.40000e-002
        params['SwChezy']= 2.00000e+000 
        params['SwTau']= 2.00000e+000  
        params['Rough']= 5.00000e-002  
        params['CLOSE_ERR']= 1 
        params['ScaleVDisp']= 1.0 
        params['MaxIter']= 100  
        params['Tolerance']= 1.0e-7
        params['Iteration Report']= 0.00000e+000 
        params['ONLY_ACTIVE']= 1.0
        
        params['ACTIVE_Res_Pickup']= 1.0 
        params['ACTIVE_Sed_IM1']= 1.0 
        params['ACTIVE_Res_IM1']= 1.0 
        params['ACTIVE_Sed_IM2']= 1.0 
        params['ACTIVE_Res_IM2']= 1.0 
        params['ACTIVE_CalVS_IM1']= 1.0
        params['ACTIVE_CalVS_IM2']= 1.0
        params['ACTIVE_Tau']= 1.0 
        params['ACTIVE_DynDepth']= 1.0 
        params['ACTIVE_Res_DM']= 1.0 
        params['ACTIVE_S1_Comp']= 1.0
        params['ACTIVE_S2_Comp']= 1.0
        params['ACTIVE_Compos']= 1.0 
        params['ACTIVE_Veloc']= 1.0 
        params['ACTIVE_Chezy']= 1.0  
        params['ACTIVE_TotDepth']= 1.0

        if self.hydro.n_exch_z>0:
            # if we don't have vertical diffusion specified, this will
            # raise an error, but that's good.  3D runs should have it.
            params['ACTIVE_VertDisp']= 1.0
            
        params['ACTIVE_Tau']= 1.0

        # also need Tau, maybe surf, bottomdept, vertdisper
        # tau, vertdiper, surf should all be populated from HydroFiles automatically.
        return params

    def init_substances(self):
        subs=super(Scen,self).init_substances()

        subs['IM1']=Sub(10)
        subs['IM2']=Sub(10)
        subs['IM1S1']=Sub(100,active=False)
        subs['IM1S2']=Sub(100,active=False)
        subs['IM2S1']=Sub(100,active=False)
        subs['IM2S2']=Sub(100,active=False)
        subs['Zsand']=Sub(100,active=False)
        
        return subs
    
    #def init_bcs(self):
    #def init_loads(self):

    
scen=Scen(hydro)
scen.map_output += ('tau','ss','Depth',
                    'fSedIM1' ,
                    'fResS1IM1',
                    'fResS1IM2' ,
                    'FrIM1S1',
                    'fResS2DM',
                    'FrIM1S2',
                    'fSedIM2',
                    'FrIM2S1',
                    'FrIM2S2',
                    'DMS1',
                    'DMS2',
                    'ZResDM','Velocity')

os.environ['DELFT_BIN']="/Users/rusty/src/dfm/r52184-opt/bin"
os.environ['DELFT_SRC']="/Users/rusty/src/dfm/r52184-opt/build/dwaq/delft3d/src"

os.path.exists(scen.base_path) and shutil.rmtree(scen.base_path) # for dev

scen.cmd_default()
scen.cmd_delwaq1()
scen.cmd_delwaq2()

