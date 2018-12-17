# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 16:17:17 2018

@author: rustyh
"""

import sys, os

if sys.platform=='win32':
    dimr="C:\Program Files (x86)\Deltares\Delft3D FM Suite 2018.01 HMWQ (1.4.4.39490)\plugins\DeltaShell.Dimr\kernels\x64"
    dfm_bin_dir=os.path.join(dimr,'dflowfm','bin')
    share_bin_dir=os.path.join(dimr,'share','bin')
    dflowfm=os.path.join(dfm_bin_dir,'dflowfm-cli.exe')
    mpiexec=os.path.join(share_bin_dir,'mpiexec.exe')
    mapmerge=os.path.join(dfm_bin_dir,'dfmoutput.exe')
    delwaq1=os.path.join(dimr,'dwaq','bin','delwaq1.exe')
    delwaq2=os.path.join(dimr,'dwaq','bin','delwaq2.exe')
    dwaq_proc=os.path.join(dimr,'dwaq','default','proc_def')

    dimr_paths=";".join([dfm_bin_dir, share_bin_dir])

else:
    dfm_bin_dir=os.path.join(os.environ['HOME'],
                             "src/dfm/r53925-opt/bin")


def install():
    from stompy.model.delft import dflow_model
    dflow_model.DFlowModel.dfm_bin_dir=dfm_bin_dir
    dflow_model.DFlowModel.mpi_bin_dir=share_bin_dir
 
    
    
    
    
