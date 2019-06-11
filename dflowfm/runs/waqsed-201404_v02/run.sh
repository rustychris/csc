#!/bin/bash

DFMROOT=$HOME/src/dfm/1.5.2/lnx64
DFMBIN=$DFMROOT/bin
export LD_LIBRARY_PATH=$DFMROOT/lib

[ -d DFM_OUTPUT_flowm ] && rm -r DFM_OUTPUT_flowfm

$DFMBIN/dflowfm --partition:ndomains=12:icgsolver=6 flowfm.mdu || exit 1 
	
/opt/anaconda3/bin/mpiexec -n 12 $DFMBIN/dflowfm --autostartstop flowfm.mdu --processlibrary $DFMROOT/share/delft3d/proc_def.def
