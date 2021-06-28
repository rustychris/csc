#!/bin/bash -l
#SBATCH --job-name csc
#SBATCH -o slurm_out-%j.output
#SBATCH -e slurm_out-%j.output
#SBATCH --partition high2
#SBATCH --mem-per-cpu 2G
#SBATCH --time 1-00:00:00

# there is no -n option given in this script!
# specify -n <num cores> on the command line to sbatch.
DFMROOT=/home/rustyh/src/dfm/delft3dfm_2021.03/lnx64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DFMROOT/lib

echo 'working directory:'
echo $PWD

srun --mpi=pmi2 $DFMROOT/bin/dflowfm --autostartstop -t 1 flowfm.mdu --processlibrary $DFMROOT/share/delft3d/proc_def.def

