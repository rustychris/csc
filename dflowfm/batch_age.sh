#!/bin/bash -l
#SBATCH --job-name csc
#SBATCH -o slurm_out-%j.output
#SBATCH -e slurm_out-%j.output
#SBATCH --partition high2
#SBATCH --mem-per-cpu 2G
#SBATCH --time 1-00:00:00
#SBATCH -n 16
#SBATCH -N 1

# srun --mpi=pmi2 $DFMROOT/bin/dflowfm --autostartstop -t 1 flowfm.mdu --processlibrary $DFMROOT/share/delft3d/proc_def.def
conda activate general
python csc_age_2014.py


