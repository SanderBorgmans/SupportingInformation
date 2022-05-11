#!/bin/bash
#
#PBS -N adaptive_WHAM
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=2

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID

# Activate custom Yaff conda environment that circumvents the Yaff logger
# For an equivalent Yaff version contact sander.borgmans@ugent.be
# or go to https://github.com/SanderBorgmans/yaff/tree/logger_overhaul
eval "$(/scratch/gent/419/vsc41947/nanoscale/anaconda3/bin/conda shell.bash hook)"
conda activate yaff_log
ml WHAM/2.0.10.2-intel-2020a-kj_mol

cd $ORIGDIR
python WHAM_calc_FES.py
