#!/bin/bash
#
#PBS -N US_MATERIAL_FOLD_TEMPERATURE
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=all

ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID

# Activate custom Yaff conda environment that circumvents the Yaff logger
# For an equivalent Yaff version contact sander.borgmans@ugent.be
# or go to https://github.com/SanderBorgmans/yaff/tree/logger_overhaul
eval "$(/scratch/gent/419/vsc41947/nanoscale/anaconda3/bin/conda shell.bash hook)"
conda activate yaff_log

mkdir -p $WORKDIR
cp $ORIGDIR/*.* $WORKDIR/.

cd $WORKDIR
python US_simulate.py

cp -r $WORKDIR/logs/ $ORIGDIR/.
cp -r $WORKDIR/trajs/ $ORIGDIR/.

cd $ORIGDIR
rm -rf $WORKDIR
