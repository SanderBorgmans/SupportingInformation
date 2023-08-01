#!/bin/bash
#
#PBS -N ogre_post
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=64
#PBS -A 2023_007

ORIGDIR=$PBS_O_WORKDIR

eval "$(/dodrio/scratch/users/vsc41947/data/anaconda3/bin/conda shell.bash hook)"
conda activate ogre

cd $ORIGDIR
python database_fes.py

