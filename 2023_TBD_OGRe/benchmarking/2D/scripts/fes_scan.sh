#!/bin/bash
#PBS -N ogre_post
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=12

# Script to scan several hyper parameters for the OGRe framework

eval "$(/scratch/gent/419/vsc41947/nanoscale/anaconda3/bin/conda shell.bash hook)"
conda activate ogre

# This will allow the worker module to work while there is still a problem with the qsub wrapper
for x in $(env | grep ^SLURM | cut -f1 -d= | sort | egrep -v '^SLURM_CLUSTERS$|^SLURM_CONF$'); do unset $x; done


ORIGDIR=$PBS_O_WORKDIR
cd $ORIGDIR

# init_spacing - init was 2.0
spacing=1.0

# init_kappa - init was 1.0
k=1.0

# kappa_growth_factor
kg=2

# Create a seperate directory if it does not exist
scan_dir=scans_${spacing}_${k}_${kg}


for n in 1 2 3; do
    for i in 0.00 0.10 0.25 0.33 0.50 0.66 0.75 0.90 1.00; do
        for j in 0.00 0.10 0.25 0.33 0.50 0.66 0.75 0.90 1.00; do
            MAX_LAYERS=$n
            CONFINEMENT_THR=$i
            OVERLAP_THR=$j
            KAPPA_GROWTH_FACTOR=${kg}
            INIT_KAPPA=$k
            INIT_SPACING=${spacing}
            
            directory=${scan_dir}/${CONFINEMENT_THR}_${OVERLAP_THR}_${MAX_LAYERS}/
            echo "Looking at ${directory}"
            
            
            if [ -d "$directory" ]; then
                cd ${directory}
                python database_fes.py
            fi
            
            # Reset
            cd $ORIGDIR
        done
    done
done
