#!/bin/bash

#
#PBS -N ogre_fes
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1

# Script to scan several hyper parameters for the OGRe framework

eval "$(/scratch/gent/419/vsc41947/nanoscale/anaconda3/bin/conda shell.bash hook)"
conda activate ogre

ORIGDIR=$PBS_O_WORKDIR
cd $ORIGDIR


# iterate over scan directories
for i in "1.0 0.1 2" "1.0 1.0 2" ; do
    set -- $i # convert the "tuple" into the param args $1 $2...

    # init_spacing - init was 2.0
    spacing=$1

    # init_kappa - init was 1.0
    k=$2

    # kappa_growth_factor
    kg=$3

    # Load the correct directory
    scan_dir=scans_${spacing}_${k}_${kg}
    
    
    if [ -d ${scan_dir} ]; then
        for n in 1 2 3 4 5; do
            for i in 0.00 0.10 0.25 0.33 0.50 0.66 0.75 0.90 1.00; do
                for j in 0.00 0.10 0.25 0.33 0.50 0.66 0.75 0.90 1.00; do
                    MAX_LAYERS=$n
                    CONFINEMENT_THR=$i
                    OVERLAP_THR=$j
                    KAPPA_GROWTH_FACTOR=${kg}
                    INIT_KAPPA=$k
                    
                    directory=${scan_dir}/${CONFINEMENT_THR}_${OVERLAP_THR}_${KAPPA_GROWTH_FACTOR}_${MAX_LAYERS}/
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
    fi
done
