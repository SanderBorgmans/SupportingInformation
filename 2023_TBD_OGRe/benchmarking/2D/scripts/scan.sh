#!/bin/bash
#PBS -N ogre_scan_ackley
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=4

# Script to scan several hyper parameters for the OGRe framework

ml worker 

eval "$(/scratch/gent/419/vsc41947/nanoscale/anaconda3/bin/conda shell.bash hook)"
conda activate ogre

# This will allow the worker module to work while there is still a problem with the qsub wrapper
for x in $(env | grep ^SLURM | cut -f1 -d= | sort | egrep -v '^SLURM_CLUSTERS$|^SLURM_CONF$'); do unset $x; done


ORIGDIR=$PBS_O_WORKDIR
cd $ORIGDIR


# iterate over scan directories
for i in  "2.0 1.0 2"; do
    set -- $i # convert the "tuple" into the param args $1 $2...

    # init_spacing - init was 2.0
    spacing=$1

    # init_kappa - init was 1.0
    k=$2

    # kappa_growth_factor
    kg=$3
    
    
    # Load the correct directory
    scan_dir=scans_${spacing}_${k}_${kg}

    if [ ! -d ${scan_dir} ]; then
        mkdir ${scan_dir}
    fi


    if [ -d ${scan_dir} ]; then
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

                    
                    if [ ! -d "$directory" ]; then
                        mkdir ${directory}
                        cp scripts/*.* ${directory}
                        
                        # initialize
                        cd ${directory}
                        python database_init.py $MAX_LAYERS $CONFINEMENT_THR $OVERLAP_THR $INIT_SPACING $INIT_KAPPA $KAPPA_GROWTH_FACTOR
                        
                        
                        # while the run file contains more entries than the header
                        while [[ $(find run.txt -type f -size +25c 2>/dev/null) ]]; do
                            python database_simulate_load_log.py
                            
                            if [ ! -f "grid_restart.txt" ]; then
                                python database_post.py
                            else
                                # if additional simulations needs to be performed it could be that some are missing 
                                # (due to e.g. insufficient simulation time)
                                while [ -f "grid_restart.txt" ]; do
                                    echo "Performing additional simulations! --------------------"
                                    python adapt_cores.py
                                    wsub -batch run.pbs -data grid_restart.txt -epilog epilog.sh --master
                                    
                                    # while loop to wait for these simulations to finish
                                    while [[ ! -f done ]]; do
                                        sleep 10
                                    done
                                    
                                    # Remove the 'done' file from the epilog script to reset for the next iteration
                                    rm done
                                    rm grid_restart.txt
                                    
                                    python database_post.py
                                done
                            fi
                        done
                    fi
                    
                    # Reset
                    cd $ORIGDIR
                done
            done
        done
    fi
done
