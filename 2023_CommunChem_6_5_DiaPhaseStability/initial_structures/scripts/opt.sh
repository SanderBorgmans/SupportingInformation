#!/bin/sh
#
#PBS -N _opt_XXX
#PBS -l walltime=02:00:00
#PBS -l nodes=1:ppn=1
#PBS -m n

date

# Set up input
ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID

if [ ! -d $WORKDIR ]; then mkdir -p $WORKDIR; fi
cd $WORKDIR

cp ${ORIGDIR}/opt_uff_straindof.py $WORKDIR
cp ${ORIGDIR}/../output/COF/MAT/XXX.chk $WORKDIR
cp ${ORIGDIR}/../output/COF/pars_cluster.txt $WORKDIR
cp ${ORIGDIR}/../output/COF/pars_uff.txt $WORKDIR

# Load modules
module load LAMMPS/3Mar2020-foss-2019b-Python-3.7.4-kokkos # Load LAMMPS, also loads yaff

# Run
python opt_uff_straindof.py XXX.chk > opt.log

# Copy back results
cp ${WORKDIR}/XXX_opt.chk $ORIGDIR
cp ${WORKDIR}/opt.log $ORIGDIR

# Finalize
rm -rf $WORKDIR

date



