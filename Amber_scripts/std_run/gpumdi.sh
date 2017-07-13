#!/bin/bash
# MOAB/Torque submission script for SciNet GPC (OpenMP)
#
#PBS -l nodes=2:ppn=12:gpus=2,walltime=12:00:00
#PBS -q gravity
#PBS -N ugrMD_LF-3

N=${1}

module purge
module load intel/15.0.2 intelmpi/5.0.3.048 gcc/4.8.1  cuda/6.0 amber/14.0

#cd /scratch/p/pmkim/ccorbi/MD_fyn/1-fyn_wt/MD_1

echo PBS_O_WORKDIR $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
#mpirun -np 24 pmemd.cuda  -O -i md.1.inp -o md.1.out -p protein.prm -c  equil.4.rst -x md.1.mdcrd  -r md.1.rst -ref equil.4.rst
mpirun -np 24  pmemd.cuda.MPI  -O -i md.2.inp -o md.$N.out -p protein.prm -c  md.$(($N-1)).rst -x md.$N.mdcrd  -r md.$N.rst

