#!/bin/sh
#PBS -l nodes=1:ppn=16
#PBS relax
cd $PBS_O_WORKDIR
export I_MPI_COMPATIBILITY=4
NPROCS=$(cat $PBS_NODEFILE | wc -l)

python test.py
