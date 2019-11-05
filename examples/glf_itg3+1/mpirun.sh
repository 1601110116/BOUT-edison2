#!/bin/sh
#PBS -l nodes=4:ppn=32   
#PBS -N Beer_n20

NPROCS=`wc -l < $PBS_NODEFILE`

hostname
date

cd $PBS_O_WORKDIR

cp $PBS_NODEFILE nodefile

mpirun -np $NPROCS -hostfile $PBS_NODEFILE ./itg_Beer3+1 > shotlog