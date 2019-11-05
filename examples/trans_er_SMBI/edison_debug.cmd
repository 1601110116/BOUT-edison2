#!/bin/bash
#SBATCH -N 16
#SBATCH -p debug
#SBATCH -J data
#SBATCH -o data/output
#SBATCH -e data/error
#SBATCH -t 00:30:00

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#run the application:
srun -n 512 ./trans_er_smbi -d data 
