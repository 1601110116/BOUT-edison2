#!/bin/bash -l
#SBATCH -J data
#SBATCH -p debug
#SBATCH -o data/output
#SBATCH -e data/error
#SBATCH -N 3
#SBATCH -t 00:30:00

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#srun -n 64 ./trans_er_Nn -d data restart append   
srun -n 64 ./trans_er_Nn -d data 
