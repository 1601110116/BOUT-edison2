#!/bin/bash
#SBATCH -N 8
#SBATCH -C knl,quad,cache
#SBATCH -p debug
#SBATCH -J jobname
#SBATCH --ntasks-per-node=64
#SBATCH -o %j.out
#SBATCH -L SCRATCH
#SBATCH -t 00:30:00

#OpenMP settings:
export OMP_NUM_THREADS=4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


#run the application:
srun -n 512 -c 4 --cpu_bind=cores ./elm_pb

# NOTE:
# $ module switch craype-haswell craype-mic-knl
# try script generator:
#   https://my.nersc.gov/script_generator.php
# You May Change Options:
#   -c: cpus per task --> default: 4
#   -N: nodes
#   -n: tasks
# if c=4, then N*64 = n

