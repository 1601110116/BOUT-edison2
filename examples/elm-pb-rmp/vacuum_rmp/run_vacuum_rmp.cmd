#PBS -q regular
#PBS -l mppwidth=32
#PBS -l mppnppn=24
#PBS -l walltime=0:15:00
#PBS -N my_job
#PBS -e my_job.$PBS_JOBID.err
#PBS -e my_job.$PBS_JOBID.out
#PBS -V

cd $PBS_O_WORKDIR
aprun -n 32 -N 24 ./vacuum_rmp 




















