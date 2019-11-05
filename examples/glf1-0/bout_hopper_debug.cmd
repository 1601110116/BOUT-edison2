#PBS -q debug
#PBS -l mppwidth=256
#PBS -l mppnppn=16
#PBS -l walltime=0:15:00
#PBS -N my_job
#PBS -e my_job.$PBS_JOBID.err
#PBS -e my_job.$PBS_JOBID.out
#PBS -V

cd $PBS_O_WORKDIR
aprun -n 256 -N 16 ./elm_pb
