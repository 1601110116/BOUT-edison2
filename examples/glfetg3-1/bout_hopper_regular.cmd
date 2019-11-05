#PBS -q regular
#PBS -l mppwidth=512
#PBS -l mppnppn=24
#PBS -l walltime=0:30:00
#PBS -N my_job
#PBS -e my_job.$PBS_JOBID.err
#PBS -e my_job.$PBS_JOBID.out
#PBS -V

cd $PBS_O_WORKDIR
aprun -n 512 -N 24 ./glf_etg
