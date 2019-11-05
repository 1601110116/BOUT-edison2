#PBS -q regular
#PBS -l mppwidth=512
#PBS -l walltime=2:00:00
#PBS -N my_job
#PBS -j oe
#PBS -e my_job.$PBS_JOBID.err
#PBS -o my_job.$PBS_JOBID.out
#PBS -V

cd $PBS_O_WORKDIR
qsub -W depend=afterany:$PBS_JOBID@edique02 bout_edison_regular_loop.cmd
aprun -n 1024 -j 2 ./kbm3-1 restart append

