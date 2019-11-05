#PBS -q debug
#PBS -l mppwidth=512
#PBS -l walltime=0:30:00
#PBS -N my_job
#PBS -j oe
#PBS -e my_job.$PBS_JOBID.err
#PBS -o my_job.$PBS_JOBID.out
#PBS -V

cd $PBS_O_WORKDIR
cp data/BOUT.inp data/BOUT.inp~
sed 's/wall_limit = [0-9.]*/wall_limit = 0.15/' < data/BOUT.inp~ > data/BOUT.inp
qsub -W depend=afterany:$PBS_JOBID@edique02 bout_edison_debug_loop.cmd
aprun -n 1024 -j 2 ./kbm3-1 restart append

