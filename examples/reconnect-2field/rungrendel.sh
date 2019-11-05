#!/bin/bash

MPIEXEC=mpirun
NP=4

EXEC=2field2 
DATA=phi5+2_etab0x-4
LOG=run_${DATA}.log

echo Begin: $DATA

$MPIEXEC -np $NP ./$EXEC -d $DATA | tee $LOG
#$MPIEXEC -np $NP ./$EXEC -d $DATA restart append | tee $LOG

echo End: $DATA
