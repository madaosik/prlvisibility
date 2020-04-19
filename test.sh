#!/bin/bash

ARG_STRING=$1
mpic++ --prefix /usr/local/share/OpenMPI -o vid vid.cpp

# Get number of processors
PROC_CNT=$(mpirun --prefix /usr/local/share/OpenMPI -np 1 vid -gp $ARG_STRING)
#echo "Count of processors: $PROC_CNT"

# Run the visibility algorithm
mpirun --prefix /usr/local/share/OpenMPI -np $PROC_CNT vid $ARG_STRING