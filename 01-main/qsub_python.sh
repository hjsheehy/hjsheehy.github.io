#!/bin/bash

touch qsub_python.pbs

echo '#PBS -S /bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=00:05:00
#PBS -l mem=63gb

MAIN=multiorbital-superconductivity/01-main/

cd $MAIN

module load python3
'$@>qsub_python.pbs

exec qsub qsub_python.pbs
