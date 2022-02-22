#!/bin/bash

touch qsub_python.pbs

echo '#PBS -S /bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=00:05:00
#PBS -l mem=63gb

module load python3

cd multiorbital-superconductivity/01-main/

'$@>qsub_python.pbs

exec qsub qsub_python.pbs
