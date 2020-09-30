#!/bin/bash 
#PBS -P wj97
#PBS -l walltime=01:00:00 
#PBS -l ncpus=16
#PBS -l mem=64GB
#PBS -l jobfs=10GB
#PBS -l storage=scratch/wj97+gdata/v10+gdata/xu18+gdata/if87+gdata/fk4+gdata/rs0
#PBS -l software=python3
#PBS -l wd
#PBS -M amos.j.bennett@gmail.com
#PBS -m abe

module load dea
python3 dNBR_dask_tester.py