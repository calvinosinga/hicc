#!/bin/bash

#SBATCH --share
#SBATCH --job-name=hisubhalo_pk
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=150000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=hisubhalo_pk.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/hicc/power_spectrum/pk.py auto 99 100 1 0 hisubhalo
python3 /lustre/cosinga/hicc/power_spectrum/pk.py auto 67 100 1 0 hisubhalo
