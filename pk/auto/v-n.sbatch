#!/bin/bash

#SBATCH --share
#SBATCH --job-name=v-n_pk
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=25000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=v-n_pk.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto 99 100 1 0 v-n
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto 99 100 0 0 v-nrs
