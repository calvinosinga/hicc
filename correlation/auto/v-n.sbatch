#!/bin/bash

#SBATCH --share
#SBATCH --job-name=v-n_corr
#SBATCH --time=1-20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=54000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=v-n_corr.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/hicc/correlation/2pcf.py auto 99 100 0 v-n
python3 /lustre/cosinga/hicc/correlation/2pcf.py auto 99 100 0 v-nrs