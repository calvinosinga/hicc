#!/bin/bash

#SBATCH --share
#SBATCH --job-name=ptl-galaxy_corr
#SBATCH --time=1-20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=54000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=ptl-galaxy_corr.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/hicc/correlation/2pcf.py cross 99 100 0 ptl nelson_mid
python3 /lustre/cosinga/hicc/correlation/2pcf.py cross 99 100 0 ptlrs nelsonrs_mid