#!/bin/bash

#SBATCH --share
#SBATCH --job-name=ptl_corr
#SBATCH --time=1-20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=54000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=ptl_corr.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/hicc/correlation/2pcf.py auto 99 100 0 ptl
python3 /lustre/cosinga/hicc/correlation/2pcf.py auto 99 100 0 ptlrs