#!/bin/bash

#SBATCH --share
#SBATCH --job-name=ptl-hiptl_corr
#SBATCH --time=1-20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=ptl-hiptl_corr.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/correlation/2pcf.py cross 99 100 0 ptl hiptl
python3 /lustre/cosinga/HI-color/hicc/correlation/2pcf.py cross 99 100 0 ptlrs hiptlrs
