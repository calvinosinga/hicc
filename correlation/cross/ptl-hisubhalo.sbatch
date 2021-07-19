#!/bin/bash

#SBATCH --share
#SBATCH --job-name=ptl-hisubhalo_corr
#SBATCH --time=1-20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=54000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=ptl-hisubhalo_corr.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/hicc/correlation/2pcf.py cross 99 100 0 ptl hisubhalo
python3 /lustre/cosinga/hicc/correlation/2pcf.py cross 99 100 0 ptlrs hisubhalors