#!/bin/bash

#SBATCH --share
#SBATCH --job-name=hiptl_corr
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=15000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=hiptl_corr.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/correlation/2pcf.py auto 99 100 0 hiptl
python3 /lustre/cosinga/HI-color/hicc/correlation/2pcf.py auto 99 100 0 hiptlrs
