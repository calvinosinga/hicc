#!/bin/bash

#SBATCH --share
#SBATCH --job-name=galaxy_corr
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=15000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=galaxy_corr.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/correlation/2pcf.py auto 99 100 0 nelson_mid
python3 /lustre/cosinga/HI-color/hicc/correlation/2pcf.py auto 99 100 0 nelsonrs_mid
