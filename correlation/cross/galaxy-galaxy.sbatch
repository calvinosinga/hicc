#!/bin/bash

#SBATCH --share
#SBATCH --job-name=galaxy-galaxy_corr
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=galaxy-galaxy_corr.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/correlation/2pcf.py cross 99 100 0 nelson_mid nelson_mid
python3 /lustre/cosinga/HI-color/hicc/correlation/2pcf.py cross 99 100 0 nelsonrs_mid nelsonrs_mid
