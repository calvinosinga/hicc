#!/bin/bash

#SBATCH --share
#SBATCH --job-name=hiptlrspk
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=150000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=hiptlrspk.dat

. ~/.profile
python3 /lustre/cosinga/hicc/pk/pk.py auto 99 100 0 0 hiptlrs
#python3 /lustre/cosinga/hicc/power_spectrum/pk.py auto 99 100 2 0 hiptlrs
