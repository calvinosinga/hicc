#!/bin/bash

#SBATCH --share
#SBATCH --job-name=hiptl_respk
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=150000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=hiptl_respk.dat

. ~/.profile
python3 /lustre/cosinga/hicc/pk/pk.py auto 99 100 1 0 hiptl_res
#python3 /lustre/cosinga/hicc/power_spectrum/pk.py auto 67 100 1 0 hiptl
