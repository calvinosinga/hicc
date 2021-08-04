#!/bin/bash

#SBATCH --share
#SBATCH --job-name=nelson-hiptl_xpk
#SBATCH --time=1-10:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=50000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=nelson-hiptl_xpk.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py cross 99 100 1 0 nelson_mid hiptl
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py cross 99 100 0 0 nelsonrs_mid hiptlrs
