#!/bin/bash

#SBATCH --share
#SBATCH --job-name=nelson-ptl_xpk
#SBATCH --time=1-10:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=nelson-ptl_xpk.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py cross 99 100 1 0 nelson_mid ptl
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py cross 99 100 0 0 nelsonrs_mid ptlrs
