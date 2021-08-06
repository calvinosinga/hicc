#!/bin/bash

#SBATCH --share
#SBATCH --job-name=ptl-v-n_xpk
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=210000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=ptl-v-n_xpk.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py cross 99 100 1 0 ptl v-n
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py cross 99 100 0 0 ptlrs v-nrs
