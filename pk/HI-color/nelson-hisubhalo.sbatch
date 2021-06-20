#!/bin/bash

#SBATCH --share
#SBATCH --job-name=nelson-hisubhalo_xpk
#SBATCH --time=2-10:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=100000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=nelson-hisubhalo_xpk.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/hicc/power_spectrum/pk.py cross 99 100 1 0 nelson_mid hisubhalo
python3 /lustre/cosinga/hicc/power_spectrum/pk.py cross 99 100 1 0 nelsonrs_mid hisubhalors