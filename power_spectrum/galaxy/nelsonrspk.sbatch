#!/bin/bash

#SBATCH --share
#SBATCH --job-name=nelsonrs_pk
#SBATCH --time=1-20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=150000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=nelsonrs_pk.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/hicc/power_spectrum/pk.py auto 99 100 1 0 nelsonrs_mid
python3 /lustre/cosinga/hicc/power_spectrum/pk.py auto 67 100 1 0 nelsonrs_mid
python3 /lustre/cosinga/hicc/power_spectrum/pk.py auto 99 100 2 0 nelsonrs_mid
python3 /lustre/cosinga/hicc/power_spectrum/pk.py auto 67 100 2 0 nelsonrs_mid