#!/bin/bash

#SBATCH --share
#SBATCH --job-name=nelson_pk
#SBATCH --time=1-20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=25000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=nelson_pk.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto 99 100 1 0 nelson_low
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto 99 100 1 0 nelson_high
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto 99 100 1 0 nelson_mid
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto 99 100 0 0 nelsonrs_mid
#python3 /lustre/cosinga/hicc/power_spectrum/pk.py auto 67 100 1 0 nelson_low
#python3 /lustre/cosinga/hicc/power_spectrum/pk.py auto 67 100 1 0 nelson_high
#python3 /lustre/cosinga/hicc/power_spectrum/pk.py auto 67 100 1 0 nelson_mid
