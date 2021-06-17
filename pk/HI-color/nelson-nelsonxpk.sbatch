#!/bin/bash

#SBATCH --share
#SBATCH --mem-per-cpu=60000
#SBATCH --job-name=nelson-nelson_xpk
#SBATCH --time=2-10:00:00
#SBATCH --ntasks=1
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=nelson-nelson_xpk.dat
#SBATCH --error=nelson-nelson_xpk.err
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/hicc/pk/pk.py cross 99 100 1 0 nelson_mid nelson_mid 0
python3 /lustre/cosinga/hicc/pk/pk.py cross 99 100 0 0 nelsonrs_mid nelsonrs_mid 0
