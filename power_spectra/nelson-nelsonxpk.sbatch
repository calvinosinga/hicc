#!/bin/bash

#SBATCH --share
#SBATCH --job-name=nelson-nelson_xpk
#SBATCH --time=2-10:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=80000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=./outlogs/nelson-nelson_xpk.dat
#SBATCH --error=./errors/nelson-nelson_xpk.err
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/hicc/pk/pk.py cross $SNAP $BOX 1 $AXIS nelson_mid nelson_mid 1
python3 /lustre/cosinga/hicc/pk/pk.py cross $SNAP $BOX 0 $AXIS nelsonrs_mid nelsonrs_mid 1