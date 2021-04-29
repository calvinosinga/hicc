#!/bin/bash

#SBATCH --share
#SBATCH --job-name=hiptl-nelson_xpk
#SBATCH --time=2-10:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=230000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=./outlogs/hiptl-nelson_xpk.dat
#SBATCH --error=./errors/hiptl-nelson_xpk.err
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/hicc/pk/pk.py cross $SNAP $BOX 1 $AXIS nelson_mid hiptl