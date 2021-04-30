#!/bin/bash

#SBATCH --share
#SBATCH --job-name=nelson_pk
#SBATCH --error=./errors/nelson_pk.err
#SBATCH --time=1-20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=150000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=./outlogs/nelson_pk.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/hicc/pk/pk.py auto $SNAP $BOX 1 $AXIS nelson_low
python3 /lustre/cosinga/hicc/pk/pk.py auto $SNAP $BOX 1 $AXIS nelson_high
python3 /lustre/cosinga/hicc/pk/pk.py auto $SNAP $BOX 1 $AXIS nelson_mid
python3 /lustre/cosinga/hicc/pk/pk.py auto $SNAP $BOX 0 $AXIS nelsonrs_mid