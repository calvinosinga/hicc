#!/bin/bash

#SBATCH --share
#SBATCH --job-name=hiptl_pk
#SBATCH --error=./errors/hiptl_pk.err
#SBATCH --time=1-18:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=150000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=./outlogs/hiptl_pk.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/hicc/pk/pk.py auto $SNAP $BOX 1 $AXIS hiptl
python3 /lustre/cosinga/hicc/pk/pk.py auto $SNAP $BOX 0 $AXIS hiptlrs