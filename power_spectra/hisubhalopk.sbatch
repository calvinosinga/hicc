#!/bin/bash

#SBATCH --share
#SBATCH --job-name=hisubhalo_pk
#SBATCH --error=./errors/hisubhalo_pk.err
#SBATCH --time=18:00:00
#SBATCH --ntasks=1
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=./outlogs/hisubhalo_pk.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/hicc/pk/pk.py auto $SNAP $BOX 1 $AXIS hisubhalo
python3 /lustre/cosinga/hicc/pk/pk.py auto $SNAP $BOX 0 $AXIS hisubhalors