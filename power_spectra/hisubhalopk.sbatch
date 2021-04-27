#!/bin/bash

#SBATCH --share
#SBATCH --job-name=hisubhalo_pk
#SBATCH --error=./errors/hisubhalo_pk.err
#SBATCH --time=2-8:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=150000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --output=./outlogs/hisubhalo_pk.dat
#SBATCH --mail-user=cosinga@umd.edu

. ~/.profile
python3 /lustre/cosinga/hicc/pk/pk.py auto $SNAP $BOX 1 $AXIS hisubhalo
python3 /lustre/cosinga/hicc/pk/pk.py auto $SNAP $BOX 1 $AXIS hisubhalors
python3 /lustre/cosinga/hicc/pk/pk.py auto $SNAP $BOX 2 $AXIS hisubhalors