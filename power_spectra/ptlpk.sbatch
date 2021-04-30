#!/bin/bash

#SBATCH --share
#SBATCH --job-name=ptlpk
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=150000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=./outlogs/ptlpk.dat
#SBATCH --error=./errors/ptlpk.err

. ~/.profile
python3 /lustre/cosinga/hicc/pk/pk.py auto $SNAP $BOX 1 $AXIS ptl
python3 /lustre/cosinga/hicc/pk/pk.py auto $SNAP $BOX 0 $AXIS ptlrs