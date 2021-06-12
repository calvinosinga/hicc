#!/bin/bash

#SBATCH --share
#SBATCH --job-name=nelson
#SBATCH --time=15:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=100000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=./outlogs/nelson.dat
#SBATCH --error=./errors/nelson.dat

. ~/.profile
python3 /lustre/cosinga/hicc/galaxy/nelson.py $SNAP $BOX 'mid' $RES
python3 /lustre/cosinga/hicc/galaxy/nelson.py $SNAP $BOX 'high'
python3 /lustre/cosinga/hicc/galaxy/nelson.py $SNAP $BOX 'low'
# really only want the main one in redshift space
python3 /lustre/cosinga/hicc/galaxy/nelson.py $SNAP $BOX 'mid' $AXIS