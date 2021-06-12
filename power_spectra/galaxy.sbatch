#!/bin/bash

#SBATCH --share
#SBATCH --job-name=nelson
#SBATCH --time=15:00:00
#SBATCH --ntasks=1
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=./outlogs/nelson.dat
#SBATCH --error=./errors/nelson.dat

. ~/.profile
python3 /lustre/cosinga/hicc/galaxy/nelson.py $SNAP $BOX 'mid' $RES -1
python3 /lustre/cosinga/hicc/galaxy/nelson.py $SNAP $BOX 'high' $RES -1
python3 /lustre/cosinga/hicc/galaxy/nelson.py $SNAP $BOX 'low' $RES -1
# really only want the main one in redshift space
python3 /lustre/cosinga/hicc/galaxy/nelson.py $SNAP $BOX 'mid' $RES $AXIS