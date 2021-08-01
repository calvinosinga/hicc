#!/bin/bash

#SBATCH --share
#SBATCH --job-name=slices
#SBATCH --time=2-8:00:00
#SBATCH --ntasks=1
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=./outlogs/slices.dat
#SBATCH --error=./errors/slices.err

python3 /lustre/cosinga/HI-color/hicc/slices/slices.py ptl $SNAP $BOX
python3 /lustre/cosinga/HI-color/hicc/slices/slices.py ptlrs $SNAP $BOX
python3 /lustre/cosinga/HI-color/hicc/slices/slices.py hiptl $SNAP $BOX
python3 /lustre/cosinga/HI-color/hicc/slices/slices.py hiptlrs $SNAP $BOX
python3 /lustre/cosinga/HI-color/hicc/slices/slices.py hisubhalo $SNAP $BOX
python3 /lustre/cosinga/HI-color/hicc/slices/slices.py hisubhalors $SNAP $BOX
python3 /lustre/cosinga/HI-color/hicc/slices/slices.py v-n $SNAP $BOX
python3 /lustre/cosinga/HI-color/hicc/slices/slices.py v-nrs $SNAP $BOX
python3 /lustre/cosinga/HI-color/hicc/slices/slices.py nelson_mid $SNAP $BOX
python3 /lustre/cosinga/HI-color/hicc/slices/slices.py nelsonrs_mid $SNAP $BOX