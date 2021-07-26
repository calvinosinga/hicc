#!/bin/bash

#SBATCH --share
#SBATCH --job-name=nelson
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=nelson.dat

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/galaxy/nelson.py 99 100 'mid' 800 0
python3 /lustre/cosinga/HI-color/hicc/galaxy/nelson.py 99 100 'high' 800 0
python3 /lustre/cosinga/HI-color/hicc/galaxy/nelson.py 99 100 'low' 800 0
echo 'end of nelson.sbatch file'

