#!/bin/bash

#SBATCH --share
#SBATCH --job-name=nelsonrs
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=100000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=nelsonrs.dat

. ~/.profile
python3 /lustre/cosinga/hicc/galaxy/nelson.py 99 100 'mid' 0
python3 /lustre/cosinga/hicc/galaxy/nelson.py 99 100 'high' 0
python3 /lustre/cosinga/hicc/galaxy/nelson.py 99 100 'low' 0
python3 /lustre/cosinga/hicc/galaxy/nelson.py 67 100 'mid' 0
python3 /lustre/cosinga/hicc/galaxy/nelson.py 67 100 'high' 0
python3 /lustre/cosinga/hicc/galaxy/nelson.py 67 100 'low' 0
