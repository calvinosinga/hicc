#!/bin/bash

#SBATCH --share
#SBATCH --job-name=nelson
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=100000
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=nelson.dat

. ~/.profile
python3 /lustre/cosinga/hicc/galaxy/nelson.py 99 100 'mid'
python3 /lustre/cosinga/hicc/galaxy/nelson.py 99 100 'high'
python3 /lustre/cosinga/hicc/galaxy/nelson.py 99 100 'low'
# python3 /lustre/cosinga/hicc/galaxy/nelson.py 67 100 'mid'
# python3 /lustre/cosinga/hicc/galaxy/nelson.py 67 100 'high'
# python3 /lustre/cosinga/hicc/galaxy/nelson.py 67 100 'low'
