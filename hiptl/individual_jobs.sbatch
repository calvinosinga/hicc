#!/bin/bash

#SBATCH --share
#SBATCH --output=./logs/ind.log
#SBATCH --time=2-1:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=80000
#SBATCH --mail-type=ALL

. ~/.profile
python3 /lustre/cosinga/hicc/hiptl/redshift_space/hiptlrs.py 37 99 100 0
python3 /lustre/cosinga/hicc/hiptl/redshift_space/hiptlrs.py 26 99 100 0
