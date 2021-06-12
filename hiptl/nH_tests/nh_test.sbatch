#!/bin/bash

#SBATCH --share
#SBATCH --job-name=nh_test
#SBATCH --error=nh_test.err
#SBATCH --output=nh_test.log
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=50000
#SBATCH --mail-type=ALL

. ~/.profile
python3 /lustre/cosinga/hicc/hiptl/nH_tests/nH_bins.py 125 99 100 1600 -1
