#!/bin/bash

#SBATCH --share
#SBATCH --array=0-447
#SBATCH --job-name=nh
#SBATCH --error=nh%a.err
#SBATCH --output=nh%a.log
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=50000
#SBATCH --mail-type=ALL

. ~/.profile
python3 /lustre/cosinga/hicc/hiptl/nH_tests/nH_bins.py $SLURM_ARRAY_TASK_ID 99 100 1600 -1
python3 /lustre/cosinga/hicc/hiptl/nH_tests/nH_bins.py $SLURM_ARRAY_TASK_ID 99 100 1600 0
