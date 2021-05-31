#!/bin/bash

#SBATCH --share
#SBATCH --array=0-447
#SBATCH --job-name=hiptlres
#SBATCH --error=../logs/hiptlres%a.err
#SBATCH --output=../logs/hiptlres%a.log
#SBATCH --time=1-1:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=80000
#SBATCH --mail-type=ALL

. ~/.profile
python3 /lustre/cosinga/hicc/hiptl/res_tests/hiptlres.py $SLURM_ARRAY_TASK_ID 99 100
