#!/bin/bash

#SBATCH --share
#SBATCH --array=0-447
#SBATCH --job-name=numdens
#SBATCH --error=./logs/numdens%a.err
#SBATCH --output=./logs/numdens%a.log
#SBATCH --time=2-1:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=20000
#SBATCH --mail-type=ALL

. ~/.profile
python3 /lustre/cosinga/hicc/hiptl/numdens.py $SLURM_ARRAY_TASK_ID 99 100
