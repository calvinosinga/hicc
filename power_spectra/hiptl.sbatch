#!/bin/bash

#SBATCH --share
#SBATCH --job-name=hiptl
#SBATCH --error=./errors/hiptl%a.err
#SBATCH --output=./outlogs/hiptl%a.log
#SBATCH --time=2-1:00:00
#SBATCH --ntasks=1
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --account=astronomy-hi
#SBATCH --mem-per-cpu=80000
#SBATCH --mail-type=ALL

. ~/.profile
python3 /lustre/cosinga/hicc/hiptl/hiptl.py $SLURM_ARRAY_TASK_ID $SNAP $BOX
python3 /lustre/cosinga/hicc/hiptl/hiptl.py $SLURM_ARRAY_TASK_ID $SNAP $BOX $AXIS