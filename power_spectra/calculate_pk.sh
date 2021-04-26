#!/bin/bash

# submit hisubhalo jobs
hisubgrid=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 hisubhalo.sbatch)
hisubgrid="${hisubgrid##* }"
sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 --dependency=afterok:$hisubgrid hisubhalopk.sbatch

# submit hiptl jobs

# split up the combine procedure into two equal chunks

hiptlgrid=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 hiptl.sbatch)

job1="${job1##* }"

job2=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 --dependency=afterok:$job1 combine.sbatch)
job2="${job2##* }"
sbatch --dependency=afterok:$job2 final.sbatch
job1=$(sbatch hiptlrs.sbatch)
job1="${job1##* }"
job2=$(sbatch --dependency=afterok:$job1 combiners.sbatch)
job2="${job2##* }"
sbatch --dependency=afterok:$job2 finalrs.sbatch