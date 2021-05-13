#!/bin/bash

job1=$(sbatch hiptl.sbatch)
job1="${job1##* }"

job2=$(sbatch --dependency=afterok:$job1 combine.sbatch)
job2="${job2##* }"
sbatch --dependency=afterok:$job2 final.sbatch

job1=$(sbatch hiptlrs.sbatch)
job1="${job1##* }"
job2=$(sbatch --dependency=afterok:$job1 combiners.sbatch)
job2="${job2##* }"
sbatch --dependency=afterok:$job2 finalrs.sbatch