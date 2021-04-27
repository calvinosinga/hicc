#!/bin/bash
job1=$(sbatch hiptl.sbatch)
job1="${job1##* }"
sbatch --dependency=afterok:$job1 combine.sbatch

job1=$(sbatch hiptlrs.sbatch)
job1="${job1##* }"
sbatch --dependency=afterok:$job1 combiners.sbatch