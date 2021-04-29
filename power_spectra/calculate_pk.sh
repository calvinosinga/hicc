#!/bin/bash
# TODO: some way for rerunning individual jobs that fail in array in slurm
# TODO: error logs for array into one file - to find individual jobs that fail

# create directories to store output
mkdir outlogs
mkdir errors

# check to make sure input is correct
if [ -z "$1" ]
then
    echo "$1 is not provided: needs snapshot 99,67,50,33"
    exit 125
else
    echo "snapshot: $1"
    
fi
if [ -z "$2" ]
then
    echo "$2 is not provided: needs box length: 100,300,50"
    exit 125
else
    echo "box length: $2"
fi
if [ -z "$3" ]
then
    echo "$3 is not provided: needs axis for line-of-sight - either 0,1,2 although typically 0"
    exit 125
else
    echo "axis: $3"
fi
if [ -z "$4" ]
then
    echo "$4 is not provided: needs number of files: 448-TNG100, 680-TNG50, 600-TNG300"
    exit 125
else
    echo "numfiles: $4"
fi

# submit hisubhalo jobs
hisubgrid=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 hisubhalo.sbatch)
hisubgrid="${hisubgrid##* }"

sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 --dependency=afterok:$hisubgrid hisubhalopk.sbatch

# submit galaxy jobs
galgrid=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 galaxy.sbatch)
galgrid="${galgrid##* }"

sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 --dependency=afterok:$hisubgrid galaxypk.sbatch

# submit hiptl jobs
hiptlgrid=$(sbatch --array=0-$4 --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 hiptl.sbatch)
hiptlgrid="${hiptlgrid##* }"

hiptlcomb=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,NUMFILES=$4 --dependency=afterok:$hiptlgrid hiptl_combine.sbatch)
hiptlcomb="${hiptlcomb##* }"

sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 --dependency=afterok:$hiptlcomb hiptlpk.sbatch

# submit ptl jobs
ptlgrid=$(sbatch --array=0-$4 --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 ptl.sbatch)
ptlgrid="${ptlgrid##* }"

ptlcomb=$(sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3,NUMFILES=$4 --dependency=afterok:$ptlgrid ptl_combine.sbatch)
ptlcomb="${ptlcomb##* }"

sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 --dependency=afterok:$ptlcomb ptlpk.sbatch

# calculate pk for paco
# TODO: analysis step, still need redshift space stuff from paco

sbatch --export=ALL,SNAP=$1,BOX=$2,AXIS=$3 pacopk.sbatch

# nelson-nelson xpk

# hisubhalo-nelson xpk

# hiptl-nelson xpk

# paco-nelson xpk

# hiptl-ptl xpk

# paco-ptl xpk

# hisubhalo-ptl xpk

# hydrotools visualization step?

# make gr-stmass plots
python /lustre/cosinga/hicc/color-stmass/gr-stmass.py $1 $2 &