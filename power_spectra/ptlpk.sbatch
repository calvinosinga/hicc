#!/bin/bash

#SBATCH --share
#SBATCH --job-name=ptlpk
#SBATCH --time=2-8:00:00
#SBATCH --ntasks=1
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=./outlogs/ptlpk.dat
#SBATCH --error=./errors/ptlpk.err

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto $SNAP $BOX 1 $AXIS ptl
printf "finished ptl pk calculation \n"

python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto $SNAP $BOX 0 $AXIS ptlrs
printf "finished ptlrs pk calculation \n"

printf "now starting correlation calculations\n"

python3 /lustre/cosinga/HI-color/hicc/correlation/2pcf.py auto $SNAP $BOX $AXIS ptl
printf "finished ptl correlation\n"

python3 /lustre/cosinga/HI-color/hicc/correlation/2pcf.py auto $SNAP $BOX $AXIS ptlrs
printf "finished ptlrs correlation\n"

printf "now making slice plots\n"

python3 /lustre/cosinga/HI-color/hicc/slices/slice.py ptl $SNAP $BOX
printf "made ptl slice plot\n"

python3 /lustre/cosinga/HI-color/hicc/slices/slice.py ptlrs $SNAP $BOX
printf "made ptlrs slice plot\n"

printf "finished ptl analysis!"
