#!/bin/bash

#SBATCH --share
#SBATCH --job-name=v-npk
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --account=astronomy-hi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosinga@umd.edu
#SBATCH --output=./outlogs/v-npk.dat
#SBATCH --error=./errors/v-npk.err

. ~/.profile
python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto $SNAP $BOX 1 $AXIS v-n
printf "finished v-n pk calculation \n"

python3 /lustre/cosinga/HI-color/hicc/pk/pk.py auto $SNAP $BOX 0 $AXIS v-nrs
printf "finished v-nrs pk calculation \n"

printf "now starting correlation calculations\n"

python3 /lustre/cosinga/HI-color/hicc/correlation/2pcf.py auto $SNAP $BOX $AXIS v-n
printf "finished v-n correlation\n"

python3 /lustre/cosinga/HI-color/hicc/correlation/2pcf.py auto $SNAP $BOX $AXIS v-nrs
printf "finished v-nrs correlation\n"

printf "now making slice plots\n"

python3 /lustre/cosinga/HI-color/hicc/slices/slice.py v-n $SNAP $BOX
printf "made v-n slice plot\n"

python3 /lustre/cosinga/HI-color/hicc/slices/slice.py v-nrs $SNAP $BOX
printf "made v-nrs slice plot\n"

printf "finished v-n analysis!"
