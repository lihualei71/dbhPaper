#!/bin/bash

#SBATCH --job-name=expr
#SBATCH --output=../log/dbh_mvt_job_%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH --array=1-400

ml R

LINE=$(sed -n ${SLURM_ARRAY_TASK_ID}p "dBH_mvt_params.txt")
seed=$(echo $LINE | cut -d ' ' -f 1)
n=$(echo $LINE | cut -d ' ' -f 2)
df=$(echo $LINE | cut -d ' ' -f 3)
pi1=$(echo $LINE | cut -d ' ' -f 4)
mutype=$(echo $LINE | cut -d ' ' -f 5)
side=$(echo $LINE | cut -d ' ' -f 6)
nreps=$(echo $LINE | cut -d ' ' -f 7)
dBH2=$(echo $LINE | cut -d ' ' -f 8)

output=$(echo "../results/dBH-mvt-n${n}-df${df}-pi1${pi1}-mutype${mutype}-side${side}-nreps${nreps}-dBH2${dBH2}-seed${seed}.out")

cd ../R/
Rscript dBH_mvt_expr.R --n "$n" --df "$df" --pi1 "$pi1" --mutype "$mutype" --side "$side" --nreps "$nreps" --dBH2 "$dBH2" --seed "$seed" > "$output" 2>&1
