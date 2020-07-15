#!/bin/bash

#SBATCH --job-name=expr
#SBATCH --output=../log/dbh_mcc_job_%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH --array=1-400

ml R

LINE=$(sed -n ${SLURM_ARRAY_TASK_ID}p "dBH_mcc_params.txt")
seed=$(echo $LINE | cut -d ' ' -f 1)
ng=$(echo $LINE | cut -d ' ' -f 2)
nr=$(echo $LINE | cut -d ' ' -f 3)
pi1=$(echo $LINE | cut -d ' ' -f 4)
mutype=$(echo $LINE | cut -d ' ' -f 5)
side=$(echo $LINE | cut -d ' ' -f 6)
nreps=$(echo $LINE | cut -d ' ' -f 7)
dBH2=$(echo $LINE | cut -d ' ' -f 8)
knockoff=$(echo $LINE | cut -d ' ' -f 9)

output=$(echo "../results/dBH-mcc-ng${ng}-nr${nr}-pi1${pi1}-mutype${mutype}-side${side}-nreps${nreps}-dBH2${dBH2}-knockoff${knockoff}-seed${seed}.out")

cd ../R/
Rscript dBH_mcc_expr.R --ng "$ng" --nr "$nr" --pi1 "$pi1" --mutype "$mutype" --side "$side" --nreps "$nreps" --dBH2 "$dBH2" --knockoff "$knockoff" --seed "$seed" > "$output" 2>&1
