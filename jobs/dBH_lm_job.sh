#!/bin/bash

#SBATCH --job-name=expr
#SBATCH --output=../log/dbh_lm_job_%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH --array=1-400

ml R

LINE=$(sed -n ${SLURM_ARRAY_TASK_ID}p "dBH_lm_params.txt")
seed=$(echo $LINE | cut -d ' ' -f 1)
Xseed=$(echo $LINE | cut -d ' ' -f 2)
n=$(echo $LINE | cut -d ' ' -f 3)
p=$(echo $LINE | cut -d ' ' -f 4)
pi1=$(echo $LINE | cut -d ' ' -f 5)
mutype=$(echo $LINE | cut -d ' ' -f 6)
side=$(echo $LINE | cut -d ' ' -f 7)
nreps=$(echo $LINE | cut -d ' ' -f 8)
dBH2=$(echo $LINE | cut -d ' ' -f 9)
knockoff=$(echo $LINE | cut -d ' ' -f 10)

output=$(echo "../results/dBH-lm-n${n}-p${p}-pi1${pi1}-mutype${mutype}-side${side}-nreps${nreps}-dBH2${dBH2}-knockoff${knockoff}-seed${seed}-Xseed${Xseed}.out")

cd ../R/
Rscript dBH_lm_expr.R --n "$n" --p "$p" --pi1 "$pi1" --mutype "$mutype" --side "$side" --nreps "$nreps" --dBH2 "$dBH2" --knockoff "$knockoff" --seed "$seed" --Xseed "$Xseed" > "$output" 2>&1
