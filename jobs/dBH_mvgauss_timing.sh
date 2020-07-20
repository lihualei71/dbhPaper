#!/bin/bash

#SBATCH --job-name=expr
#SBATCH --output=../log/dbh_mvgauss_timing_%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH --array=1-200

ml R

LINE=$(sed -n ${SLURM_ARRAY_TASK_ID}p "dBH_mvgauss_timing_params.txt")
seed=$(echo $LINE | cut -d ' ' -f 1)
nalt=$(echo $LINE | cut -d ' ' -f 2)
rho=$(echo $LINE | cut -d ' ' -f 3)
fac=$(echo $LINE | cut -d ' ' -f 4)
nreps=$(echo $LINE | cut -d ' ' -f 5)
dBH2=$(echo $LINE | cut -d ' ' -f 6)

output=$(echo "../results/dBH-mvgauss-timing-nalt${nalt}-rho${rho}-fac${fac}-nreps${nreps}-dBH2${dBH2}-seed${seed}.out")

cd ../R/
Rscript dBH_mvgauss_timing.R --nalt "$nalt" --rho "$rho" --fac "$fac" --nreps "$nreps" --dBH2 "$dBH2" --seed "$seed" > "$output" 2>&1
