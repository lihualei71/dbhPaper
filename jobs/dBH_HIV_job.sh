#!/bin/bash

#SBATCH --job-name=expr
#SBATCH --output=../log/dbh_HIV.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-16:00:00

ml R
cd ../R/
Rscript dBH_lm_HIV_expr.R
