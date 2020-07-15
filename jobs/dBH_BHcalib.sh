#!/bin/bash

#SBATCH --job-name=expr
#SBATCH --output=../log/dbh_BHcalib.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00

./dBH_mvgauss_BHcalib.sh ../jobs/dBH_mvgauss_params.txt

./dBH_mvt_BHcalib.sh ../jobs/dBH_mvt_params.txt

./dBH_lm_BHcalib.sh ../jobs/dBH_lm_params.txt

./dBH_mcc_BHcalib.sh ../jobs/dBH_mcc_params.txt
