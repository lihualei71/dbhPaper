#!/bin/bash

sbatch -p candes,hns,stat dBH_mvt_job.sh

sbatch -p candes,hns,stat dBH_mvgauss_job.sh

sbatch -p candes,hns,stat dBH_mcc_job.sh

sbatch -p candes,hns,stat dBH_lm_job.sh

sbatch -p candes,hns,stat dBH_HIV_job.sh
