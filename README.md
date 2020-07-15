# Paper Repository

This repository contains the code to implement all examples in our paper: Conditional calibration for false discovery rate control under dependence. A new version will be released in a few days with all dBH/dSU functions implemented in the dBH package. 

## Introduction
All R scripts are included in the folder `R/`. The bash files to submit jobs to the cluster are included in the folder `jobs/` (note that this depends on your cluster and the bash file might need to be changed accordingly). The outputs and the plots are included in the folder `data/` and the folder `figs/`, respectively. 

## R scripts
The folder `R/` contains all R scripts:

- `dBH_mvgauss.R`, `dBH_mvgauss_qc.R`, `dBH_mvgauss_qc_grid.R` implement dBH and dBH$$^2$$ for multivariate z-statistics; `dBH_mvt.R`, `dBH_mvt_qc.R`, `dBH_mvt_qc_grid.R` implement dBH and dBH$$^2$$ for multivariate t-statistics; `dBH_lm.R` implements dBH and dBH$$^2$$ for linear models. All these files will be integrated into the dBH package and removed from this repo.

- `dBH_mvgauss_expr.R` is an executable R script that produces the result of a single run of the simulation study in the multivariate case in Section 5 and Appendix D. It takes five inputs: `--n` for the number of hypotheses, `--pi1` for the fraction of non-nulls, `--mutype` for the type of signal strengths (1 for constant means and 2 for random means, but only `mutype = 1` is used in the paper), `--side` for the type of tests ("right" for the one-sided test $$\mu_i \le 0$$, "left" for the one-sided test $$\mu_i \ge 0$$, and "two" for the two-sided test), `--nreps` for the number of replicates in each run, `--dBH2` for whether to include dBH$$^2$$ in the experiments (TRUE for all our experiments), and `--seed` for the random seed.

## Submitting jobs, Postprocessing results, and generating plots

