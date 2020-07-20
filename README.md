# Paper Repository

This repository contains the code to implement all examples in our paper: Conditional calibration for false discovery rate control under dependence. A new version will be released in a few days with all dBH/dSU functions implemented in the dBH package. 

## Introduction
All R scripts are included in the folder `R/`. The bash files to submit jobs to the cluster are included in the folder `jobs/` (note that this depends on your cluster and the bash file might need to be changed accordingly). The outputs and the plots are included in the folder `data/` and the folder `figs/`, respectively. 

## Installing the package
The [dbh](https://github.com/lihualei71/dbh) package needs to be installed.
```
if (!require("devtools")){
    install.packages("devtools")
}
devtools::install_github("lihualei71/dbh")
```

The following packages are required to be installed as well: [knockoff](https://cran.r-project.org/web/packages/knockoff/index.html), [tidyverse](https://www.tidyverse.org/), [ggplot2](https://ggplot2.tidyverse.org/) and [argparse](https://cran.r-project.org/web/packages/argparse/index.html).

## R scripts
The folder `R/` contains all R scripts:

<!-- - `dBH_mvgauss.R`, `dBH_mvgauss_qc.R`, `dBH_mvgauss_qc_grid.R` implement dBH and dBH$$^2$$ for multivariate z-statistics; `dBH_mvt.R`, `dBH_mvt_qc.R`, `dBH_mvt_qc_grid.R` implement dBH and dBH2 for multivariate t-statistics; `dBH_lm.R` implements dBH and dBH2 for linear models. `dBH_mvgauss_gf.R`, `dBH_mvgauss_gf_grid.R`, `dBH_mvt_gf.R` and `dBH_mvt_gf_grid.R` implement another version of dBH using a different \tau(c; X) which is not discussed in the paper. `RBH_homotopy.cpp`, `compute_knots_mvgauss.R`, `compute_knots_mvt.R` and `dBH_utils` implement intermediate functions for the aforementioned ones. All these files will be integrated into the dBH package and removed from this repo.  -->

- `dBH_mvgauss_expr.R` is an executable R script that produces the result of a single run of the simulation study for the multivariate z-statistics in Section 5 and Appendix D. It takes five inputs: `--n` for the number of hypotheses, `--pi1` for the fraction of non-nulls, `--mutype` for the type of signal strengths (1 for constant means and 2 for random means, but only `mutype = 1` is used in the paper), `--side` for the type of tests ("right" for the one-sided test H0: \mu_i <= 0$$, "left" for the one-sided test H0: \mu_i >= 0, and "two" for the two-sided test), `--nreps` for the number of replicates in each run, `--dBH2` for whether to include dBH^2 in the experiments (TRUE for all our experiments), and `--seed` for the random seed. The result will be stored in the `cluster_raw_data/` folder.

- `dBH_mvt_expr.R` is an executable R script that produces the result of a single run of the simulation study for multivariate t-statistics in Section 5 and Appendix D. The inputs include `--df` for the degree-of-freedom (n - d in Section 3.2), as well as `--n`, `--pi1`, `--mutype`, `--side`, `--nreps`, `--dBH2`, `--seed`, each with the same meaning as in `dBH_mvgauss_expr.R`. The result will be stored in the `cluster_raw_data/` folder. 

- `dBH_lm_expr.R` is an executable R script that produces the result of a single run of the simulation study for linear models in Section 5 and Appendix D. The inputs include `--n` for the sample size (instead of the number of hypotheses), `--p` for the number of covariates (which is the number of hypotheses), `--knockoff` for whether to include the fixed-X knockoffs (TRUE for all experiments in the paper), `--Xseed` for generating the fixed-design matrix, as well as `--pi1`, `--mutype`, `--side`, `--nreps`, `--dBH2`, `--seed`, each with the same meaning as in `dBH_mvgauss_expr.R`. The result will be stored in the `cluster_raw_data/` folder.

- `dBH_mcc_expr.R` is an executable R script that produces the result of a single run of the simulation study for multiple comparisons to control in Section 5 and Appendix D. The inputs include `--ng` for the number of treatment groups (which is the number of hypotheses), `--nr` for the number of replicates in each group, as well as `--pi1`, `--mutype`, `--side`, `--nreps`, `--dBH2`, `--knockoff`m `--seed`, each with the same meaning as in `dBH_lm_expr.R`. The result will be stored in the `cluster_raw_data/` folder.

- `dBH_lm_HIV_expr.R` is an executable R script that produces the results on the HIV datasets in Section 6. The result will be stored in the `data/` folder. `HIV_preprocess.R` generates the data that `dBH_lm_HIV_expr.R` requires.

- `dBH_mvgauss_expr_BHcalib.R`, `dBH_mvt_expr_BHcalib.R`, `dBH_lm_expr_BHcalib.R`, and `dBH_mcc_expr_BHcalib.R` are executable R scripts that tune the signal strength such that BH(0.05) has approximately 30% power. `BHcalib.R` implements the helpers for BH calibration.

- `expr_functions.R` implements helpers for the above experiments and `utils.R`implements all other helpers.

- `postprocess_dBH_mvgauss.R`, `postprocess_dBH_mvt.R`, `postprocess_dBH_lm.R`, `postprocess_dBH_mcc.R` and `postprocess_dBH_HIV.R` process the results obtained from the cluster and store an aggregated result in the `data/` folder.

- `plot_dBH_mvgauss_paper.R`, `plot_dBH_mvt_paper.R`, `plot_dBH_lm_paper.R`, `plot_dBH_mcc_paper.R` generate all plots in Section 5. Their counterparts without `_paper` generate all plots in Appendix D. `plot_dBH_HIV.R` generates all plots in Section 6.

- `Rcurve.R` implements functions to compute the R curves for a hypothesis and `plot_Rcurve.R` implements a function to visualize the curves. `dBH_mvgauss_Rcurve_illustration.R` generates two examples for illustration, with the figures generated by `plot_dBH_mvgauss_Rcurve.R`.

- `knockoffs.R` modifies the functions in [`knockoff` package https://cran.r-project.org/web/packages/knockoff/index.html] to generate a fixed-X knockoff matrix that is orthogonal to the intercept term. 

## Submitting jobs, Postprocessing results, and generating plots

The folder `jobs/` contains all bash scripts to submit jobs to the cluster. To start with, create the following folders.
```
mkdir log results cluster_raw_data
```
`log/` stores the system reports, `results/` stores the R stdouts, and `cluster_raw_data` stores the output/results of each job.

The next step is to calibrate the signal strengths for each experiment. To do this, run `./dBH_BHcalib.sh` or `sbatch dBH_BHcalib.sh` (to the cluster) in the `jobs/folder`. After this step, run `sbatch dBH_XXX_job.sh` for each experiment where XXX can be `mvgauss`, `mvt`, `lm`, `mcc`, and `HIV`.

Upon finishing all jobs, run `postprocess_dBH_XXX.R`. They generate `dBH_XXX_aggregate.RData` in the `data/` folder. Finally, run `plot_dBH_XXX.R` to generate all figures. 
