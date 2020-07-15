#!/usr/bin/env Rscript

source("dBH_mvt.R")
source("utils.R")
source("expr_functions.R")

if (!interactive()){
    suppressPackageStartupMessages(library("argparse"))

    parser <- ArgumentParser()

    parser$add_argument("--n", type = "integer", default = 100, help = "Number of hypotheses")
    parser$add_argument("--df", type = "integer", default = 50, help = "Degree of freedom")
    parser$add_argument("--pi1", type = "double", default = 0.1, help = "Proportion of non-nulls")
    parser$add_argument("--mutype", type = "integer", default = 1, help = "Type of signal strength")
    parser$add_argument("--side", type = "character", default = "two", help = "Side of the tests")
    parser$add_argument("--nreps", type = "integer", default = 20, help = "Number of repeats")
    parser$add_argument("--dBH2", type = "character", default = "TRUE", help = "Whether to include dBH2")
    parser$add_argument("--seed", type = "integer", default = 1, help = "Random seed")

    
    args <- parser$parse_args()

    n <- args$n
    df <- args$df
    pi1 <- args$pi1
    mu_type <- args$mutype
    side <- args$side
    nreps <- args$nreps
    skip_dBH2 <- !as.logical(args$dBH2)
    seed <- args$seed
} else {
    n <- 100
    df <- 5
    pi1 <- 0.1
    mu_type <- 1
    nreps <- 5
    side <- "two"
    skip_dBH2 <- FALSE
    seed <- 601
}

set.seed(seed)
file_root <- paste0("../cluster_raw_data/dBH_mvt",
                    "_n", n,
                    "_df", df,
                    "_pi1", pi1,
                    "_mutype", mu_type,
                    "_side", side,
                    "_nreps", nreps,
                    "_dBH2", !skip_dBH2,
                    "_seed", seed)
if (side == "two"){
    ## alpha_fac <- c(0.8, NA)
    alpha_fac <- 1
} else {
    alpha_fac <- c(1, 0.8, NA)
}
gamma <- NA
alphas <- 0.05
tautype <- "QC"

mu1 <- BH_mvt_calib(n, df, pi1,
                    "fix", mu_type,
                    0.5, "equi",
                    side, nreps,
                    0.05, 0.3)

#### Expr 1
print("Expr 1\n")
rho <- 0
Sigma_type <- "iid" 
mu_posit_type <- "fix"

res <- dBH_mvt_expr(n, df, mu1, pi1,
                    mu_posit_type, mu_type,
                    rho, Sigma_type,
                    side,
                    alphas, nreps,
                    alpha_fac = alpha_fac,
                    gamma = gamma,
                    tautype = tautype,
                    skip_dBH2 = skip_dBH2,
                    if_thresh_exp = FALSE)
print(postprocess(res))
