#!/usr/bin/env Rscript

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
    df <- 50
    pi1 <- 0.1
    mu_type <- 1
    nreps <- 5
    side <- "right"
    skip_dBH2 <- FALSE
    seed <- 0
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
geom_fac <- c(NA, 2)
alphas <- 0.05
tautype <- "QC"

mu1_file <- paste0("../data/BH_calib_mvt", 
                   "_n", n,
                   "_df", df,
                   "_pi1", pi1,
                   "_mutype", mu_type,
                   "_side", side,
                   ".RData")
if (!file.exists(mu1_file)){
    command <- paste0("Rscript dBH_mvt_expr_BHcalib.R --n \"", n, "\" --df \"", df, "\" --pi1 \"", pi1, "\" --mutype \"", mu_type, "\" --side \"", side, "\"")
    system(command)
}
load(mu1_file)

#### Expr 1
print("Expr 1\n")
rho <- 0.8
Sigma_type <- "AR" 
mu_posit_type <- "fix"
gamma <- c(0.9, NA)

res <- dBH_mvt_expr(n, df, mu1$ARplus_fix, pi1,
                    mu_posit_type, mu_type,
                    rho, Sigma_type,
                    side,
                    alphas, nreps,
                    gamma = gamma,
                    geom_fac = geom_fac,
                    tautype = tautype,
                    skip_dBH2 = skip_dBH2)
print(postprocess(res))
filename <- paste0(file_root, "_", Sigma_type, "(", rho, ")_", mu_posit_type, ".RData")
save(res, file = filename)

#### Expr 2
print("Expr 2\n")
rho <- 0
Sigma_type <- "iid" 
mu_posit_type <- "fix"
gamma <- c(1, NA)

res <- dBH_mvt_expr(n, df, mu1$iid_fix, pi1,
                    mu_posit_type, mu_type,
                    rho, Sigma_type,
                    side,                        
                    alphas, nreps,
                    gamma = gamma,
                    geom_fac = geom_fac,
                    tautype = tautype,
                    skip_dBH2 = skip_dBH2)
print(postprocess(res))
filename <- paste0(file_root, "_", Sigma_type, "(", rho, ")_", mu_posit_type, ".RData")
save(res, file = filename)

#### Expr 3
print("Expr 3\n")
rho <- 0.5
Sigma_type <- "block" 
mu_posit_type <- "fix"
gamma <- c(0.9, NA)

res <- dBH_mvt_expr(n, df, mu1$block_fix, pi1,
                    mu_posit_type, mu_type,
                    rho, Sigma_type,
                    side,                        
                    alphas, nreps,
                    gamma = gamma,
                    geom_fac = geom_fac,
                    tautype = tautype,
                    skip_dBH2 = skip_dBH2)
print(postprocess(res))
filename <- paste0(file_root, "_", Sigma_type, "(", rho, ")_", mu_posit_type, ".RData")
save(res, file = filename)
