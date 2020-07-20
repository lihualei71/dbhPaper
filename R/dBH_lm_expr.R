#!/usr/bin/env Rscript

source("utils.R")
source("BHcalib.R")
source("expr_functions.R")

if (!interactive()){
    suppressPackageStartupMessages(library("argparse"))

    parser <- ArgumentParser()

    parser$add_argument("--n", type = "integer", default = 100, help = "Sample size")
    parser$add_argument("--p", type = "integer", default = 100, help = "Number of hypotheses")
    parser$add_argument("--pi1", type = "double", default = 0.2, help = "Proportion of non-nulls")
    parser$add_argument("--mutype", type = "integer", default = 1, help = "Type of signal strength")
    parser$add_argument("--side", type = "character", default = "two", help = "Side of the tests")
    parser$add_argument("--nreps", type = "integer", default = 20, help = "Number of repeats")
    parser$add_argument("--dBH2", type = "character", default = "TRUE", help = "Whether to include dBH2")
    parser$add_argument("--knockoff", type = "character", default = "TRUE", help = "Whether to include knockoff")
    parser$add_argument("--seed", type = "integer", default = 1, help = "Random seed")
    parser$add_argument("--Xseed", type = "integer", default = 2020, help = "Random seed to generate X")

    args <- parser$parse_args()

    n <- args$n
    p <- args$p
    pi1 <- args$pi1
    mu_type <- args$mutype    
    side <- args$side
    nreps <- args$nreps
    skip_dBH2 <- !as.logical(args$dBH2)
    skip_knockoff <- !as.logical(args$knockoff)
    seed <- args$seed
    Xseed <- args$Xseed
} else {
    n <- 150
    p <- 100
    pi1 <- 0.1
    mu_type <- 1
    nreps <- 5
    side <- "right"
    skip_dBH2 <- FALSE
    skip_knockoff <- FALSE
    seed <- 0
    Xseed <- 2020
}

file_root <- paste0("../cluster_raw_data/dBH_lm", 
                    "_n", n,
                    "_p", p,
                    "_pi1", pi1,
                    "_mutype", mu_type,
                    "_side", side,
                    "_nreps", nreps,
                    "_dBH2", !skip_dBH2,
                    "_knockoff", !skip_knockoff,
                    "_seed", seed,
		    "_Xseed", Xseed)
gamma <- c(0.9, NA)
geom_fac <- c(NA, 2)
alphas <- c(0.05, 0.2)
tautype <- "QC"

mu1_file <- paste0("../data/BH_calib_lm", 
    	           "_n", n,
                   "_p", p,
                   "_pi1", pi1,
                   "_mutype", mu_type,
                   "_side", side,
		   "_Xseed", Xseed,
                   ".RData")
if (!file.exists(mu1_file)){
    command <- paste0("Rscript dBH_lm_expr_BHcalib.R --n \"", n, "\" --p \"", p, "\" --pi1 \"", pi1, "\" --mutype \"", mu_type, "\" --side \"", side, "\" --Xseed \"", Xseed, "\"")
    system(command)
}
load(mu1_file)

#### Expr 1
print("Expr 1\n")
rho <- 0.8
Sigma_type <- "AR"
mu_posit_type <- "fix"

set.seed(Xseed)
SigmaX <- genSigma(p, rho, Sigma_type)
X <- mvtnorm::rmvnorm(n, sigma = SigmaX)
X <- scale(X, center = TRUE)
X <- scale(X, scale = sqrt(colMeans(X^2)))

set.seed(seed * 2020 + 1)
res <- dBH_lm_expr(X, mu1$ARplus_fix, pi1,
                   mu_posit_type, mu_type,
                   side,
                   alphas, nreps,
                   gamma = gamma,
                   geom_fac = geom_fac,
                   skip_knockoff = skip_knockoff,
                   skip_dBH2 = skip_dBH2)
print(postprocess(res))
filename <- paste0(file_root, "_", Sigma_type, "(", rho, ")_", mu_posit_type, ".RData")
save(res, file = filename)

#### Expr 2
print("Expr 2\n")
rho <- 0
Sigma_type <- "iid"
mu_posit_type <- "fix"

set.seed(Xseed)
SigmaX <- genSigma(p, rho, Sigma_type)
X <- mvtnorm::rmvnorm(n, sigma = SigmaX)
X <- scale(X, center = TRUE)
X <- scale(X, scale = sqrt(colMeans(X^2)))

set.seed(seed * 2020 + 2)
res <- dBH_lm_expr(X, mu1$iid_fix, pi1,
                   mu_posit_type, mu_type,
                   side,
                   alphas, nreps,
                   gamma = gamma,
                   geom_fac = geom_fac,
                   skip_knockoff = skip_knockoff,
                   skip_dBH2 = skip_dBH2)
print(postprocess(res))
filename <- paste0(file_root, "_", Sigma_type, "(", rho, ")_", mu_posit_type, ".RData")
save(res, file = filename)

#### Expr 3
print("Expr 3\n")
rho <- 0.5
Sigma_type <- "block"
mu_posit_type <- "fix"

set.seed(Xseed)
SigmaX <- genSigma(p, rho, Sigma_type)
X <- mvtnorm::rmvnorm(n, sigma = SigmaX)
X <- scale(X, center = TRUE)
X <- scale(X, scale = sqrt(colMeans(X^2)))

set.seed(seed * 2020 + 3)
res <- dBH_lm_expr(X, mu1$block_fix, pi1,
                   mu_posit_type, mu_type,
                   side,
                   alphas, nreps,
                   gamma = gamma,
                   geom_fac = geom_fac,
                   skip_knockoff = skip_knockoff,
                   skip_dBH2 = skip_dBH2)
print(postprocess(res))
filename <- paste0(file_root, "_", Sigma_type, "(", rho, ")_", mu_posit_type, ".RData")
save(res, file = filename)
