#!/usr/bin/env Rscript

source("dBH_lm.R")
source("utils.R")
source("BHcalib.R")
source("expr_functions.R")
library("knockoff")

if (!interactive()){
    suppressPackageStartupMessages(library("argparse"))

    parser <- ArgumentParser()

    parser$add_argument("--ng", type = "integer", default = 100, help = "Number of groups")
    parser$add_argument("--nr", type = "integer", default = 3, help = "Number of replicates")
    parser$add_argument("--pi1", type = "double", default = 0.2, help = "Proportion of non-nulls")
    parser$add_argument("--mutype", type = "integer", default = 1, help = "Type of signal strength")
    parser$add_argument("--side", type = "character", default = "two", help = "Side of the tests")
    parser$add_argument("--nreps", type = "integer", default = 20, help = "Number of repeats")
    parser$add_argument("--dBH2", type = "character", default = "TRUE", help = "Whether to include dBH2")
    parser$add_argument("--knockoff", type = "character", default = "TRUE", help = "Whether to include knockoff")
    parser$add_argument("--seed", type = "integer", default = 1, help = "Random seed")

    args <- parser$parse_args()

    ng <- args$ng
    nr <- args$nr
    pi1 <- args$pi1
    mu_type <- args$mutype    
    side <- args$side
    nreps <- args$nreps
    skip_dBH2 <- !as.logical(args$dBH2)
    skip_knockoff <- !as.logical(args$knockoff)
    seed <- args$seed
} else {
    ng <- 100
    nr <- 30
    pi1 <- 0.1
    mu_type <- 1
    nreps <- 5
    side <- "two"
    skip_dBH2 <- FALSE
    skip_knockoff <- FALSE
    seed <- 0
}

set.seed(seed)
file_root <- paste0("../cluster_raw_data/dBH_mcc", 
                    "_ng", ng,
                    "_nr", nr,
                    "_pi1", pi1,
                    "_mutype", mu_type,
                    "_side", side,
                    "_nreps", nreps,
                    "_dBH2", !skip_dBH2,
                    "_knockoff", !skip_knockoff,
                    "_seed", seed)
gamma <- c(0.9, NA)
beta <- c(NA, 2)
alphas <- c(0.05, 0.2)
tautype <- "QC"

mu1_file <- paste0("../data/BH_calib_mcc",
                   "_ng", ng,
                   "_nr", nr,
                   "_pi1", pi1,
                   "_mutype", mu_type,
                   "_side", side,
                   ".RData")
if (!file.exists(mu1_file)){
    command <- paste0("Rscript dBH_mcc_expr_BHcalib.R --ng \"", ng, "\" --nr \"", nr, "\" --pi1 \"", pi1, "\" --mutype \"", mu_type, "\" --side \"", side, "\"")
    system(command)
}
load(mu1_file)

#### Expr 1
print("Expr 1\n")

res <- dBH_mcc_expr(ng, nr,
                    mu1, pi1,
                    "fix", mu_type,
                    side,
                    alphas, nreps,
                    gamma = gamma,
                    beta = beta,
                    skip_knockoff = skip_knockoff,
                    skip_dBH2 = skip_dBH2)
print(postprocess(res))
filename <- paste0(file_root, ".RData")
## save(res, file = filename)
