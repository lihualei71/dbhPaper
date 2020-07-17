source("expr_functions.R")
source("BHcalib.R")
nreps <- 1000
alpha <- 0.05
target <- 0.3
set.seed(1)

if (!interactive()){
    suppressPackageStartupMessages(library("argparse"))

    parser <- ArgumentParser()

    parser$add_argument("--ng", type = "integer", default = 100, help = "Number of groups")
    parser$add_argument("--nr", type = "integer", default = 50, help = "Number of replicates")
    parser$add_argument("--pi1", type = "double", default = 0.1, help = "Proportion of non-nulls")
    parser$add_argument("--mutype", type = "integer", default = 1, help = "Type of signal strength")
    parser$add_argument("--side", type = "character", default = "two", help = "Side of the tests")

    
    args <- parser$parse_args()

    ng <- args$ng
    nr <- args$nr
    pi1 <- args$pi1
    mu_type <- args$mutype
    side <- args$side
    params <- data.frame(ng, nr, pi1, mu_type, side)
} else {
    ng <- 100
    nr <- 30
    pi1 <- 0.1
    mu_type <- 2
    side <- "two"
    params <- data.frame(ng, nr, pi1, mu_type, side)
}

for (i in 1:nrow(params)){
    ng <- params[i, 1]
    nr <- params[i, 2]
    pi1 <- params[i, 3]
    mutype <- params[i, 4]
    side <- params[i, 5]
    filename <- paste0("../data/BH_calib_mcc",
                       "_ng", ng,
                       "_nr", nr,
                       "_pi1", pi1,
                       "_mutype", mutype,
                       "_side", side,
                       ".RData")
    mu1 <- BH_mcc_calib(ng, nr, pi1,
                        mutype,
                        side, nreps,
                        alpha, target)
    print(mu1)
    save(mu1, file = filename)
}
