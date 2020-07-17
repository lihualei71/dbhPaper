source("expr_functions.R")
source("BHcalib.R")
nreps <- 1000
alpha <- 0.05
target <- 0.3
set.seed(1)

if (!interactive()){
    suppressPackageStartupMessages(library("argparse"))

    parser <- ArgumentParser()

    parser$add_argument("--n", type = "integer", default = 100, help = "Number of hypotheses")
    parser$add_argument("--df", type = "integer", default = 50, help = "Degree of freedom")
    parser$add_argument("--pi1", type = "double", default = 0.1, help = "Proportion of non-nulls")
    parser$add_argument("--mutype", type = "integer", default = 1, help = "Type of signal strength")
    parser$add_argument("--side", type = "character", default = "two", help = "Side of the tests")

    
    args <- parser$parse_args()
    n <- args$n
    df <- args$df
    pi1 <- args$pi1
    mu_type <- args$mutype
    side <- args$side
    params <- data.frame(n, df, pi1, mu_type, side)
} else {
    n <- 100
    df <- 50
    pi1 <- 0.1
    mu_type <- 1
    side <- "two"
    params <- data.frame(n, df, pi1, mu_type, side)
}

for (i in 1:nrow(params)){
    n <- params[i, 1]
    df <- params[i, 2]
    pi1 <- params[i, 3]
    mutype <- params[i, 4]
    side <- params[i, 5]
    filename <- paste0("../data/BH_calib_mvt", 
    	               "_n", n,
                       "_df", df,
                       "_pi1", pi1,
                       "_mutype", mutype,
                       "_side", side,
                       ".RData")
    mu1_ARplus_fix <- BH_mvt_calib(n, df, pi1,
                                   "fix", mutype,
                                   0.8, "AR",
                                   side, nreps,
                                   alpha, target)
    mu1_iid_fix <- BH_mvt_calib(n, df, pi1,
                                "fix", mutype,
                                0, "iid",
                                side, nreps,
                                alpha, target)
    mu1_block_fix <- BH_mvt_calib(n, df, pi1,
                                  "fix", mutype,
                                  0.5, "block",
                                  side, nreps,
                                  alpha, target)
    mu1 <- list(ARplus_fix = mu1_ARplus_fix,
                iid_fix = mu1_iid_fix,
                block_fix = mu1_block_fix)
    print(mu1)
    save(mu1, file = filename)
}
