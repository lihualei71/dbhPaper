source("expr_functions.R")
source("BHcalib.R")
nreps <- 1000
alpha <- 0.05
target <- 0.3
set.seed(1)

if (!interactive()){
    suppressPackageStartupMessages(library("argparse"))

    parser <- ArgumentParser()

    parser$add_argument("--n", type = "integer", default = 100, help = "Sample size")
    parser$add_argument("--p", type = "integer", default = 100, help = "Number of hypotheses")
    parser$add_argument("--pi1", type = "double", default = 0.2, help = "Proportion of non-nulls")
    parser$add_argument("--mutype", type = "integer", default = 1, help = "Type of signal strength")
    parser$add_argument("--side", type = "character", default = "two", help = "Side of the tests")
    parser$add_argument("--Xseed", type = "integer", default = 2020, help = "Random seed to generate X")
    
    args <- parser$parse_args()

    n <- args$n
    p <- args$p
    pi1 <- args$pi1
    mu_type <- args$mutype    
    side <- args$side
    Xseed <- args$Xseed
    params <- data.frame(n, p, pi1, mu_type, side, Xseed)
} else {
    n <- 105
    p <- 100
    pi1 <- 0.3
    mu_type <- 2
    side <- "two"
    Xseed <- 2020
    params <- data.frame(n, p, pi1, mu_type, side, Xseed)
}

for (i in 1:nrow(params)){
    n <- params[i, 1]
    p <- params[i, 2]
    pi1 <- params[i, 3]
    mutype <- params[i, 4]
    side <- params[i, 5]
    Xseed <- params[i, 6]
    filename <- paste0("../data/BH_calib_lm", 
    	               "_n", n,
                       "_p", p,
                       "_pi1", pi1,
                       "_mutype", mutype,
                       "_side", side,
		       "_Xseed", Xseed,
                       ".RData")

    set.seed(Xseed)
    SigmaX <- genSigma(p, 0.8, "AR")
    X <- mvtnorm::rmvnorm(n, sigma = SigmaX)
    X <- scale(X, center = TRUE)
    X <- scale(X, scale = sqrt(colMeans(X^2)))
    mu1_ARplus_fix <- BH_lm_calib(X, pi1,
                                  "fix", mutype,
                                  side, nreps,
                                  alpha, target)

    set.seed(Xseed)
    SigmaX <- genSigma(p, 0, "iid")
    X <- mvtnorm::rmvnorm(n, sigma = SigmaX)
    X <- scale(X, center = TRUE)
    X <- scale(X, scale = sqrt(colMeans(X^2)))
    mu1_iid_fix <- BH_lm_calib(X, pi1,
                               "fix", mutype,
                               side, nreps,
                               alpha, target)

    set.seed(Xseed)
    SigmaX <- genSigma(p, 0.5, "block")
    X <- mvtnorm::rmvnorm(n, sigma = SigmaX)
    X <- scale(X, center = TRUE)
    X <- scale(X, scale = sqrt(colMeans(X^2)))
    mu1_block_fix <- BH_lm_calib(X, pi1,
                                 "fix", mutype,
                                 side, nreps,
                                 alpha, target)

    mu1 <- list(ARplus_fix = mu1_ARplus_fix,
                iid_fix = mu1_iid_fix,
                block_fix = mu1_block_fix)
    print(mu1)
    save(mu1, file = filename)
}
