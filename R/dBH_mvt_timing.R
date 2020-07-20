source("utils.R")
source("expr_functions.R")
library("dbh")

if (!interactive()){
    suppressPackageStartupMessages(library("argparse"))

    parser <- ArgumentParser()
    parser$add_argument("--nalt", type = "integer", default = 10, help = "Number of alternative hypotheses")
    parser$add_argument("--df", type = "integer", default = 50, help = "Degree of freedom")
    parser$add_argument("--rho", type = "double", default = 0.1, help = "Autocorrelation")
    parser$add_argument("--fac", type = "double", default = 1, help = "Mean factor (mu * sqrt(2log n))")
    parser$add_argument("--nreps", type = "integer", default = 20, help = "Number of repeats")
    parser$add_argument("--dBH2", type = "character", default = "TRUE", help = "Whether to include dBH2")
    parser$add_argument("--seed", type = "integer", default = 1, help = "Random seed")

    args <- parser$parse_args()

    nalt <- args$nalt
    df <- args$df
    rho <- args$rho
    fac <- args$fac
    nreps <- args$nreps
    skip_dBH2 <- !as.logical(args$dBH2)
    seed <- args$seed
} else {
    nalt <- 30
    df <- 50
    rho <- 0.8
    fac <- 1
    nreps <- 1
    skip_dBH2 <- FALSE
    seed <- 0
}

set.seed(seed)
filename <- paste0("../cluster_raw_data/dBH_mvt_timing", 
                   "_nalt", nalt,
                   "_df", df,
                   "_rho", rho,
                   "_fac", fac,
                   "_nreps", nreps,
                   "_dBH2", !skip_dBH2,
                   "_seed", seed,
                   ".RData")
alpha <- 0.05
gamma_list <- 1
avals_type_list <- c("BH", "geom")
side_list <- c("right", "two")

n1_list <- 10^(1:4)
n1_list <- n1_list[n1_list > nalt]
expr1_params <- expand.grid(n = n1_list,
                            side = side_list,
                            gamma = gamma_list,
                            avals_type = avals_type_list,
                            stringsAsFactors = FALSE)

res <- data.frame()

for (i in 1:nrow(expr1_params)){
    n <- expr1_params$n[i]
    gamma <- expr1_params$gamma[i]
    if (is.na(gamma)){
        gamma <- NULL
    }
    avals_type <- expr1_params$avals_type[i]
    side <- expr1_params$side[i]
    pi1 <- nalt / n
    mu1 <- fac * sqrt(2 * log(n))
    mu <- genmu(n, pi1, mu1, "fix")    
    zvals <- gen_fast_AR(n, rho) + mu
    Sigmafun <- function(i){
        rho^(abs(1:n - i))
    }
    sigmahat <- sqrt(rchisq(1, df = df) / df)
    tvals <- zvals / sigmahat
    runtime <- system.time(
        dBH_mvt(tvals = tvals,
                df = df,
                Sigmafun = Sigmafun,
                side = side,
                alpha = alpha,
                gamma = gamma,
                avals_type = avals_type,
                niter = 1,
                verbose = TRUE)
    )
    method <- ifelse(is.null(gamma), "dBY", "dBH")
    time <- as.numeric(runtime[3])
    tmp <- data.frame(n = n, niter = 1,
                      method = method,
                      side = side,
                      avals_type = avals_type,
                      time = time,
                      nalt = nalt,
                      df = df,
                      model = paste0("AR(", rho, ")"), 
                      fac = fac)
    print(tmp)
    res <- rbind(res, tmp)
    gc()	
    gc()
    save(res, file = filename)    
}

if (!skip_dBH2){
    n2_list <- 10^(1:4)
    n2_list <- n2_list[n2_list > nalt]
    expr2_params <- expand.grid(n = n2_list,
                                side = side_list,
                                gamma = gamma_list,
                                avals_type = avals_type_list,
                                stringsAsFactors = FALSE)
    
    for (i in 1:nrow(expr2_params)){
        n <- expr2_params$n[i]
        gamma <- expr2_params$gamma[i]
        if (is.na(gamma)){
            gamma <- NULL
        }
        avals_type <- expr2_params$avals_type[i]
        side <- expr2_params$side[i]
        pi1 <- nalt / n
        mu1 <- fac * sqrt(2 * log(n))
        mu <- genmu(n, pi1, mu1, "fix")    
        zvals <- gen_fast_AR(n, rho) + mu
        Sigmafun <- function(i){
            rho^(abs(1:n - i))
        }
        sigmahat <- sqrt(rchisq(1, df = df) / df)
        tvals <- zvals / sigmahat
        runtime <- system.time(
            dBH_mvt(tvals = tvals,
                    df = df,
                    Sigmafun = Sigmafun,
                    side = side,
                    alpha = alpha,
                    gamma = gamma,
                    avals_type = avals_type,
                    niter = 2,
                    verbose = TRUE)
        )
        method <- ifelse(is.null(gamma), "dBY", "dBH")
        time <- as.numeric(runtime[3])        
        tmp <- data.frame(n = n, niter = 2,
                         method = method,
                         side = side,
                         avals_type = avals_type,
                         time = time,
                         nalt = nalt,
                         df = df,
                         model = paste0("AR(", rho, ")"), 
                         fac = fac)
        print(tmp)
        res <- rbind(res, tmp)
	gc()	
	gc()
        save(res, file = filename)        
    }
}
