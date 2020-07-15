source("dBH_lm.R")
source("utils.R")
source("expr_functions.R")
library("knockoff")

if (!file.exists("../data/HIV_data.RData")){
    source("HIV_preprocess.R")
}
load("../data/HIV_data.RData")

HIV_expr <- function(X, y, 
                     side,
                     alphas, 
                     alpha_fac = 0.8,
                     gamma = 2,
                     tautype = c("GF", "QC"),
                     skip_knockoff = TRUE,
                     skip_dBH2 = TRUE,
                     ...){
    ## Log-transform the drug resistance measurements.
    y <- log(y)

    ## Remove patients with missing measurements.
    missing <- is.na(y)
    y <- y[!missing]
    X <- X[!missing,]

    ## Remove predictors that appear less than 3 times.
    X <- X[,colSums(X) >= 3]

    ## Remove duplicate predictors.
    X <- X[,colSums(abs(cor(X)-1) < 1e-4) == 1]

    ## Get names
    genes <- colnames(X)

    ## Get stats
    nalphas <- length(alphas)
    n <- nrow(X)
    p <- ncol(X)
    Sigma <- solve(t(X) %*% X)
    H <- X %*% Sigma %*% t(X)
    df <- n - p
    zvals <- Sigma %*% (t(X) %*% y)
    tmp <- as.numeric(t(y) %*% H %*% y)
    sigmahat <- sqrt((sum(y^2) - tmp) / df)
    pvals <- pvals_mvt(zvals, Sigma, sigmahat, df, side)
    
    ## Get methods and experimental settings for dBH
    methods <- gen_methods(alpha_fac, gamma, tautype,
                           skip_knockoff, skip_dBH2)
    expr_params <- expand.grid(
        alpha_fac = alpha_fac,
        gamma = gamma,
        tautype = tautype
    )

    if (!skip_knockoff){
        if (n < 2 * p){
            obj <- lm(y ~ X)
            sigma <- summary(obj)$sigma
            Xnew <- matrix(0, nrow = 2 * p - n, ncol = p)
            ynew <- rnorm(2 * p - n) * sigma
            X <- rbind(X, Xnew)
            y <- c(y, ynew)
        }
        print("n < 2p. Approximate knockoff is used.")
        Xk_equi <- knockoff::create.fixed(X, "equi")$Xk
        Xk_sdp <- knockoff::create.fixed(X, "sdp")$Xk
    }
    
    res <- list()
    pb <- txtProgressBar(style=3)
    for (k in 1:nalphas){
        alpha <- alphas[k]
        res[[k]] <- list(alpha = alpha, rejs = list())
        obj <- list()

        ## BH rejections
        for (x in union(NA, gamma)){
            if (is.na(x)){
                avals <- 1:p
            } else {
                avals <- geom_avals(x, p)
            }
            rejs_BH <- BH(pvals, alpha, avals, FALSE)
            rejs_BH_safe <- BH(pvals, alpha, avals, TRUE)
            obj <- c(obj, list(rejs_BH, rejs_BH_safe))
        }

        ## BC rejections
        rejs_BC <- BC(pvals, alpha)
        obj <- c(obj, list(rejs_BC))

        ## Knockoff rejections
        if (!skip_knockoff) {        
            rejs_knockoff_equi <- knockoff_lm(X, Xk_equi, y, alpha)
            rejs_knockoff_sdp <- knockoff_lm(X, Xk_sdp, y, alpha)
            obj <- c(obj, list(rejs_knockoff_equi, rejs_knockoff_sdp))
        }
        
        ## Number of methods so far
        nBHBCkn <- length(obj)
        
        ## dBH rejections
        for (j in 1:nrow(expr_params)){
            fac <- expr_params[j, 1]
            x <- expr_params[j, 2]
            type <- expr_params[j, 3]
            if (is.na(x)){
                avals_type <- "BH"
            } else {
                avals_type <- "geom"
            }
            if (is.na(fac)){
                alpha0 <- NULL
            } else {
                alpha0 <- fac * alpha
            }
            rejs_dBH <- dBH_mvt(
                zvals, Sigma, sigmahat, df,
                side, alpha,
                alpha0 = alpha0, 
                niter = 1,
                tautype = type,
                avals_type = avals_type,
                gamma = x, ...)
            rejs_dBH_init <- list(rejs = rejs_dBH$initrejs)
            obj <- c(obj, list(rejs_dBH, rejs_dBH_init))
        }

        if (!skip_dBH2){
            ## dBH2 rejections
            for (j in 1:nrow(expr_params)){
                fac <- expr_params[j, 1]
                x <- expr_params[j, 2]
                type <- expr_params[j, 3]
                if (is.na(x)){
                    avals_type <- "BH"
                } else {
                    avals_type <- "geom"
                }
                if (is.na(fac)){
                    alpha0 <- NULL
                } else {
                    alpha0 <- fac * alpha
                }
                rejs_dBH2 <- dBH_mvt(
                    zvals, Sigma, sigmahat, df,
                    side, alpha,
                    alpha0 = alpha0,
                    niter = 2,
                    tautype = type,
                    avals_type = avals_type,
                    gamma = x, ...)
                rejs_dBH2_init <- list(rejs = rejs_dBH2$initrejs)
                obj <- c(obj, list(rejs_dBH2, rejs_dBH2_init))
                setTxtProgressBar(pb, ((k - 1) * nrow(expr_params) + j) / (nalphas * nrow(expr_params)))
            }
        }
        names(obj) <- methods
        res[[k]]$rejs <- obj
    }
    
    return(res)
}

set.seed(20200711)
alpha_fac <- c(0.9, NA)
gamma <- c(NA, 2)
alphas <- c(0.05, 0.1, 0.2)
tautype <- "QC"
side <- "two"
skip_knockoff <- FALSE
skip_dBH2 <- FALSE
res <- list()

for (drug_class in names(data)[-1]){
    Y <- data[[drug_class]]$Y
    X <- data[[drug_class]]$X
    res[[drug_class]] <- lapply(1:length(Y), function(j){
        HIV_expr(X, Y[[j]],
                 side, alphas,
                 alpha_fac,
                 gamma,
                 tautype,
                 skip_knockoff,
                 skip_dBH2)
        print(j)
    })
}
save(file = "../data/HIV_res.RData", res)
