source("dBH_utils.R")
source("utils.R")

BH_mvgauss_calib <- function(n, pi1,
                             mu_posit_type, mu_size_type,
                             rho, Sigma_type,
                             side,
                             nreps = 1000,
                             alpha = 0.05,
                             target = 0.3){
    Sigma <- genSigma(n, rho, Sigma_type)
    eigSigma <- eigen(Sigma)
    sqrtSigma <- with(eigSigma, vectors %*% (sqrt(values) * t(vectors)))

    mu_list <- lapply(1:nreps, function(i){
        mu <- genmu(n, pi1, 1, mu_posit_type, mu_size_type)
        if (side == "right"){
            mu <- abs(mu)
        } else if (side == "left"){
            mu <- -abs(mu)
        }
        return(mu)
    })
    null_zvals_list <- lapply(1:nreps, function(i){
        as.numeric(sqrtSigma %*% rnorm(n))
    })
    
    BH_power <- function(mu1){
        power <- sapply(1:nreps, function(i){
            H0 <- mu_list[[i]] == 0            
            mu <- mu_list[[i]] * mu1
            zvals <- null_zvals_list[[i]] + mu
            pvals <- pvals_mvgauss(zvals, Sigma, side)
            rejs_BH <- BH(pvals, alpha, 1:n, FALSE)$rejs
            tmp <- FDPpower(rejs_BH, H0)
            tmp[2]
        })
        mean(power) - target
    }

    lower <- 0
    upper <- 10
    while (TRUE & upper < 1000){
        tmp <- try(uniroot(BH_power, c(lower, upper))$root)
        if (class(tmp) == "try-error"){
            upper <- upper * 2
        } else {
            return(tmp)
        }
    }
    return(NA)
}

BH_mvt_calib <- function(n, df, pi1,
                         mu_posit_type, mu_size_type,
                         rho, Sigma_type,
                         side,
                         nreps = 1000,
                         alpha = 0.05,
                         target = 0.3){
    Sigma <- genSigma(n, rho, Sigma_type)
    eigSigma <- eigen(Sigma)
    sqrtSigma <- with(eigSigma, vectors %*% (sqrt(values) * t(vectors)))

    mu_list <- lapply(1:nreps, function(i){
        mu <- genmu(n, pi1, 1, mu_posit_type, mu_size_type)
        if (side == "right"){
            mu <- abs(mu)
        } else if (side == "left"){
            mu <- -abs(mu)
        }
        return(mu)
    })
    null_zvals_list <- lapply(1:nreps, function(i){
        as.numeric(sqrtSigma %*% rnorm(n))
    })
    
    BH_power <- function(mu1){
        power <- sapply(1:nreps, function(i){
            H0 <- mu_list[[i]] == 0            
            mu <- mu_list[[i]] * mu1
            zvals <- null_zvals_list[[i]] + mu
            sigmahat <- sqrt(rchisq(1, df = df) / df)
            tvals <- zvals / sigmahat
            pvals <- pvals_mvt(tvals, Sigma, df, side)
            rejs_BH <- BH(pvals, alpha, 1:n, FALSE)$rejs
            tmp <- FDPpower(rejs_BH, H0)
            tmp[2]
        })
        mean(power) - target
    }

    lower <- 0
    upper <- 10
    while (TRUE & upper < 1000){
        tmp <- try(uniroot(BH_power, c(lower, upper))$root)
        if (class(tmp) == "try-error"){
            upper <- upper * 2
        } else {
            return(tmp)
        }
    }
    return(NA)
}

BH_lm_calib <- function(X, pi1,
                        mu_posit_type, mu_size_type,
                        side,
                        nreps = 1000,
                        alpha = 0.05,
                        target = 0.3){
    n <- nrow(X)
    p <- ncol(X)
    Sigma <- solve(t(X) %*% X)
    H <- X %*% Sigma %*% t(X)
    df <- n - p

    beta_list <- lapply(1:nreps, function(i){
        beta <- genmu(p, pi1, 1, mu_posit_type, mu_size_type)
        if (side == "right"){
            beta <- abs(beta)
        } else if (side == "left"){
            beta <- -abs(beta)
        }
        return(beta)
    })
    eps_list <- lapply(1:nreps, function(i){
        rnorm(n)
    })
    
    BH_power <- function(mu1){
        power <- sapply(1:nreps, function(i){
            H0 <- beta_list[[i]] == 0            
            beta <- beta_list[[i]] * mu1
            eps <- eps_list[[i]]
            y <- X %*% beta + eps
            zvals <- Sigma %*% (t(X) %*% y)
            tmp <- as.numeric(t(y) %*% H %*% y)
            sigmahat <- sqrt((sum(y^2) - tmp) / df)
            tvals <- zvals / sigmahat
            pvals <- pvals_mvt(tvals, Sigma, df, side)
            rejs_BH <- BH(pvals, alpha, 1:n, FALSE)$rejs
            tmp <- FDPpower(rejs_BH, H0)
            tmp[2]
        })
        mean(power) - target
    }

    lower <- 0
    upper <- 10
    while (TRUE & upper < 1000){
        tmp <- try(uniroot(BH_power, c(lower, upper))$root)
        if (class(tmp) == "try-error"){
            upper <- upper * 2
        } else {
            return(tmp)
        }
    }
    return(NA)
}

BH_mcc_calib <- function(ng, nr, pi1,
                         mu_size_type,
                         side,
                         nreps = 1000,
                         alpha = 0.05,
                         target = 0.3){
    Sigma <- (diag(1, ng) + 1) / 2    
    df <- (ng + 1) * (nr - 1)

    mu_list <- lapply(1:nreps, function(i){
        mu <- genmu(ng, pi1, 1, "fix", mu_size_type)
        if (side == "right"){
            mu <- abs(mu)
        } else if (side == "left"){
            mu <- -abs(mu)
        }
        return(mu)
    })
    null_zvals_list <- lapply(1:nreps, function(i){
        matrix(rnorm((ng + 1) * nr), nrow = ng + 1)
    })
    
    BH_power <- function(mu1){
        power <- sapply(1:nreps, function(i){
            H0 <- mu_list[[i]] == 0            
            mu <- c(mu_list[[i]] * mu1, 0)
            tmp <- null_zvals_list[[i]] + mu
            zvals <- rowMeans(tmp)
            sigmahat <- sqrt(mean(rowMeans(tmp^2) - zvals^2) * nr / (nr - 1))
            zvals <- head(zvals, -1) - tail(zvals, 1)
            zvals <- zvals * sqrt(nr / 2)
            tvals <- zvals / sigmahat
            pvals <- pvals_mvt(tvals, Sigma, df, side)
            rejs_BH <- BH(pvals, alpha, 1:ng, FALSE)$rejs
            tmp <- FDPpower(rejs_BH, H0)
            tmp[2]
        })
        mean(power) - target
    }

    lower <- 0
    upper <- 10
    while (TRUE & upper < 1000){
        tmp <- try(uniroot(BH_power, c(lower, upper))$root)
        if (class(tmp) == "try-error"){
            upper <- upper * 2
        } else {
            return(tmp)
        }
    }
    return(NA)
}
