source("dBH_utils.R")
source("utils.R")

genSigma <- function(n, rho = 0,
                     type = c("AR", "MA", "equi", "block", "iid"),
                     bsize = 10){
    type <- type[1]
    if (type == "AR"){
        Sigma <- rho^(abs(outer(1:n, 1:n, "-")))
    } else if (type == "MA"){
        Sigma <- diag(n)
        Sigma[cbind(1:(n-1), 2:n)] <- rho
        Sigma[cbind(2:n, 1:(n-1))] <- rho
    } else if (type == "equi"){
        Sigma <- matrix(rho, n, n)
        diag(Sigma) <- 1
    } else if (type == "block"){
        m <- floor(n / bsize)
        if (n > bsize * m){
            warning("n is not divisible by bsize")
        }
        Sigma <- diag(n)
        blockSigma <- matrix(rho, bsize, bsize)
        diag(blockSigma) <- 1
        ## inds1 <- 1:ceiling(bsize / 2)
        ## inds2 <- ceiling(bsize / 2 + 1):bsize
        ## blockSigma[inds1, inds2] <- -rho
        ## blockSigma[inds2, inds1] <- -rho        
        for (i in 1:m){
            inds <- ((i - 1) * bsize + 1):(i * bsize)
            Sigma[inds, inds] <- blockSigma
        }
    } else if (type == "iid"){
        Sigma <- diag(n)
    }
    return(Sigma)
}


genmu <- function(n, pi1, mu1,
                  posit_type = c("random", "fix"),
                  mu_type = 1:3){
    m <- ceiling(n * pi1)
    posit_type <- posit_type[1]
    mu_type <- mu_type[1]
    if (posit_type == "random"){
        ## inds <- sample(n, m)
        inds <- seq(1, n, floor(1 / pi1))[1:m]
    } else if (posit_type == "fix"){
        inds <- 1:m
    }
    mu <- rep(0, n)
    altmu <- switch(mu_type,
                    `1` = rep(1, m),
                    `2` = rnorm(m),
                    `3` = rep(1, m) + 0.15 * (2 * rbinom(m, 1, 0.5) - 1))
    mu[inds] <- mu1 * altmu
    mu
}

gen_methods <- function(gamma,
                        beta,
                        tautype,
                        skip_knockoff,
                        skip_dBH2){
    expr_params <- expand.grid(
        gamma = gamma,
        beta = beta,
        tautype = tautype
    )

    BH_methods <- sapply(union(NA, beta), function(x){
        if (is.na(x)){
            beta <- "full"
        } else {
            beta <- paste0("sparse(", x, ")")
        }
        tmp <- paste0("BH_", beta)
        c(tmp, paste0(tmp, "_safe"))
    })
    methods <- c(BH_methods, "BC")
    if (!skip_knockoff){
        methods <- c(methods,
                     "Knockoff_equi",
                     "Knockoff_sdp")
    }
    dBH_methods <- apply(expr_params, 1, function(x){
        if (is.na(x[1])){
            gamma <- "safe"
        } else {
            gamma <- x[1]
        }
        if (is.na(x[2])){
            beta <- "full"
        } else {
            beta <- paste0("sparse(", as.numeric(x[2]), ")")
        }
        tautype <- x[3]
        
        method1 <- paste0("dBH_", tautype,
                          "_", beta,
                          "_", gamma)
        method2 <- paste0("dBH_init_", tautype,
                          "_", beta,
                          "_", gamma)
        c(method1, method2)
    })
    methods <- c(methods, as.character(dBH_methods))
    if (!skip_dBH2){
        dBH2_methods <- apply(expr_params, 1, function(x){
            if (is.na(x[1])){
                gamma <- "safe"
            } else {
                gamma <- x[1]
            }
            if (is.na(x[2])){
                beta <- "full"
            } else {
                beta <- paste0("sparse(", as.numeric(x[2]), ")")
            }
            tautype <- x[3]
            
            method1 <- paste0("dBH2_", tautype,
                              "_", beta,
                              "_", gamma)
            method2 <- paste0("dBH2_init_", tautype,
                              "_", beta,
                              "_", gamma)
            c(method1, method2)
        })
        methods <- c(methods,
                     as.character(dBH2_methods))
    }
    return(methods)
}

dBH_mvgauss_expr <- function(n, mu1, pi1,
                             mu_posit_type, mu_size_type,
                             rho, Sigma_type,
                             side,
                             alphas, nreps,
                             gamma = 0.9,
                             beta = 2,
                             tautype = "QC",
                             skip_dBH2 = TRUE,
                             ...){
    Sigma <- genSigma(n, rho, Sigma_type)

    nalphas <- length(alphas)
    eigSigma <- eigen(Sigma)
    sqrtSigma <- with(eigSigma, vectors %*% (sqrt(values) * t(vectors)))

    methods <- gen_methods(gamma, beta, tautype,
                           TRUE, skip_dBH2)
    expr_params <- expand.grid(
        gamma = gamma,
        beta = beta,
        tautype = tautype
    )

    results <- lapply(1:nalphas, function(k){
        tmp <- matrix(NA, length(methods), nreps)
        rownames(tmp) <- methods
        return(list(alpha = alphas[k],
                    FDP = tmp,
                    power = tmp,
                    secBH = tmp))
    })

    pb <- txtProgressBar(style=3)
    for (i in 1:nreps){
        mu <- genmu(n, pi1, mu1, mu_posit_type, mu_size_type)
        if (side == "right"){
            mu <- abs(mu)
        } else if (side == "left"){
            mu <- -abs(mu)
        }
        H0 <- mu == 0
        zvals <- as.numeric(mu + sqrtSigma %*% rnorm(n))
        pvals <- pvals_mvgauss(zvals, Sigma, side)

        for (k in 1:nalphas){
            obj <- list()
            alpha <- alphas[k]
            
            ## BH rejections
            for (x in union(NA, beta)){
                if (is.na(x)){
                    avals <- 1:n
                } else {
                    avals <- geom_avals(x, n)
                }
                rejs_BH <- BH(pvals, alpha, avals, FALSE)
                rejs_BH_safe <- BH(pvals, alpha, avals, TRUE)
                obj <- c(obj, list(rejs_BH, rejs_BH_safe))
            }

            ## BC rejections
            rejs_BC <- BC(pvals, alpha)
            obj <- c(obj, list(rejs_BC))

            ## Number of methods so far
            nBHBC <- length(obj)
            
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
                    gamma <- NULL
                } else {
                    gamma <- fac
                }
                rejs_dBH <- dBH_mvgauss(
                    zvals, Sigma, side, alpha,
                    gamma = gamma, 
                    niter = 1,
                    tautype = type,
                    avals_type = avals_type,
                    beta = x, ...)
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
                        gamma <- NULL
                    } else {
                        gamma <- fac
                    }
                    rejs_dBH2 <- dBH_mvgauss(
                        zvals, Sigma, side, alpha,
                        gamma = gamma, 
                        niter = 2,
                        tautype = type,
                        avals_type = avals_type,
                        beta = x, ...)
                    rejs_dBH2_init <- list(rejs = rejs_dBH2$initrejs)
                    obj <- c(obj, list(rejs_dBH2, rejs_dBH2_init))
                }
            }

            res <- sapply(obj, function(output){
                FDPpower(output$rejs, H0)
            })
            results[[k]]$FDP[, i] <- as.numeric(res[1, ])
            results[[k]]$power[, i] <- as.numeric(res[2, ])
            inds <- seq(nBHBC + 1, length(methods), 2)
            results[[k]]$secBH[inds, i] <- sapply(obj[inds], function(output){
                output$secBH
            })
            setTxtProgressBar(pb, ((i-1)*nalphas + k)/(nreps*nalphas))
        }
    }

    close(pb)
    return(results)
}

dBH_mvt_expr <- function(n, df, mu1, pi1,
                         mu_posit_type, mu_size_type,
                         rho, Sigma_type,
                         side,
                         alphas, nreps,
                         gamma = 0.9,
                         beta = 2,
                         tautype = "QC",
                         skip_dBH2 = TRUE,
                         ...){
    nalphas <- length(alphas)
    Sigma <- genSigma(n, rho, Sigma_type)
    eigSigma <- eigen(Sigma)
    sqrtSigma <- with(eigSigma, vectors %*% (sqrt(values) * t(vectors)))

    methods <- gen_methods(gamma, beta, tautype,
                           TRUE, skip_dBH2)
    expr_params <- expand.grid(
        gamma = gamma,
        beta = beta,
        tautype = tautype
    )

    results <- lapply(1:nalphas, function(k){
        tmp <- matrix(NA, length(methods), nreps)
        rownames(tmp) <- methods
        return(list(alpha = alphas[k],
                    FDP = tmp,
                    power = tmp,
                    secBH = tmp))
    })
    
    pb <- txtProgressBar(style=3)
    for (i in 1:nreps){
        mu <- genmu(n, pi1, mu1, mu_posit_type, mu_size_type)
        if (side == "right"){
            mu <- abs(mu)
        } else if (side == "left"){
            mu <- -abs(mu)
        }
        H0 <- mu == 0
        zvals <- as.numeric(mu + sqrtSigma %*% rnorm(n))
        sigmahat <- sqrt(rchisq(1, df = df) / df)
        pvals <- pvals_mvt(zvals, Sigma, sigmahat, df, side)

        for (k in 1:nalphas){
            obj <- list()
            alpha <- alphas[k]            

            ## BH rejections
            for (x in union(NA, beta)){
                if (is.na(x)){
                    avals <- 1:n
                } else {
                    avals <- geom_avals(x, n)
                }
                rejs_BH <- BH(pvals, alpha, avals, FALSE)
                rejs_BH_safe <- BH(pvals, alpha, avals, TRUE)
                obj <- c(obj, list(rejs_BH, rejs_BH_safe))
            }

            ## BC rejections
            rejs_BC <- BC(pvals, alpha)
            obj <- c(obj, list(rejs_BC))

            ## Number of methods so far
            nBHBC <- length(obj)

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
                    gamma <- NULL
                } else {
                    gamma <- fac
                }
                rejs_dBH <- dBH_mvt(
                    zvals, Sigma, sigmahat, df,
                    side, alpha,
                    gamma = gamma, 
                    niter = 1,
                    tautype = type,
                    avals_type = avals_type,
                    beta = x, ...)
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
                        gamma <- NULL
                    } else {
                        gamma <- fac
                    }
                    rejs_dBH2 <- dBH_mvt(
                        zvals, Sigma, sigmahat, df,
                        side, alpha,
                        gamma = gamma,
                        niter = 2,
                        tautype = type,
                        avals_type = avals_type,
                        beta = x, ...)
                    rejs_dBH2_init <- list(rejs = rejs_dBH2$initrejs)
                    obj <- c(obj, list(rejs_dBH2, rejs_dBH2_init))
                }
            }

            res <- sapply(obj, function(output){
                FDPpower(output$rejs, H0)
            })
            results[[k]]$FDP[, i] <- as.numeric(res[1, ])
            results[[k]]$power[, i] <- as.numeric(res[2, ])
            inds <- seq(nBHBC + 1, length(methods), 2)
            results[[k]]$secBH[inds, i] <- sapply(obj[inds], function(output){
                output$secBH
            })
            setTxtProgressBar(pb, ((i-1)*nalphas + k)/(nreps*nalphas))
        }
    }

    close(pb)
    return(results)
}

dBH_lm_expr <- function(X, mu1, pi1,
                        mu_posit_type, mu_size_type,
                        side,
                        alphas, nreps,
                        gamma = 0.9,
                        beta = 2,
                        tautype = "QC",
                        skip_knockoff = TRUE,
                        skip_dBH2 = TRUE,
                        ...){
    nalphas <- length(alphas)    
    n <- nrow(X)
    p <- ncol(X)
    if (n < 2 * p){
        skip_knockoff <- TRUE
    }
    if (!skip_knockoff){
        Xk_equi <- knockoff::create.fixed(X, "equi")$Xk
        Xk_sdp <- knockoff::create.fixed(X, "sdp")$Xk
    }

    Sigma <- solve(t(X) %*% X)
    H <- X %*% Sigma %*% t(X)
    df <- n - p

    methods <- gen_methods(gamma, beta, tautype,
                           skip_knockoff, skip_dBH2)
    expr_params <- expand.grid(
        gamma = gamma,
        beta = beta,
        tautype = tautype
    )

    results <- lapply(1:nalphas, function(k){
        tmp <- matrix(NA, length(methods), nreps)
        rownames(tmp) <- methods
        return(list(alpha = alphas[k],
                    FDP = tmp,
                    power = tmp,
                    secBH = tmp))
    })

    pb <- txtProgressBar(style=3)
    for (i in 1:nreps){
        beta <- genmu(p, pi1, 1, mu_posit_type, mu_size_type)
        if (side == "right"){
            beta <- abs(beta)
        } else if (side == "left"){
            beta <- -abs(beta)
        }
        beta <- beta * mu1
        H0 <- beta == 0

        eps <- rnorm(n)
        y <- X %*% beta + eps
        zvals <- Sigma %*% (t(X) %*% y)
        tmp <- as.numeric(t(y) %*% H %*% y)
        sigmahat <- sqrt((sum(y^2) - tmp) / df)
        pvals <- pvals_mvt(zvals, Sigma, sigmahat, df, side)

        for (k in 1:nalphas){
            obj <- list()
            alpha <- alphas[k]            

            ## BH rejections
            for (x in union(NA, beta)){
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
                    beta = x, ...)
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
                        beta = x, ...)
                    rejs_dBH2_init <- list(rejs = rejs_dBH2$initrejs)
                    obj <- c(obj, list(rejs_dBH2, rejs_dBH2_init))
                }
            }

            res <- sapply(obj, function(output){
                FDPpower(output$rejs, H0)
            })
            results[[k]]$FDP[, i] <- as.numeric(res[1, ])
            results[[k]]$power[, i] <- as.numeric(res[2, ])
            inds <- seq(nBHBCkn + 1, length(methods), 2)
            results[[k]]$secBH[inds, i] <- sapply(obj[inds], function(output){
                output$secBH
            })
            setTxtProgressBar(pb, ((i-1)*nalphas + k)/(nreps*nalphas))
        }
    }

    close(pb)
    return(results)
}

dBH_mcc_expr <- function(ng, nr,
                         mu1, pi1,
                         mu_posit_type, mu_size_type,
                         side,
                         alphas, nreps,
                         gamma = 0.9,
                         beta = 2,
                         tautype = "QC",
                         skip_knockoff = TRUE,
                         skip_dBH2 = TRUE,
                         ...){
    df <- (ng + 1) * (nr - 1)
    Sigma <- (diag(1, ng) + 1) / 2
    
    X <- lapply(1:ng, function(i){diag(rep(1, nr))})
    X <- do.call(rbind, X)
    X <- rbind(X, matrix(0, nr, nr))
    nalphas <- length(alphas)
    n <- nrow(X)
    p <- ncol(X)
    if (!skip_knockoff){
        Xk_equi <- knockoff::create.fixed(X, "equi")$Xk
        Xk_sdp <- knockoff::create.fixed(X, "sdp")$Xk
    }

    methods <- gen_methods(gamma, beta, tautype,
                           skip_knockoff, skip_dBH2)
    expr_params <- expand.grid(
        gamma = gamma,
        beta = beta,
        tautype = tautype
    )

    results <- lapply(1:nalphas, function(k){
        tmp <- matrix(NA, length(methods), nreps)
        rownames(tmp) <- methods
        return(list(alpha = alphas[k],
                    FDP = tmp,
                    power = tmp,
                    secBH = tmp))
    })

    pb <- txtProgressBar(style=3)
    for (i in 1:nreps){
        mu <- genmu(ng, pi1, 1, "fix", mu_size_type) * mu1
        if (side == "right"){
            mu <- abs(mu)
        } else if (side == "left"){
            mu <- -abs(mu)
        }
        H0 <- mu == 0
        y <- matrix(rnorm((ng + 1) * nr), nrow = ng + 1) + c(mu, 0)
        zvals <- rowMeans(y)
        sigmahat <- sqrt(mean(rowMeans(y^2) - zvals^2) * nr / (nr - 1))
        zvals <- head(zvals, -1) - tail(zvals, 1)
        zvals <- zvals * sqrt(nr / 2)
        tvals <- zvals / sigmahat
        pvals <- pvals_mvt(zvals, Sigma, sigmahat, df, side)

        for (k in 1:nalphas){
            obj <- list()
            alpha <- alphas[k]            

            ## BH rejections
            for (x in union(NA, beta)){
                if (is.na(x)){
                    avals <- 1:ng
                } else {
                    avals <- geom_avals(x, ng)
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
                y <- as.numeric(y)
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
                    zvals, Sigma,
                    sigmahat, df,
                    side, alpha,
                    alpha0 = alpha0, 
                    niter = 1,
                    tautype = type,
                    avals_type = avals_type,
                    beta = x, ...)
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
                        zvals, Sigma,
                        sigmahat, df,
                        side, alpha,
                        alpha0 = alpha0,
                        niter = 2,
                        tautype = type,
                        avals_type = avals_type,
                        beta = x, ...)
                    rejs_dBH2_init <- list(rejs = rejs_dBH2$initrejs)
                    obj <- c(obj, list(rejs_dBH2, rejs_dBH2_init))
                }
            }

            res <- sapply(obj, function(output){
                FDPpower(output$rejs, H0)
            })
            results[[k]]$FDP[, i] <- as.numeric(res[1, ])
            results[[k]]$power[, i] <- as.numeric(res[2, ])
            inds <- seq(nBHBCkn + 1, length(methods), 2)
            results[[k]]$secBH[inds, i] <- sapply(obj[inds], function(output){
                output$secBH
            })
            setTxtProgressBar(pb, ((i-1)*nalphas + k)/(nreps*nalphas))
        }
    }

    close(pb)
    return(results)
}

postprocess <- function(res){
    summaryres <- lapply(res, function(re){
        FDR <- as.numeric(rowMeans(re$FDP))
        FDR <- round(FDR, 4)
        power <- as.numeric(rowMeans(re$power))
        secBH <- as.numeric(rowMeans(re$secBH))
        methods <- rownames(re$power)
        df <- data.frame(method = methods,
                         FDR = FDR,
                         power = power,
                         secBH = secBH)
        inds1 <- grep("^BH", methods)
        inds2 <- grep("^BC", methods)
        inds3 <- grep("^Knockoff", methods)
        inds <- c(inds1, inds2, inds3)
        df1 <- df[inds, ]
        df1[, 5:6] <- NA
        names(df1)[5:6] <- c("FDR (init)", "power (init)")
        df2 <- df[-inds, ]
        m <- nrow(df2)
        df2_1 <- df2[seq(1, m, 2), ]
        df2_2 <- df2[seq(2, m, 2), ][, 2:3]
        names(df2_2) <- c("FDR (init)", "power (init)")
        df2 <- cbind(df2_1, df2_2)

        df <- rbind(df1, df2)
        df <- df[, c(1, 2, 5, 3, 6, 4)]
        df$alpha <- re$alpha
        return(df)
    })
    do.call(rbind, summaryres)
}
