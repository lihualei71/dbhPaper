source("dBH_mvgauss.R")
source("dBH_mvt.R")
source("utils.R")
source("dBH_utils.R")
source("expr_functions.R")

Rcurve_mvgauss <- function(id, zvals, Sigma,
                           side = c("right", "left", "two"),
                           alpha = 0.05,
                           alpha0 = NULL,
                           tautype = c("GF", "QC"),
                           avals = NULL,
                           avals_type = c("BH", "geom", "bonf", "manual"),
                           gamma = 2,
                           low = NULL, high = NULL,
                           niter = 0,
                           gridsize = 100){
    side <- side[1]
    avals_type <- avals_type[1]
    n <- length(zvals) 
    avals_obj <- list(avals = avals,
                      avals_type = avals_type,
                      gamma = gamma)
   
    if (is.null(avals)){
        if (avals_type == "manual"){
            stop("avals must be inputted when avals_type = \"manual\"")
        } else if (avals_type == "geom" && gamma <= 1){
            stop("gamma must be larger than 1 when avals_type = \"geom\"")
        }
        avals <- switch(avals_type,
                        BH = 1:n,
                        geom = geom_avals(gamma, n),
                        bonf = 1)
    } else {
        if (avals[1] != 1){
            stop("The first element of avals must be 1.")
        }
        avals_type <- "manual"
        warning("avals is inputted and avals_type is set to be \"manual\" by default. This may slow down the code. Use the built-in avals_type (\"BH\", \"geom\" or \"bonf\") instead unless there is a good reason to use the inputted avals.")
    }

    if (is.null(alpha0)){
        alpha0 <- alpha / normalize(avals)
    }

    if (any(diag(Sigma) != 1)){
        vars <- diag(Sigma)
        zvals <- zvals / sqrt(vars)
        Sigma <- cov2cor(Sigma)
    }
    if (side == "left"){
        zvals <- -zvals
        side <- "one"
    } else if (side == "right"){
        side <- "one"
    }

    ntails <- ifelse(side == "two", 2, 1)
    pvals <- zvals_pvals(zvals, side)

    if (niter == 0){
        ## RBH function with alpha = alpha0
        res_alpha0 <- compute_knots_mvgauss(
            zstat = zvals[id],
            zminus = zvals[-id],
            cor = Sigma[-id, id],
            alpha = alpha0,
            side = side,            
            low = low, 
            high = high,
            avals = avals,
            avals_type = avals_type,
            gamma = gamma)
        res_alpha0 <- lapply(res_alpha0, function(re){
            RBH <- RejsBH(re$posit, re$sgn, re$RCV, avals)
            knots <- c(re$low, re$knots)
            RBH <- rle(RBH)
            nrejs <- RBH$values
            cutinds <- c(1, cumsum(RBH$lengths) + 1)
            knots <- c(knots, re$high)        
            knots <- knots[cutinds]
            if (knots[1] < 0){
                knots <- rev(abs(knots))
                ## This requires the null distribution to be symmetric
                nrejs <- rev(nrejs)
            }
            if (avals_type == "BH"){
                thra <- nrejs
            } else if (avals_type == "geom"){
                thra <- find_ind_geom_avals(gamma, nrejs, "max")
                ## 0 rejection should return aval = 0
                thra[thra == 0] <- NA
                thra <- avals[thra]
                thra[is.na(thra)] <- 0
            } else if (avals_type == "manual"){
                thra <- fill_int_general(nrejs, avals)
                ## 0 rejection should return aval = 0
                thra[thra == 0] <- NA
                thra <- avals[thra]
                thra[is.na(thra)] <- 0
            } else if (avals_type == "bonf"){
                thra <- rep(1, length(nrejs))
            }
            thr <- qnorm(thra * alpha0 / n / ntails, lower.tail = FALSE)
            knots_lo <- head(knots, -1)
            knots_hi <- tail(knots, -1)
            nrejs <- nrejs + ((knots_lo + knots_hi) / 2 < thr)
            list(knots = knots, nrejs = nrejs)
        })
    } else if (niter >= 1) {
        cor <- Sigma[-id, id]
        s <- zvals[-id] - cor * zvals[id]

        params_root <- list(Sigma = Sigma, side = side,
                            alpha = alpha, alpha0 = alpha0,
                            avals = avals_obj$avals,
                            avals_type = avals_obj$avals_type,
                            gamma = avals_obj$gamma)
        knots <- seq(low, high,
                     length.out = tail(gridsize, 1) + 1)
        knots_fun <- function(...){
            dBH_mvgauss(...,
                        niter = niter,
                        tautype = tautype,
                        gridsize = gridsize[1])
        }
        
        res_alpha0 <- lapply(1:ntails, function(k){
            nrejs <- sapply(2:(gridsize + 1), function(j){
                tmp <- recover_stats_mvgauss(zvals[id], knots[j] * (-1)^(k-1), s, cor)
                zvals_tmp <- rep(0, n)
                zvals_tmp[id] <- tmp[1]
                zvals_tmp[-id] <- tmp[-1]
                params <- c(params_root, list(zvals = zvals_tmp))
                res <- do.call(knots_fun, params)
                nrejs <- length(res$initrejs) + !(id %in% res$initrejs)
            })
            list(knots = knots, nrejs = nrejs)            
        })
    }
    res_alpha0 <- lapply(1:ntails, function(k){
        knots <- pnorm(res_alpha0[[k]]$knots, lower.tail = FALSE)
        list(knots = rev(knots), nrejs = rev(res_alpha0[[k]]$nrejs))
    })
    return(c(res_alpha0, list(alpha0 = alpha0)))
}

Rcurve_mvt <- function(id, zvals, Sigma, sigmahat, df,
                       side = c("right", "left", "two"),
                       alpha = 0.05,
                       alpha0 = NULL,
                       tautype = c("GF", "QC"),
                       avals = NULL,
                       avals_type = c("BH", "geom", "bonf", "manual"),
                       gamma = 2,
                       low = NULL, high = NULL,
                       niter = 0,
                       gridsize = 100){
    side <- side[1]    
    n <- length(zvals)
    avals_type <- avals_type[1]
    avals_obj <- list(avals = avals,
                      avals_type = avals_type,
                      gamma = gamma)

    if (is.null(avals)){
        if (avals_type == "manual"){
            stop("avals must be inputted when avals_type = \"manual\"")
        } else if (avals_type == "geom" && gamma <= 1){
            stop("gamma must be larger than 1 when avals_type = \"geom\"")
        }
        avals <- switch(avals_type,
                        BH = 1:n,
                        geom = geom_avals(gamma, n),
                        bonf = 1)
    } else {
        if (avals[1] != 1){
            stop("The first element of avals must be 1.")
        }
        avals_type <- "manual"
        warning("avals is inputted and avals_type is set to be \"manual\" by default. This may slow down the code. Use the built-in avals_type (\"BH\", \"geom\" or \"bonf\") instead unless there is a good reason to use the inputted avals.")
    }

    if (is.null(alpha0)){
        alpha0 <- alpha / normalize(avals)
    }

    if (any(diag(Sigma) != 1)){
        vars <- diag(Sigma)
        zvals <- zvals / sqrt(vars)
        Sigma <- cov2cor(Sigma)
    }
    if (side == "left"){
        zvals <- -zvals
        side <- "one"
    } else if (side == "right"){
        side <- "one"
    }

    ntails <- ifelse(side == "two", 2, 1)        
    tvals <- zvals / sigmahat
    high <- qt(alpha * eps / n / ntails, df = df, lower.tail = FALSE)
    pvals <- tvals_pvals(tvals, df, side)
    
    if (niter == 0){
        ## RBH function with alpha = alpha0
        res_alpha0 <- compute_knots_mvgauss(
            tstat = tvals[id],
            tminus = tvals[-id],
            df = df, 
            cor = Sigma[-id, id],
            alpha = alpha0,
            side = side,            
            low = low, 
            high = high,
            avals = avals,
            avals_type = avals_type,
            gamma = gamma)
        res_alpha0 <- lapply(res_alpha0, function(re){
            RBH <- RejsBH(re$posit, re$sgn, re$RCV, avals)
            knots <- c(re$low, re$knots)
            RBH <- rle(RBH)
            nrejs <- RBH$values
            cutinds <- c(1, cumsum(RBH$lengths) + 1)
            knots <- c(knots, re$high)        
            knots <- knots[cutinds]
            if (knots[1] < 0){
                knots <- rev(abs(knots))
                ## This requires the null distribution to be symmetric
                nrejs <- rev(nrejs)
            }
            if (avals_type == "BH"){
                thra <- nrejs
            } else if (avals_type == "geom"){
                thra <- find_ind_geom_avals(gamma, nrejs, "max")
                ## 0 rejection should return aval = 0
                thra[thra == 0] <- NA
                thra <- avals[thra]
                thra[is.na(thra)] <- 0
            } else if (avals_type == "manual"){
                thra <- fill_int_general(nrejs, avals)
                ## 0 rejection should return aval = 0
                thra[thra == 0] <- NA
                thra <- avals[thra]
                thra[is.na(thra)] <- 0
            } else if (avals_type == "bonf"){
                thra <- rep(1, length(nrejs))
            }
            thr <- qt(thra * alpha0 / n / ntails, df = df, lower.tail = FALSE)
            knots_lo <- head(knots, -1)
            knots_hi <- tail(knots, -1)
            nrejs <- nrejs + ((knots_lo + knots_hi) / 2 < thr)
            list(knots = knots, nrejs = nrejs)
        })
    } else if (niter >= 1) {
        cor <- Sigma[-id, id]
        s <- zvals[-id] - cor * zvals[id]

        params_root <- list(Sigma = Sigma, df = df, side = side,
                            alpha = alpha, alpha0 = alpha0,
                            avals = avals_obj$avals,
                            avals_type = avals_obj$avals_type,
                            gamma = avals_obj$gamma)
        knots <- seq(low, high, length.out = gridsize + 1)
        knots_fun <- function(...){
            dBH_mvt(..., niter = niter, tautype = tautype)
        }
        
        res_alpha0 <- lapply(1:ntails, function(k){
            nrejs <- sapply(2:(gridsize + 1), function(j){
                tmp <- recover_stats_mvt(zvals[id], knots[j] * (-1)^(k-1), s, cor)
                tvals_tmp <- rep(0, n)
                tvals_tmp[id] <- tmp[1]
                tvals_tmp[-id] <- tmp[-1]
                sigmahat_tmp <- sigmahat * sqrt((df + tvals[i]^2) / (df + knots[j]^2))
                zvals_tmp <- tvals_tmp * sigmahat_tmp
                params <- c(params_root, list(zvals = zvals_tmp, sigmahat = sigmahat_tmp))
                res <- do.call(knots_fun, params)
                nrejs <- length(res$initrejs) + !(id %in% res$initrejs)
            })
            list(knots = knots, nrejs = nrejs)            
        })
    }
    res_alpha0 <- lapply(1:ntails, function(k){
        knots <- pt(res_alpha0[[k]]$knots, lower.tail = FALSE, df = df)
        list(knots = rev(knots), nrejs = rev(res_alpha0[[k]]$nrejs))
    })
    return(c(res_alpha0, list(alpha0 = alpha0)))
}

Rcurve_lm <- function(id, y, X,
                      side = c("right", "left", "two"),
                      alpha = 0.05,
                      alpha0 = NULL,
                      tautype = c("GF", "QC"),
                      avals = NULL,
                      avals_type = c("BH", "geom", "bonf", "manual"),
                      gamma = 2,
                      low = NULL, high = NULL,
                      niter = 0,
                      gridsize = 100){
    stats <- lm_mvt(y, X)
    Rcurve_mvt(id, stats$zvals, stats$Sigma,
               stats$sigmahat, stats$df,
               side,
               alpha, alpha0,
               tautype,
               avals, avals_type, gamma,
               low, high,
               niter,
               gridsize)
}
