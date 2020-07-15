source("dBH_mvgauss.R")
source("dBH_mvt.R")
source("utils.R")
source("dBH_utils.R")
source("expr_functions.R")

## check_knots_mvgauss/mvt/lm tests the knots
## Here R can either be R_BH(alpha0) or R_dBH using
## R_BH(alpha0) as the initializer. R_BH is computed exactly
## using the homotopy algorithm. R_dBH is computed approximately
## based on a grid of points.
##
## Inputs:
##    id: the index of the test statistic
##    zvals, Sigma, sigmahat, df, X, y: inputs for
##        dBH_mvgauss, dBH_mvt and dBH_lm
##    side: "right" or "left" or "two"
##    alpha: target FDR level
##    avals: a values in the reshaping function
##    alpha0: FDR level used in dBH
##    low, high: left and right endpoints
##    dBH2: TRUE if R_dBH is computed; FALSE if R_BH is computed
##    gridfun: the function to compute a grid; see dBH_XXX_grid
##    gridsize: the size of the grid
##    knots_test: indicating whether the knots are checked
##        numerically; used only for R_BH
##    visual_test: indicating whether a visual test of knots
##        by comparing them with brute-force computation;
##        used only for R_BH
check_knots_mvgauss <- function(id, zvals, Sigma,
                                side = c("right", "left", "two"),
                                alpha = 0.05,
                                alpha0 = NULL,
                                avals = NULL,
                                avals_type = c("BH", "geom", "bonf", "manual"),
                                gamma = 2,

                                low = NULL, high = NULL){
    n <- length(zvals)
    avals_type <- avals_type[1]
    
    if (is.null(avals)){
        avals <- switch(avals_type,
                        BH = 1:n,
                        geom = geom_avals(gamma, n),
                        bonf = 1,
                        manual = NULL)
    }
    if (is.null(alpha0)){
        alpha0 <- alpha / normalize(avals)
    }
    side <- side[1]
    ntails <- ifelse(side == "two", 2, 1)
    
    vars <- diag(Sigma)
    zvals <- zvals / sqrt(vars)
    Sigma <- cov2cor(Sigma)
    thresh <- qnorm(1 - alpha0 * avals / n / ntails)
    if (is.null(low)){
        low <- qnorm(1 - alpha0 * max(avals) / n / ntails)
    }
    if (is.null(high)){
        high <- qnorm(1 - alpha / n / ntails * 0.01)
    }

    cor <- Sigma[-id, id]
    s <- zvals[-id] - zvals[id] * cor

    res <- compute_knots_mvgauss(
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
    
    ## Test the validity of each knot
    flag <- TRUE
    for (j in 1:ntails){
        knots <- res[[j]]$knots
        sgn <- res[[j]]$sgn
        hyp <- res[[j]]$hyp
        posit <- res[[j]]$posit + 1
        RCV <- res[[j]]$RCV
        m <- length(knots)
        if (m == 0){
            next
        }
        for (i in 1:m){
            knot <- knots[i]
            hyid <- hyp[i]
            pos <- posit[i]
            zvals0 <- recover_stats_mvgauss(zvals[id], knot, s, cor)
            if (side == "two"){
                zvals0 <- abs(zvals0)
            }
            if (abs(zvals0[hyid] - thresh[pos]) > 1e-10){
                flag <- FALSE
                break
            }
            zvals1 <- recover_stats_mvgauss(zvals[id], knot-1e-6, s, cor)
            zvals2 <- recover_stats_mvgauss(zvals[id], knot+1e-6, s, cor)
            if (side == "two"){
                zvals1 <- abs(zvals1)
                zvals2 <- abs(zvals2)
            }
            if ((zvals2[hyid] - zvals1[hyid]) * sgn[i] < 0){
                flag <- FALSE
                break
            }
        }
    }
    if (flag){
        print("The computed R curve is correct")
    } else {
        data <- list(zvals = zvals, Sigma = Sigma,
                     id = id, side = side,
                     alpha = alpha, avals = avals,
                     alpha0 = alpha0,
                     low = low, high = high)
        save(res, file = "../data/bad_data_mvgauss.RData")
        stop("The computed R curve is incorrect! This data is stored in ../data/bad_data_mvgauss.RData. Please debug compute_knots_mvgauss before proceeding.")
    }

    ## Visual test
    res <- lapply(res, function(re){
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
    
    par(mfrow = c(ntails, 1))
    for (j in 1:ntails){
        re <- res[[j]]
        RBH_dumb <- rep(NA, 10000)
        zlist <- seq(min(re$knots), max(re$knots), length.out = 10000) * (-1)^(j-1)

        ## Brute-force computation of number of rejections
        for (i in 1:length(zlist)){
            zvals1 <- recover_stats_mvgauss(zvals[id], zlist[i], s, cor)
            if (ntails == 2){
                pvals1 <- 2 * pnorm(abs(zvals1), lower.tail = FALSE)
            } else if (ntails == 1){
                pvals1 <- pnorm(zvals1, lower.tail = FALSE)
            }
            BH_res <- BH(pvals1, alpha0, avals, FALSE)
            RBH_dumb[i] <- BH_res$nrejs + !(1 %in% BH_res$rejs)
        }

        ## Compare two curves
        m <- length(re$nrejs)
        rex <- rep(re$knots, each = 2)
        rex[seq(1, length(rex), 2)] <- rex[seq(1, length(rex), 2)] - 1e-6
        rex[seq(2, length(rex), 2)] <- rex[seq(1, length(rex), 2)] + 1e-6
        rex <- rex[-c(1, length(rex))]
        rex <- rex * (-1)^(j-1)
        rey <- rep(re$nrejs, each = 2)

        plot(zlist, RBH_dumb, type = "l", lwd = 2)
        lines(rex, rey, col = "red", lty = 2, lwd = 2)
    }
}

check_knots_mvt <- function(id, zvals, Sigma, sigmahat, df,
                            side = c("right", "left", "two"),
                            alpha = 0.05,
                            alpha0 = NULL,
                            avals = NULL,
                            avals_type = c("BH", "geom", "bonf", "manual"),
                            gamma = 2,
                            low = NULL, high = NULL){
    n <- length(zvals)
    avals_type <- avals_type[1]
    
    if (is.null(avals)){
        avals <- switch(avals_type,
                        BH = 1:n,
                        geom = geom_avals(gamma, n),
                        bonf = 1,
                        manual = NULL)
    }
    if (is.null(alpha0)){
        alpha0 <- alpha / normalize(avals)
    }
    side <- side[1]
    ntails <- ifelse(side == "two", 2, 1)
    
    vars <- diag(Sigma)
    zvals <- zvals / sqrt(vars)
    Sigma <- cov2cor(Sigma)
    thresh <- qt(1 - alpha0 * avals / n /ntails, df = df)    
    if (is.null(low)){
        low <- qt(1 - alpha0 * max(avals) / n / ntails, df = df)        
    }
    if (is.null(high)){
        high <- qt(1 - alpha / n / ntails * 0.01, df = df)
    }
    tvals <- zvals / sigmahat
    
    cor <- Sigma[-id, id]
    s <- tvals[-id] - tvals[id] * cor

    res <- compute_knots_mvt(
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

    ## Test the validity of each knot    
    flag <- TRUE
    for (j in 1:ntails){
        knots <- res[[j]]$knots
        sgn <- res[[j]]$sgn
        hyp <- res[[j]]$hyp
        posit <- res[[j]]$posit + 1
        RCV <- res[[j]]$RCV
        m <- length(knots)
        if (m == 0){
            next
        }
        for (i in 1:m){
            knot <- knots[i]
            hyid <- hyp[i]
            pos <- posit[i]
            tvals0 <- recover_stats_mvt(tvals[id], knot, s, cor, df)
            if (side == "two"){
                tvals0 <- abs(tvals0)
            }
            if (abs(tvals0[hyid] - thresh[pos]) > 1e-10){
                flag <- FALSE
                break
            }
            tvals1 <- recover_stats_mvt(tvals[id], knot-1e-6, s, cor, df)
            tvals2 <- recover_stats_mvt(tvals[id], knot+1e-6, s, cor, df)
            if (side == "two"){
                tvals1 <- abs(tvals1)
                tvals2 <- abs(tvals2)
            }
            if ((tvals2[hyid] - tvals1[hyid]) * sgn[i] < 0){
                flag <- FALSE
                break
            }
        }
    }
    if (flag){
        print("The computed R curve is correct")
    } else {
        data <- list(zvals = zvals, Sigma = Sigma,
                     id = id, side = side,
                     alpha = alpha, avals = avals,
                     alpha0 = alpha0,
                     low = low, high = high)
        ## save(res, file = "../data/bad_data_mvt.RData")
        stop("The computed R curve is incorrect! This data is stored in ../data/bad_data_mvt.RData. Please debug compute_knots_mvt before proceeding.")
    }

    ## Visual test
    res <- lapply(res, function(re){
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
        ## All intervals either below thr or above thr
        ## For stability, we compare thr with the midpoint
        nrejs <- nrejs + ((knots_lo + knots_hi) / 2 < thr)
        list(knots = knots, nrejs = nrejs)
    })
        
    par(mfrow = c(ntails, 1))
    for (j in 1:ntails){
        re <- res[[j]]
        RBH_dumb <- rep(NA, 10000)
        tlist <- seq(min(re$knots), max(re$knots), length.out = 10000) * (-1)^(j-1)

        for (i in 1:length(tlist)){
            tvals1 <- recover_stats_mvt(tvals[id], tlist[i], s, cor, df)
            if (ntails == 2){
                pvals1 <- 2 * pt(abs(tvals1), df = df, lower.tail = FALSE)
            } else if (ntails == 1){
                pvals1 <- pt(tvals1, df = df, lower.tail = FALSE)
            }
            BH_res <- BH(pvals1, alpha0, avals, FALSE)
            RBH_dumb[i] <- BH_res$nrejs + !(1 %in% BH_res$rejs)
        }
        
#### Compare two curves
        m <- length(re$nrejs)
        rex <- rep(re$knots, each = 2)
        rex[seq(1, length(rex), 2)] <- rex[seq(1, length(rex), 2)] - 1e-6
        rex[seq(2, length(rex), 2)] <- rex[seq(1, length(rex), 2)] + 1e-6
        rex <- rex[-c(1, length(rex))]
        rex <- rex * (-1)^(j-1)
        rey <- rep(re$nrejs, each = 2)

        plot(tlist, RBH_dumb, type = "l", lwd = 2)
        lines(rex, rey, col = "red", lty = 2)
    }
}

check_knots_lm <- function(id, y, X,
                           side = c("right", "left", "two"),
                           alpha = 0.05,
                           alpha0 = NULL,
                           avals = NULL,
                           avals_type = c("BH", "geom", "bonf", "manual"),
                           gamma = 2,
                           low = NULL, high = NULL){
    stats <- lm_mvt(y, X)
    check_knots_mvt(id, stats$zvals, stats$Sigma,
                    stats$sigmahat, stats$df,
                    side, alpha, alpha0, 
                    avals, avals_type, gamma,
                    low, high)
}
