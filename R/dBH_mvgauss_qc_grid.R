source("dBH_utils.R")

dBH_mvgauss_qc_grid <- function(zvals, Sigma,
                                side = c("right", "left", "two"),
                                alpha = 0.05, alpha0 = NULL,
                                avals = NULL,
                                avals_type = c("BH", "geom", "bonf", "manual"),
                                gamma = 2,
                                eps = 0.05,
                                Rfun = dBH_mvgauss_qc,
                                gridfun = lingrid,
                                gridsize = 20,
                                if_thresh_expt = TRUE,
                                thresh_expt = min(0.9 * alpha, 1.1 * alpha0),
                                ...){
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
    pvals <- zvals_pvals(zvals, side)
    qvals <- qvals_BH_reshape(pvals, avals)
    
    params_root <- list(Sigma = Sigma, side = side,
                        alpha = alpha, alpha0 = alpha0,
                        avals = avals_obj$avals,
                        avals_type = avals_obj$avals_type,
                        gamma = avals_obj$gamma,
                        eps = eps, ...)
    params <- c(params_root, list(zvals = zvals))
    res_init <- do.call(Rfun, params)
    Rinit <- rep(length(res_init$initrejs) + 1, n)
    Rinit[res_init$initrejs] <- Rinit[res_init$initrejs] - 1

    ntails <- ifelse(side == "two", 2, 1)
    ## if (res_init$safe){
    ##     low <- qnorm(alpha / ntails, lower.tail = FALSE)
    ## } else {
    ##     low <- 0
    ## }
    high <- qnorm(alpha * eps / n / ntails, lower.tail = FALSE)

    cand <- res_init$cand
    if (res_init$safe){
        init_rejlist <- res_init$initrejs
    } else {
        if (if_thresh_expt){
            init_rejlist <- cand[res_init$expt <= thresh_expt]
        } else {
            init_rejlist <- integer(0)
        }
    }
    cand <- setdiff(cand, init_rejlist)    
    
    if (length(cand) == 0){
        return(list(rejs = init_rejlist,
                    initrejs = init_rejlist,
                    cand = numeric(0),
                    expt = numeric(0),
                    safe = res_init$safe,
                    secBH = FALSE))
    }

    cand_info <- sapply(cand, function(i){
        cor <- Sigma[-i, i]
        s <- zvals[-i] - cor * zvals[i]
        low <- qnorm(qvals[i] * max(avals) / n / ntails, lower.tail = FALSE)

        ## RBH function with alpha = qi
        res_q <- compute_knots_mvgauss(
            zstat = zvals[i],
            zminus = zvals[-i],
            cor = Sigma[-i, i],
            alpha = qvals[i],
            side = side,            
            low = low,
            high = high,
            avals = avals,
            avals_type = avals_type,
            gamma = gamma)
        res_q <- lapply(res_q, function(re){
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
            thr <- qnorm(thra * qvals[i] / n / ntails, lower.tail = FALSE)
            list(knots = knots, thr = thr)
        })

        ## Get the intervals of knots at which the numerator
        ## is 1
        grids <- lapply(1:ntails, function(k){
            ints <- find_int_above_thr(res_q[[k]]$knots, res_q[[k]]$thr)
            ## Filter out very small intervals
            if (k == 2){
                lapply(ints, function(int){
                    -rev(int)
                })
            } else {
                ints
            }
        })
        grids <- do.call(c, grids)
        prob <- sapply(grids, function(int){
            diff(pnorm(int))
        })
	if (length(prob) < 1 || !is.numeric(prob) || sum(prob) * n <= alpha){
            return(c(1, NA))
        }

        ## Add extra knots
        nknots <- ceiling(prob / sum(prob) * gridsize * ntails)
        grids <- lapply(1:length(grids), function(i){
            seq(grids[[i]][1], grids[[i]][2],
                length.out = nknots[i] + 1)
        })

        ## Create grid for the denominator
        expt <- sapply(grids, function(grid){
            prob <- diff(pnorm(grid))            
            ex <- sapply(1:(length(grid) - 1), function(j){
                pr <- prob[j]
                if (any(grid > 0)){
                    j <- j + 1
                }
                tmp <- recover_stats_mvgauss(zvals[i], grid[j], s, cor)
                zvals_tmp <- rep(0, n)
                zvals_tmp[i] <- tmp[1]
                zvals_tmp[-i] <- tmp[-1]
                params <- c(params_root, list(zvals = zvals_tmp))
                res <- do.call(Rfun, params)
                nrejs <- length(res$initrejs) + !(i %in% res$initrejs)
                return(pr / nrejs)
            })
            sum(ex)
        })
        expt <- sum(expt) * n
        ifrej <- expt <= alpha
        return(c(ifrej, expt))
    })

    ifrej <- as.logical(cand_info[1, ])
    rejlist <- which(ifrej)
    rejlist <- c(init_rejlist, cand[rejlist])
    expt <- cand_info[2, ]
    
    if (length(rejlist) == 0){
        return(list(rejs = numeric(0),
                    initrejs = numeric(0),
                    cand = cand,
                    expt = expt,
                    safe = res_init$safe,
                    secBH = FALSE))
    }

    rejlist <- sort(rejlist)
    Rplus <- length(rejlist)
    if (Rplus >= max(Rinit[rejlist])){
        return(list(rejs = rejlist,
                    initrejs = rejlist, 
                    cand = cand,
                    expt = expt,
                    safe = res_init$safe,
                    secBH = FALSE))
    }

    uvec <- runif(Rplus)
    secBH_fac <- Rinit[rejlist] / Rplus
    tdp <- uvec * secBH_fac
    nrejs <- nrejs_BH(tdp, 1)
    thr <- max(nrejs, 1) / Rplus
    secrejs <- which(tdp <= thr)
    rejs <- rejlist[secrejs]
    return(list(rejs = rejs,
                initrejs = rejlist,
                cand = cand,
                expt = expt,
                safe = FALSE,
                secBH = TRUE,
                secBH_fac = secBH_fac))
}
