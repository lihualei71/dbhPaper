source("dBH_utils.R")

dBH_mvt_gf_grid <- function(zvals, Sigma,
                            sigmahat, df, 
                            side = c("right", "left", "two"),
                            alpha = 0.05, alpha0 = NULL,
                            avals = NULL,
                            avals_type = c("BH", "geom", "bonf", "manual"),
                            gamma = 2,
                            eps = 0.05,
                            Rfun = dBH_mvt_gf,
                            gridfun = lingrid,
                            gridsize = 20,
                            if_thresh_expt = TRUE,
                            thresh_expt = max(0.9 * alpha, alpha0),
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
    tvals <- zvals / sigmahat
    pvals <- tvals_pvals(tvals, df, side)
    dBH_rej0 <- BH(pvals, alpha0, avals, reshape = FALSE)$rejs

    params_root <- list(Sigma = Sigma, df = df, side = side,
                        alpha = alpha, alpha0 = alpha0,
                        avals = avals_obj$avals,
                        avals_type = avals_obj$avals_type,
                        gamma = avals_obj$gamma,
                        eps = eps, ...)
    params <- c(params_root, list(zvals = zvals, sigmahat = sigmahat))
    res_init <- do.call(Rfun, params)
    Rinit <- rep(length(res_init$initrejs) + 1, n)
    Rinit[res_init$initrejs] <- Rinit[res_init$initrejs] - 1

    ntails <- ifelse(side == "two", 2, 1)
    if (res_init$safe){
        low <- qt(alpha / ntails, df = df, lower.tail = FALSE)
    } else {
        low <- 0
    }
    high <- qt(alpha * eps / n / ntails, df = df, lower.tail = FALSE)

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
                    expt0 = numeric(0),
                    safe = res_init$safe,
                    secBH = FALSE))
    }

    cand_info <- sapply(cand, function(i){
        cor <- Sigma[-i, i]
        s <- tvals[-i] - cor * tvals[i]
        grid <- gridfun(min(low, abs(tvals[i])), high, gridsize, side)
        rejfun <- lapply(1:ntails, function(k){
            knots <- abs(grid[[k]])
            extraknot <- abs(tvals[i])
            if (!(extraknot %in% knots)){
                knots <- c(knots, extraknot)
            }
            knots <- sort(knots)
            nrejs <- sapply(2:length(knots), function(j){
                tmp <- recover_stats_mvt(tvals[i], knots[j] * (-1)^(k-1), s, cor, df)
                tvals_tmp <- rep(0, n)
                tvals_tmp[i] <- tmp[1]
                tvals_tmp[-i] <- tmp[-1]
                sigmahat_tmp <- sigmahat * sqrt((df + tvals[i]^2) / (df + knots[j]^2))
                zvals_tmp <- tvals_tmp * sigmahat_tmp

                params <- c(params_root, list(zvals = zvals_tmp, sigmahat = sigmahat_tmp))
                res <- do.call(Rfun, params)

                pvals_tmp <- tvals_pvals(tvals_tmp, df, side)
                ifrej <- (i %in% BH(pvals_tmp, alpha0, avals, FALSE)$rejs)
                nrejs <- length(res$initrejs) + !(i %in% res$initrejs)
                return(c(nrejs, ifrej))
            })
            numer1 <- as.logical(nrejs[2, ])            
            nrejs <- as.numeric(nrejs[1, ])
            numer2 <- head(knots, -1) >= abs(extraknot)
            list(knots = knots, nrejs = nrejs,
                 numer1 = numer1, numer2 = numer2)
        })

        if (res_init$safe){
            expt <- sapply(rejfun, function(fun){
                inds <- which(fun$numer1 | fun$numer2)
                sum((pt(fun$knots[inds + 1], df = df) - pt(fun$knots[inds], df = df)) / fun$nrejs[inds])
            })
            expt0 <- NA            
            expt <- sum(expt) * n
            ifrej <- expt <= alpha
            return(c(ifrej, expt0, expt))
        } else {
            expt <- sapply(rejfun, function(fun){
                inds <- which(fun$numer1)
                sum((pt(fun$knots[inds + 1], df = df) - pt(fun$knots[inds], df = df)) / fun$nrejs[inds])
            })
            expt0 <- sum(expt) * n            
            if (expt0 <= alpha){
                expt <- sapply(rejfun, function(fun){
                    inds <- which(fun$numer1 | fun$numer2)
                    sum((pt(fun$knots[inds + 1], df = df) - pt(fun$knots[inds], df = df)) / fun$nrejs[inds])
                })
                expt <- sum(expt) * n                
                ## This is the case where "or" is invoked.
                if (i %in% dBH_rej0){
                    ## If i is rejected in BH, it must be rejected here
                    ifrej <- TRUE                    
                } else {
                    ## If i is not rejected in BH, it is rejected if p_i <= \tau_i
                    ifrej <- expt <= alpha                                  }
                return(c(ifrej, expt0, expt))
            } else {
                expt <- sapply(rejfun, function(fun){
                    inds <- which(fun$numer1 & fun$numer2)
                    sum((pt(fun$knots[inds + 1], df = df) - pt(fun$knots[inds], df = df)) / fun$nrejs[inds])
                })
                expt <- sum(expt) * n                
                ## This is the case where "and" is invoked.
                if (!(i %in% dBH_rej0)){
                    ## If i is not rejected in BH, it must not be rejected here
                    ifrej <- FALSE
                } else {
                    ## If i is rejected in BH, it is rejected if p_i <= \tau_i
                    ifrej <- expt <= alpha
                }
                return(c(ifrej, expt0, expt))
            }
        }
    })

    ifrej <- as.logical(cand_info[1, ])
    rejlist <- which(ifrej == 1)
    rejlist <- c(init_rejlist, cand[rejlist])
    expt0 <- cand_info[2, ]
    expt <- cand_info[3, ]
    if (length(rejlist) == 0){
        return(list(rejs = numeric(0),
                    initrejs = numeric(0),
                    cand = cand,
                    expt = expt,
                    expt0 = expt0,
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
                    expt0 = expt0,
                    safe = res_init$safe,
                    secBH = FALSE))
    }

    u <- runif(Rplus)
    secBH_fac <- Rinit[rejlist] / Rplus    
    tdp <- u * secBH_fac
    nrejs <- nrejs_BH(tdp, 1)
    thr <- max(nrejs, 1) / Rplus
    secrejs <- which(tdp <= thr)
    rejs <- rejlist[secrejs]
    return(list(rejs = rejs,
                initrejs = rejlist,
                cand = cand,
                expt = expt,
                expt0 = expt0,
                safe = FALSE,
                secBH = TRUE,
                secBH_fac = secBH_fac))
}
