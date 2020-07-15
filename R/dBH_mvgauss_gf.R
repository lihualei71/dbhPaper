library("Rcpp")
sourceCpp("RBH_homotopy.cpp")
source("dBH_utils.R")
source("compute_knots_mvgauss.R")

dBH_mvgauss_gf <- function(zvals, Sigma,
                           side = c("right", "left", "two"),
                           alpha = 0.05, alpha0 = NULL,
                           avals = NULL, 
                           avals_type = c("BH", "geom", "bonf", "manual"),
                           gamma = 2,
                           eps = 0.05,
                           use_cap = TRUE,
                           acc_multiplier = 2,
                           acc_minrejs = 10){
    side <- side[1]
    avals_type <- avals_type[1]
    n <- length(zvals)
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
        is_safe <- TRUE
    } else {
        is_safe <- (alpha0 <= alpha / normalize(avals))
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
    low <- qnorm(alpha0 * max(avals) / n / ntails, lower.tail = FALSE)
    high <- qnorm(alpha * eps / n / ntails, lower.tail = FALSE)
    pvals <- zvals_pvals(zvals, side)
    obj <- RBH_init(pvals, alpha, avals, alpha0, is_safe,
                    use_cap, acc_multiplier, acc_minrejs)

    if (length(obj$cand) == 0){
        return(list(rejs = obj$init_rejlist,
                    initrejs = obj$init_rejlist,
                    cand = numeric(0),
                    expt = numeric(0),
                    expt0 = numeric(0),
                    safe = is_safe,
                    secBH = FALSE))
    }

    cand_info <- sapply(obj$cand, function(i){
        res <- compute_knots_mvgauss(
            zstat = zvals[i],
            zminus = zvals[-i],
            cor = Sigma[-i, i],
            alpha = alpha0,
            side = side,            
            low = min(low, abs(zvals[i])),
            high = high,
            avals = avals,
            avals_type = avals_type,
            gamma = gamma)
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
            list(knots = knots, nrejs = nrejs, thr = thr)
        })
        zstat <- abs(zvals[i])
        if (is_safe){
            expt <- sapply(res, function(re){
                thr <- pmin(zstat, re$thr)
                compute_cond_exp(zstat, re$knots, re$nrejs, thr, dist = pnorm)
            })
            expt0 <- NA            
            expt <- sum(expt) * n
            ifrej <- expt <= alpha
            return(c(ifrej, expt0, expt))
        } else {
            expt <- sapply(res, function(re){
                compute_cond_exp(zstat, re$knots, re$nrejs, re$thr, dist = pnorm)
            })
            expt0 <- sum(expt) * n
            if (expt0 <= alpha){
                expt <- sapply(res, function(re){
                    thr <- pmin(zstat, re$thr)
                    compute_cond_exp(zstat, re$knots, re$nrejs, thr, dist = pnorm)
                })
                expt <- sum(expt) * n
                ## This is the case where "or" is invoked.
                if (i %in% obj$dBH_rej0){
                    ## If i is rejected in BH, it must be rejected here
                    ifrej <- TRUE
                } else {
                    ## If i is not rejected in BH, it is rejected if p_i <= \tau_i
                    ifrej <- expt <= alpha
                }
                return(c(ifrej, expt0, expt))
            } else {
                expt <- sapply(res, function(re){
                    thr <- pmax(zstat, re$thr)
                    compute_cond_exp(zstat, re$knots, re$nrejs, thr, dist = pnorm)
                })
                expt <- sum(expt) * n
                ## This is the case where "and" is invoked.
                if (!(i %in% obj$dBH_rej0)){
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
    rejlist <- which(ifrej)
    rejlist <- c(obj$init_rejlist, obj$cand[rejlist])
    expt0 <- cand_info[2, ]
    expt <- cand_info[3, ]
    if (length(rejlist) == 0){
        return(list(rejs = numeric(0),
                    initrejs = numeric(0),
                    cand = obj$cand,
                    expt = expt,
                    expt0 = expt0,
                    safe = is_safe,
                    secBH = FALSE))
    }

    rejlist <- sort(rejlist)
    Rplus <- length(rejlist)
    if (Rplus >= max(obj$Rinit[rejlist])){
        return(list(rejs = rejlist,
                    initrejs = rejlist,
                    cand = obj$cand,
                    expt = expt,
                    expt0 = expt0,
                    safe = is_safe,
                    secBH = FALSE))
    }

    uvec <- runif(Rplus)
    secBH_fac <- obj$Rinit[rejlist] / Rplus
    tdp <- uvec * secBH_fac
    nrejs <- nrejs_BH(tdp, 1)
    thr <- max(nrejs, 1) / Rplus
    secrejs <- which(tdp <= thr)
    rejs <- rejlist[secrejs]
    return(list(rejs = rejs,
                initrejs = rejlist,
                cand = obj$cand,
                expt = expt,
                expt0 = expt0,
                safe = FALSE,
                secBH = TRUE,
                secBH_fac = secBH_fac))
}
