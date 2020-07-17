library("Rcpp")
sourceCpp("RBH_homotopy.cpp")
source("dBH_utils.R")
source("compute_knots_mvgauss.R")

dBH_mvgauss_qc <- function(zvals,
                           Sigma = NULL,
                           Sigmafun = NULL,
                           side = c("one", "two"),
                           alpha = 0.05, gamma = NULL,
                           is_safe = FALSE,
                           avals = NULL, 
                           avals_type = c("BH", "geom", "bonf", "manual"),
                           beta = 2,
                           eps = 0.05,
                           qcap = 2){
    n <- length(zvals)
    alpha0 <- gamma * alpha
    ntails <- ifelse(side == "two", 2, 1)    
    high <- qnorm(alpha * eps / n / ntails, lower.tail = FALSE)
    pvals <- zvals_pvals(zvals, side)
    qvals <- qvals_BH_reshape(pvals, avals)
    obj <- RBH_init(pvals, qvals, alpha, alpha0,
                    avals, is_safe, qcap)

    if (length(obj$cand) == 0){
        return(list(rejs = obj$init_rejlist,
                    initrejs = obj$init_rejlist,
                    cand = numeric(0),
                    expt = numeric(0),
                    safe = is_safe,
                    secBH = FALSE))
    }

    cand_info <- sapply(obj$cand, function(i){
        low <- qnorm(qvals[i] * max(avals) / n / ntails, lower.tail = FALSE)
        if (!is.null(Sigma)){
            cor <- Sigma[-i, i]
        } else {
            cor <- Sigmafun(i)[-i]
        }
        
        ## RBH function with alpha = qi
        res_q <- compute_knots_mvgauss(
            zstat = zvals[i],
            zminus = zvals[-i],
            cor = cor,
            alpha = qvals[i],
            side = side,            
            low = low,
            high = high,
            avals = avals,
            avals_type = avals_type,
            beta = beta)
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
                thra <- find_ind_geom_avals(beta, nrejs, "max")
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

        ## RBH function with alpha = alpha0
        res_alpha0 <- compute_knots_mvgauss(
            zstat = zvals[i],
            zminus = zvals[-i],
            cor = cor,
            alpha = alpha0,
            side = side,            
            low = low, 
            high = high,
            avals = avals,
            avals_type = avals_type,
            beta = beta)
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
                thra <- find_ind_geom_avals(beta, nrejs, "max")
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

        ## Combine two RBH functions
        res <- RBHfun_combine(res_q, res_alpha0)

        ## Compute conditional expectation
        expt <- sapply(res, function(re){
            compute_cond_exp(abs(zvals[i]), re$knots, re$nrejs, re$thr, dist = pnorm)
        })
        expt <- sum(expt) * n
        ifrej <- expt <= alpha
        return(c(ifrej, expt))        
    })

    ifrej <- as.logical(cand_info[1, ])
    rejlist <- which(ifrej)
    rejlist <- c(obj$init_rejlist, obj$cand[rejlist])
    expt <- cand_info[2, ]
    if (length(rejlist) == 0){
        return(list(rejs = numeric(0),
                    initrejs = numeric(0),
                    cand = obj$cand,
                    expt = expt,
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
                safe = FALSE,
                secBH = TRUE,
                secBH_fac = secBH_fac))
}
