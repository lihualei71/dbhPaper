source("dBH_utils.R")

# Reconstruct the vector of t values after changing one of them
# Inputs:
#  tstat: old value of t[1]
#  newt:  new value of t[1]
#  s:     t[-1] - cor * t[1]
#  cor:   Corr(z[1], z[-1])
#  df:    residual degrees of freedom
# Output: new value of the entire z vector holding s constant and updating z[1] <- newz
recover_stats_mvt <- function(tstat, newt, s, cor, df){
    c(newt, s * sqrt((df + newt^2) / (df + tstat^2)) +
          cor * newt)
}

# Compute a bounding box for the function
#      f(x; a, b) = a * sqrt(1 + x^2) + b * x
# Inputs:
#   coef1: vector of a coefficients
#   coef2: vector of b coefficients
#   low:   lower bound of x
#   high:  upper bound of x
thresh_bounds_mvt <- function(coef1, coef2, low, high){
    n <- length(coef1)
    bound1 <- coef1 * sqrt(1 + low^2) + coef2 * low
    bound2 <- coef1 * sqrt(1 + high^2) + coef2 * high
    bound3 <- rep(NA, n)
    inds <- (abs(coef1) > abs(coef2)) & (coef1 * coef2 < 0)
    coef1 <- coef1[inds]
    coef2 <- coef2[inds]
    bound3[inds] <- sign(coef1) * sqrt(coef1^2 - coef2^2)
    lower <- pmin(bound1, bound2, bound3, na.rm = TRUE)
    upper <- pmax(bound1, bound2, bound3, na.rm = TRUE)
    list(lower = lower, upper = upper)
}

# Return the values x for which f(x; a, b) = thresh, where
#     f(x; a, b) = a * sqrt(1 + x^2) + b * x
#     gives t[i] / sqrt(df) when t[1] = x * sqrt(df)
# Inputs:
#   a:      coefficient for f
#   b:      coefficient for f
#   thresh: n-vector of crossing points
# Outputs:
#   roots:  n-vector of values t at which f(t; a, b) = thresh
#   sgn:    n-vector: sgn = 1 means upcrossing, sgn = -1 means downcrossing (??is this right??)
#   posit:  n-vector: which of the n thresholds is crossed at each root
quadroots_mvt <- function(a, b, thresh){
    n <- length(thresh)    
    if (b == 0){
        roots <- sqrt(thresh^2 / a^2 - 1)
        roots <- c(roots, -roots)
        sgn <- c(rep(1, n), rep(-1, n))
        posit <- c(1:n, 1:n)
        return(list(roots = roots, sgn = sgn, posit = posit))        
    }
    numer1 <- b * thresh
    numer2 <- a * sqrt(b^2 + thresh^2 - a^2)
    denom <- b^2 - a^2
    if (a == 0){
        roots <- thresh / b
        sgn <- rep(sign(b), n)
        posit <- 1:n
    } else if (denom == 0){
        roots <- (thresh^2 - b^2) / 2 / b / thresh
        sgn <- rep(sign(b), n)
        posit <- 1:n
    } else if (a < abs(b) || (a > b && b > 0)){
        roots <- (numer1 - numer2 * sign(b)) / denom
        sgn <- rep(sign(b), n)
        posit <- 1:n
    } else if (a > -b && b < 0){
        roots <- (numer1 - numer2) / denom
        inds <- which(thresh <= a)
        sgn <- rep(1, n)
        posit <- 1:n
        m <- length(inds)
        if (m > 0){
            roots <- c(roots, (numer1[inds] + numer2[inds]) / denom)
            sgn <- c(sgn, rep(-1, m))
            posit <- c(posit, inds)
        }
    }
    return(list(roots = roots, sgn = sgn, posit = posit))
}


# Find locations of t[1] at which RBH changes, holding the conditioning statistic fixed
# Inputs:
#   zstat:  t[1]  (scalar)
#   zminus: t[-1] (vector)
#   cor:    cor(z[1], z[-1]) in the population
#   thresh: list of BH t thresholds
#   low:    lowest value of t[1] (or abs(t[1])) to consider, typically obs value of t[1]
#   high:   largest value of t[1] we will consider, typically t(.01 * alpha/n)
#   side:   "one" (right-tailed) or "two" (two-tailed)
# Outputs:
#   res:    List of length one or two (number of tails). First element gives info about all
#           possible knots of the RBH function as t varies from low to high, second element
#           gives same for all possible knots as t varies from -high to -low.
compute_knots_mvt <- function(tstat, tminus, df, cor,
                              alpha, side,
                              low, high,
                              avals, avals_type, gamma
                              ){
    n <- length(tminus) + 1
    navals <- length(avals)    
    if (side == "two"){
        alpha <- alpha / 2
    }
    thresh <- qt(alpha * avals / n, df = df, lower.tail = FALSE)
    s <- tminus - cor * tstat  #Conditioning statistic S
    RCV <- list()
    if (side == "one"){
        xlow <- recover_stats_mvt(tstat, low, s, cor, df)  # Full t vector at t[1] = low
        ## xlow[1] <- Inf
        ## Compute rejections as if p[1] = 0
        RCV[[1]] <- compute_RCV(xlow, thresh, avals)       # Initialize rejection-counting vector (rcv)
    } else if (side == "two"){
        xlow1 <- abs(recover_stats_mvt(tstat, low, s, cor, df))
        ## xlow1[1] <- Inf
        RCV[[1]] <- compute_RCV(xlow1, thresh, avals)              # Initialize rcv at left boundary of right tail
        xlow2 <- abs(recover_stats_mvt(tstat, -high, s, cor, df))
        ## xlow2[1] <- Inf
        RCV[[2]] <- compute_RCV(xlow2, thresh, avals)              # Initialize rcv at left boundary of left tail
    }

    # Rescale equations to be of the form a * sqrt(1 + x^2) + b * x = c,
    #        by dividing out sqrt(df)
    sqdf <- sqrt(df)
    thresh <- thresh / sqdf
    tstat <- tstat / sqdf
    tminus <- tminus / sqdf
    low <- low / sqdf
    high <- high / sqdf
    s <- c(0, s)
    cor <- c(1, cor)

    # Prepare to compute all locations where t[i] crosses a threshold as t[1] varies
    if (side == "one"){
        thrsgn <- rep(1, n)
        coef1 <- s / sqrt(1 + tstat^2) / sqdf
        coef2 <- cor
        hypid <- 1:n
    }
    if (side == "two"){
        thrsgn <- c(rep(1, 2 * n), rep(-1, 2 * n))
        tmp <- s / sqrt(1 + tstat^2) / sqdf
        coef1 <- c(tmp, -tmp, tmp, -tmp)
        coef2 <- c(cor, -cor, -cor, cor)
        hypid <- rep(1:n, 4)
    }

    # Screen out coordinates that will never cross any threshold
    ids <- which(coef1 > 0 | coef2 > abs(coef1))

    coef1 <- coef1[ids]
    coef2 <- coef2[ids]
    hypid <- hypid[ids]
    thrsgn <- thrsgn[ids]

    # Compute bounding box of z[i] (and -z[i] for 2-sided), as z[1] varies
    #      between low and high (and -high and -low, for 2-sided)    
    thr_bounds <- thresh_bounds_mvt(coef1, coef2, low, high)
    if (navals > 1){
        thrid_upper <- floor(pt(thr_bounds$lower * sqdf, df = df, lower.tail = FALSE) * n / alpha - 1e-15)
        thrid_lower <- ceiling(pt(thr_bounds$upper * sqdf, df = df, lower.tail = FALSE) * n / alpha + 1e-15)
        if (avals_type == "geom"){
            thrid_upper <- find_ind_geom_avals(gamma, thrid_upper, "max")
            thrid_lower <- find_ind_geom_avals(gamma, thrid_lower, "min")
        } else if (avals_type == "manual"){
            thrid_upper <- find_posit_vec(thrid_upper, avals, "left", FALSE)
            thrid_lower <- find_posit_vec(thrid_lower, avals, "right", FALSE)
        }
        rmids <- which(thrid_lower > navals | thrid_upper < 1 | thrid_upper < thrid_lower) # Which coordinates aren't removed
        ids <- (1:length(coef1))[-rmids]
    } else {
        ids <- (1:length(ids))[thr_bounds$lower <= thresh & thr_bounds$upper >= thresh]
    }

    res <- list()
    if (side == "one"){
        ntails <- 1
        ids <- list(ids)
        lowlist <- low
        highlist <- high
    } else if (side == "two"){
        ntails <- 2
        ids <- list(ids[thrsgn[ids] == 1],
                    ids[thrsgn[ids] == -1])
        lowlist <- c(low, -high)
        highlist <- c(high, -low)
    }

    for (tail in 1:ntails){
        if (length(ids[[tail]]) == 0){
            res[[tail]] <- list(knots = numeric(0),
                             hyp = numeric(0),
                             posit = numeric(0),
                             sgn = numeric(0),
                             low = lowlist[tail] * sqdf,
                             high = highlist[tail] * sqdf,
                             RCV = RCV[[tail]])
            next
        }

        knots <- list()
        hyp <- list()
        posit <- list()
        sgn <- list()

        ## Iterate over t coordinates to collect info on all potential knots
        for (k in 1:length(ids[[tail]])){
            i <- ids[[tail]][k]
            if (navals > 1){
                thrids <- max(1, thrid_lower[i]):min(navals, thrid_upper[i])
                thr <- thresh[thrids]                # which thresholds this coordinate crosses
            } else {
                thrids <- 1
                thr <- thresh
            }
            sol <- quadroots_mvt(coef1[i], coef2[i], thr)
            inds <- sol$roots <= high & sol$roots >= low
            m <- sum(inds)
            if (m == 0){
                next
            }
            knots[[k]] <- sol$roots[inds] * thrsgn[i]
            hyp[[k]] <- rep(hypid[i], m)
            posit[[k]] <- thrids[sol$posit[inds]] - 1
            sgn[[k]] <- thrsgn[i] * sol$sgn[inds]
        }
        knots <- .Internal(unlist(knots, F, F))
        if (length(knots) > 0){
            ord <- order(knots)
            hyp <- .Internal(unlist(hyp, F, F))[ord]
            posit <- .Internal(unlist(posit, F, F))[ord]
            sgn <- .Internal(unlist(sgn, F, F))[ord]
            knots <- knots[ord]
        } else {
            hyp <- posit <- sgn <- numeric(0)
        }

        res[[tail]] <- list(knots = knots * sqdf,
                            hyp = hyp,
                            posit = posit,
                            sgn = sgn,
                            low = lowlist[tail] * sqdf,
                            high = highlist[tail] * sqdf,
                            RCV = RCV[[tail]])
    }
    return(res)
}
