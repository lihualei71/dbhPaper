source("dBH_utils.R")

# Reconstruct the vector of z values after changing one of them
# Inputs:
#  zstat: old value of z[1]
#  newz:  new value of z[1]
#  s:     z[-1] - cor * z[1], which is the "residualized" version of z[-1] which is independent of z[1]
#  cor:   Corr(z[1], z[-1])
# Output: new value of the entire z vector holding s constant and updating z[1] <- newz
recover_stats_mvgauss <- function(zstat, newz, s, cor){
    c(newz, s + cor * newz)
}

# Compute a bounding box for (coef1 + coef2 * z: low <= z <= high)
thresh_bounds_mvgauss <- function(coef1, coef2, low, high){
    bound1 <- coef1 + coef2 * low
    bound2 <- coef1 + coef2 * high
    lower <- pmin(bound1, bound2)
    upper <- pmax(bound1, bound2)
    list(lower = lower, upper = upper)
}

# Return the values z for which a + b * z = thresh
# Inputs:
#   a:      intercept (scalar)
#   b:      slope (scalar)
#   thresh: n-vector of crossing points
# Outputs:
#   roots:  n-vector of values z at which a + b * z = thresh
#   sgn:    n-vector: sgn = 1 means upcrossing, sgn = -1 means downcrossing
#   posit:  n-vector: which of the n thresholds is crossed at each root
linroots_mvgauss <- function(a, b, thresh){
    roots <- (thresh - a) / b
    n <- length(thresh)
    sgn <- rep(sign(b), n)
    posit <- 1:n
    return(list(roots = roots, sgn = sgn, posit = posit))
}

# Find locations of z[1] at which RBH changes, holding the conditioning statistic fixed
# Inputs:
#   zstat:  z[1]  (scalar)
#   zminus: z[-1] (vector)
#   cor:    cor(z[1], z[-1]) in the population
#   thresh: list of BH z thresholds
#   low:    lowest value of z[1] (or abs(z[1])) to consider, typically obs value of z[1]
#   high:   largest value of z[1] we will consider, typically z(.01 * alpha/n)
#   side:   "one" (right-tailed) or "two" (two-tailed)
# Outputs:
#   res:    List of length one or two (number of tails). First element gives info about all
#           possible knots of the RBH function as z varies from low to high, second element
#           gives same for all possible knots as z varies from -high to -low.
compute_knots_mvgauss <- function(zstat, zminus, cor,
                                  alpha, side,
                                  low, high,
                                  avals, avals_type, beta
                                  ){
    n <- length(zminus) + 1
    navals <- length(avals)
    if (side == "two"){
        alpha <- alpha / 2
    }
    thresh <- qnorm(alpha * avals / n, lower.tail = FALSE)
    s <- zminus - cor * zstat  #Conditioning statistic S
    RCV <- list()
    if (side == "one"){
        xlow <- recover_stats_mvgauss(zstat, low, s, cor) # Full z vector at z[1] = low
        ## xlow[1] <- Inf
        ## Compute rejections as if p[1] = 0
        RCV[[1]] <- compute_RCV(xlow, thresh, avals)      # Initialize rejection-counting vector (rcv)
    } else if (side == "two"){
        xlow1 <- abs(recover_stats_mvgauss(zstat, low, s, cor))
        ## xlow1[1] <- Inf
        RCV[[1]] <- compute_RCV(xlow1, thresh, avals)     # Initialize rcv at left boundary of right tail
        xlow2 <- abs(recover_stats_mvgauss(zstat, -high, s, cor))
        ## xlow2[1] <- Inf
        RCV[[2]] <- compute_RCV(xlow2, thresh, avals)     #   Initialize rcv at left boundary of left tail
    }

    s <- c(0, s)
    cor <- c(1, cor)
    # Prepare to compute all locations where z[j] crosses a threshold as z[1] varies
    if (side == "one"){
        thrsgn <- rep(1, n)  # sign of the thresholds (+1 in R tail, -1 in L tail)
        coef1 <- s               # intercepts of z[i](z[1])
        coef2 <- cor             # slopes of z[i](z[1])
        hypid <- 1:n         # coordinate indices
    }
    if (side == "two"){
    # Same as above, but now tracking +z[i] and -z[i], for + and - values of z[1]
        thrsgn <- c(rep(1, 2 * n), rep(-1, 2 * n))
        coef1 <- c(s, -s, s, -s)
        coef2 <- c(cor, -cor, -cor, cor)
        hypid <- rep(1:n, 4)
    }

    # Compute bounding box of z[i] (and -z[i] for 2-sided), as z[1] varies
    #      between low and high (and -high and -low, for 2-sided)
    thr_bounds <- thresh_bounds_mvgauss(coef1, coef2, low, high)
    if (navals > 1){
        thrid_upper <- floor(pnorm(thr_bounds$lower, lower.tail = FALSE) * n / alpha - 1e-15)
        thrid_lower <- ceiling(pnorm(thr_bounds$upper, lower.tail = FALSE) * n / alpha + 1e-15)
        if (avals_type == "geom"){
            thrid_upper <- find_ind_geom_avals(beta, thrid_upper, "max")
            thrid_lower <- find_ind_geom_avals(beta, thrid_lower, "min")
        } else if (avals_type == "manual"){
            thrid_upper <- find_posit_vec(thrid_upper, avals, "left", FALSE)
            thrid_lower <- find_posit_vec(thrid_lower, avals, "right", FALSE)
        }
        rmids <- which(thrid_lower > navals |
                       thrid_upper < 1 |
                       thrid_upper < thrid_lower) # Which coordinates aren't removed
        ids <- (1:length(coef1))[-rmids]
    } else {
        ids <- (1:length(coef1))[thr_bounds$lower <= thresh & thr_bounds$upper >= thresh]
    }

    res <- list()
    if (side == "one"){
        ntails <- 1
        ids <- list(ids)
        tail_lbound <- low
        tail_ubound <- high
    } else if (side == "two"){
        ntails <- 2
        ids <- list(ids[thrsgn[ids] == 1],
                    ids[thrsgn[ids] == -1])
        tail_lbound <- c(low, -high)
        tail_ubound <- c(high, -low)
    }

    for (tail in 1:ntails){             # R tail first and then L tail if 2-sided
        if (length(ids[[tail]]) == 0){
            # Record that this tail has no knots
            res[[tail]] <- list(knots = numeric(0),
                                hyp = numeric(0),
                                posit = numeric(0),
                                sgn = numeric(0),
                                low = tail_lbound[tail],
                                high = tail_ubound[tail],
                                RCV = RCV[[tail]])
            next
        }

        knots <- list()
        hyp <- list()
        posit <- list()
        sgn <- list()

        # Iterate over z coordinates to collect info on all potential knots
        for (k in 1:length(ids[[tail]])){
            i <- ids[[tail]][k]
            if (navals > 1){
                thrids <- max(1, thrid_lower[i]):min(navals, thrid_upper[i])
                thr <- thresh[thrids]                # which thresholds this coordinate crosses
            } else {
                thrids <- 1
                thr <- thresh
            }
            sol <- linroots_mvgauss(coef1[i], coef2[i], thr)
            knots[[k]] <- sol$roots * thrsgn[i]  # locations of z[1] where threshold crossed
            nroots <- length(sol$roots)
            hyp[[k]] <- rep(hypid[i], nroots)    # which coordinate crossed
            posit[[k]] <- thrids[sol$posit] - 1  # which threshold crossed (R->C++ indexing)
            sgn[[k]] <- thrsgn[i] * sol$sgn      # up- or down-crossing (is that quite right?)
        }
        knots <- .Internal(unlist(knots, F, F))
        ord <- order(knots)             # This takes a long time too but is harder to get around
        hyp <- .Internal(unlist(hyp, F, F))[ord]
        posit <- .Internal(unlist(posit, F, F))[ord]
        sgn <- .Internal(unlist(sgn, F, F))[ord]
        knots <- knots[ord]

        res[[tail]] <- list(knots = knots,
                            hyp = hyp,
                            posit = posit,
                            sgn = sgn,
                            low = tail_lbound[tail],
                            high = tail_ubound[tail],
                            RCV = RCV[[tail]])
    }
    return(res)
}
