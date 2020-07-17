zvals_pvals <- function(zvals, side){
    if (side == "one"){
        pnorm(zvals, lower.tail = FALSE)
    } else if (side == "two"){
        2 * pnorm(abs(zvals), lower.tail = FALSE)
    }
}

tvals_pvals <- function(tvals, df, side){
    if (side == "one"){
        pt(tvals, df = df, lower.tail = FALSE)
    } else if (side == "two"){
        2 * pt(abs(tvals), df = df, lower.tail = FALSE)
    }        
}

nrejs_BH_reshape <- function(pvals, alpha, avals = 1:length(pvals)){
    n <- length(pvals)
    if (identical(avals, 1:n)){
        fac <- 1:n
    } else {
        fac <- fill_int(n, avals)
    }
    adjust_pvals <- sort(pvals) * n / alpha / fac
    nrejs <- max(0, which(adjust_pvals <= 1))
    if (nrejs == 0){
        aid <- 0
    } else {
        aid <- fac[nrejs]
    }
    return(list(nrejs = nrejs, aid = aid))
}

qvals_BH_reshape <- function(pvals, avals = 1:length(pvals)){
    n <- length(pvals)
    if (identical(avals, 1:n)){
        fac <- 1:n
    } else {
        fac <- fill_int(n, avals)
    }
    ord <- order(pvals)
    adjust_pvals <- cummin(rev(pvals[ord] * n / fac))
    qvals <- rep(NA, n)
    qvals[ord] <- rev(adjust_pvals)
    return(qvals)
}

BH <- function(pvals, alpha = 0.05, avals = 1:length(pvals),
               reshape = TRUE){
    n <- length(pvals)
    if (reshape){
        alpha <- alpha / normalize(avals)
    }
    tmp <- nrejs_BH_reshape(pvals, alpha, avals)
    rejs <- which(pvals <= tmp$aid * alpha / n)
    list(nrejs = tmp$nrejs, rejs = rejs)
}

# Compute the rejection-corrected vector given thresholds and a-values. Equivalent to compute_mod_nrejs when avals = 1:length(thresh)
# NOTE: Assumes thresh is in decreasing order
compute_RCV <- function(x, thresh, avals){
    ord <- order(c(x, thresh), decreasing = TRUE)
    ones.zeros <- c(rep(1, length(x)), rep(0, length(thresh)))[ord]
    cumsum(ones.zeros)[ones.zeros == 0] - avals
}

compute_cond_exp <- function(stat, knots, nrejs, thr, dist){
    right_knots <- tail(knots, -1)
    if (length(thr) != length(right_knots)){
        browser()
    }
    inds <- which(thr <= right_knots)
    left_knots <- knots[inds]
    left_knots <- pmax(thr[inds], left_knots)    
    tmp <- (dist(right_knots[inds]) - dist(left_knots)) / nrejs[inds]
    sum(tmp)
}

## Return a piecewise constant vector x with
## x[i] = max{vals[j]: vals[j] <= i}
fill_int <- function(n, vals){
    vec <- 1:n
    vec[-vals] <- 0
    cummax(vec)
}

## Need to deal with ties, in both vec and target
find_posit_vec <- function(vec, target, dir, decreasing = TRUE){
    if (decreasing){
        posit <- rank(-c(vec, target))[1:length(vec)] - rank(-vec)
    } else {
        posit <- rank(c(vec, target))[1:length(vec)] - rank(vec)
    }
    if (dir == "right"){
        posit <- posit + 1
    }
    return(posit)
}

## Return a piecewise constant vector y with
## y[i] = max{vals[j]: vals[j] <= x[i]}
fill_int_general <- function(x, vals){
    inds <- find_posit_vec(x, vals - 1e-6, "left", FALSE)
    vals[inds]
}

## a_i = ceiling((gamma^(i-1) - 1) / (gamma - 1) + 1)
## Find max{i: a_i <= target} or min{i: a_i >= target}
find_ind_geom_avals <- function(gamma, target, type){
    if (type == "max"){
        inds <- rep(0, length(target))        
        inds[target >= 1] <- floor(log((gamma - 1) * (target[target >= 1] - 1) + 1, gamma)) + 1
    } else if (type == "min"){
        inds <- rep(1, length(target))        
        inds[target >= 2] <- ceiling(log((gamma - 1) * (target[target >= 2] - 2) + 1 + 1e-10, gamma)) + 1
        inds[target < 1] <- 0
    }
    return(inds)
}

## Generate a_i = ceiling((gamma^(i-1) - 1) / (gamma - 1) + 1) for all a_i <= n
geom_avals <- function(gamma, n){
    m <- find_ind_geom_avals(gamma, n, "max")
    ceiling((gamma^(0:(m-1)) - 1) / (gamma - 1) + 1)
}

nrejs_BH <- function(pvals, alpha){
    n <- length(pvals)
    adjust_pvals <- sort(pvals) * n / alpha / 1:n
    max(0, which(adjust_pvals <= 1))
}

RBH_init <- function(pvals, qvals, alpha, alpha0,
                     avals, is_safe, qcap){
    n <- length(pvals)
    dBH_rej0 <- BH(pvals, alpha0, avals, FALSE)$rejs
    Rinit <- rep(length(dBH_rej0) + 1, n)
    Rinit[dBH_rej0] <- Rinit[dBH_rej0] - 1

    init_rejlist <- which(qvals <= alpha / max(avals))    
    if (is_safe){
        init_rejlist <- union(dBH_rej0, init_rejlist)
    }
    init_acclist <- which(qvals >= qcap * alpha)
    cand <- 1:n
    tmp <- c(init_rejlist, init_acclist)
    if (length(tmp) > 0){
        cand <- cand[-tmp]
    }
    
    return(list(Rinit = Rinit, cand = cand,
                dBH_rej0 = dBH_rej0,
                init_rejlist = init_rejlist,
                init_acclist = init_acclist))
}

expand_piecewise_const <- function(x, y, extrax){
    high <- tail(x, 1)
    x <- c(head(x, -1), extrax)
    ord <- order(x)
    tmp <- c(rep(1, length(y)), rep(0, length(extrax)))[ord]
    x <- c(x[ord], high)
    list(x = x, y = y[cumsum(tmp)])
}

RBHfun_combine <- function(res_q, res_alpha0){
    ntails <- length(res_q)
    lapply(1:ntails, function(k){
        knots_q <- res_q[[k]]$knots
        knots_alpha0 <- res_alpha0[[k]]$knots
        thr <- expand_piecewise_const(
            knots_q, res_q[[k]]$thr,
            knots_alpha0[-c(1, length(knots_alpha0))])
        nrejs <- expand_piecewise_const(
            knots_alpha0, res_alpha0[[k]]$nrejs,
            knots_q[-c(1, length(knots_q))])
        list(knots = nrejs$x, nrejs = nrejs$y, thr = thr$y)
    })
}

normalize <- function(avals){
    length(avals) - sum(head(avals, -1) / tail(avals, -1))
}

find_int_above_thr <- function(knots, thr){
    ## To avoid numerical errors when a knot is equal to thr
    midknots <- (head(knots, -1) + tail(knots, -1)) / 2
    tmp <- rle(midknots >= thr)
    ids <- cumsum(tmp$lengths)
    start <- c(1, head(ids, -1) + 1)
    end <- ids
    lapply(which(tmp$values), function(i){
        c(knots[start[i]], knots[end[i] + 1])
    })
}

lingrid <- function(low, high, gridsize, side){
    grid <- list()
    grid[[1]] <- seq(low, high, length.out = gridsize)
    if (side == "two"){
        grid[[2]] <- seq(-high, -low, length.out = gridsize)
    }
    return(grid)
}

lm_mvt <- function(y, X, subset){
    n <- nrow(X)
    p <- ncol(X)
    df <- n - p
    fit <- lm(y ~ X + 0)
    zvals <- as.numeric(coefficients(fit))
    Sigma <- solve(t(X) %*% X)
    sigmahat <- summary(fit)$sigma
    tvals <- zvals / sigmahat
    list(tvals = tvals[subset], df = df,
         Sigma = Sigma[subset, subset])
}

