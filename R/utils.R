source("dBH_utils.R")

pvals_mvgauss <- function(zvals, Sigma, side){
    zvals <- zvals / sqrt(diag(Sigma))
    if (side == "right"){
        side <- "one"
    } else if (side == "left"){
        side <- "one"
        zvals <- -zvals
    }
    zvals_pvals(zvals, side)
}

pvals_mvt <- function(zvals, Sigma, sigmahat, df, side){
    zvals <- zvals / sqrt(diag(Sigma))
    tvals <- zvals / sigmahat
    if (side == "right"){
        side <- "one"
    } else if (side == "left"){
        side <- "one"
        zvals <- -zvals
    }
    tvals_pvals(tvals, df, side)
}

BC <- function(pvals, alpha = 0.05){
    n <- length(pvals)
    sorted_mask_pvals <- sort(pmin(pvals, 1 - pvals))
    fdphat <- sapply(sorted_mask_pvals, function(thresh){
        (1 + sum(pvals >= 1 - thresh))/max(1, sum(pvals <= thresh))
    })
    khat <- which(fdphat <= alpha)
    if (length(khat) == 0){
        return(list(nrejs = 0, rejs = c()))
    } else {
        khat <- max(khat)
        phat <- sorted_mask_pvals[khat]
        rejs <- which(pvals <= phat)
        return(list(nrejs = length(rejs), rejs = rejs))
    }
}

knockoff_lm <- function(X, Xk, y, alpha = 0.05){
    ## W <- knockoff::stat.glmnet_coefdiff(X, Xk, y, nfolds=10)
    W <- knockoff::stat.glmnet_lambdasmax(X, Xk, y)
    thr <- knockoff::knockoff.threshold(W, fdr = alpha, offset = 1)
    rejs <- which(W >= thr)
    list(rejs = rejs, nrejs = length(rejs))
}

mineig <- function(A){
    min(eigen(A, symmetric = TRUE)$values)
}

solve_sdp <- function(Sigma, gaptol = 1e-06,
                      maxit = 1000, psdtol = 1e-09){
    stopifnot(isSymmetric(Sigma))
    G <- stats::cov2cor(Sigma)
    p <- dim(G)[1]
    if (mineig(G) < psdtol) {
        stop("The covariance matrix is not positive-definite: cannot solve SDP", 
            immediate. = T)
    }
    Cl1 <- rep(0, p)
    Al1 <- -Matrix::Diagonal(p)
    Cl2 <- rep(1, p)
    Al2 <- Matrix::Diagonal(p)
    d_As <- c(diag(p))
    As <- Matrix::Diagonal(length(d_As), x = d_As)
    As <- As[which(Matrix::rowSums(As) > 0), ]
    Cs <- c(G)
    A <- cbind(Al1, Al2, As)
    C <- matrix(c(Cl1, Cl2, Cs), 1)
    K <- NULL
    K$s <- p
    K$l <- 2 * p
    b <- rep(1, p)
    OPTIONS <- NULL
    OPTIONS$gaptol <- gaptol
    OPTIONS$maxit <- maxit
    OPTIONS$logsummary <- 0
    OPTIONS$outputstats <- 0
    OPTIONS$print <- 0
    sol <- Rdsdp::dsdp(A, b, C, K, OPTIONS)
    if (!identical(sol$STATS$stype, "PDFeasible")) {
        warning("The SDP solver returned a non-feasible solution")
    }
    s <- sol$y
    s[s < 0] <- 0
    s[s > 1] <- 1
    psd <- 0
    s_eps <- 1e-08
    while (mineig(G - diag(s * (1 - s_eps), length(s))) < psdtol) {
            s_eps <- s_eps * 10
    }
    s <- s * (1 - s_eps)
    if (max(s) == 0) {
        warning("In creation of SDP knockoffs, procedure failed. Knockoffs will have no power.", 
            immediate. = T)
    }
    return(s * diag(Sigma))
}

FDPpower <- function(rejs, H0){
    if (length(rejs) == 0){
        return(c(0, 0))
    }
    nrejs <- length(rejs)
    fp <- length(intersect(rejs, which(H0)))
    tp <- nrejs - fp
    FDP <- fp / nrejs
    power <- tp / sum(!H0)
    return(c(FDP, power))
}

#### The following four funtions compute the intersection of two lines given by (x1, y1) and (x2, y2) where x1 and y1 (resp. x2 and y2) are vectors of same length but length(x1) may differ from length(x2). The lines are interpolated linearly within the range and with the leftmost/rightmost value outside the range. 
interpolate_one_point <- function(x1, x2, y1, y2, newx){
    lambda <- (x2 - newx) / (x2 - x1)
    lambda * y1 + (1 - lambda) * y2
}

interpolate_two_lines <- function(x1, x2, y1, y2){
    n1 <- length(x1)
    n2 <- length(x2)
    n <- n1 + n2
    x <- c(x1, x2)
    y1left <- y1[1]
    y1right <- y1[n1]
    y2left <- y2[1]
    y2right <- y2[n2]

    y1 <- c(y1, rep(NA, n2))
    y2 <- c(rep(NA, n1), y2)
    ord <- order(x)
    y1 <- y1[ord]
    y2 <- y2[ord]
    x <- x[ord]

    if (is.na(y1[1])){
        y1[1] <- y1left
    }
    if (is.na(y1[n])){
        y1[n] <- y1right
    }
    if (is.na(y2[1])){
        y2[1] <- y2left
    }
    if (is.na(y2[n])){
        y2[n] <- y2right
    }

    for (i in 2:(n-1)){
        if (is.na(y1[i])){
            ind_left <- max(which(!is.na(y1[1:(i-1)])))
            ind_right <- i + min(which(!is.na(y1[(i+1):n])))
            y1[i] <- interpolate_one_point(x[ind_left], x[ind_right], y1[ind_left], y1[ind_right], x[i])
        }
        if (is.na(y2[i])){
            ind_left <- max(which(!is.na(y2[1:(i-1)])))
            ind_right <- i + min(which(!is.na(y2[(i+1):n])))
            y2[i] <- interpolate_one_point(x[ind_left], x[ind_right], y2[ind_left], y2[ind_right], x[i])
        }
    }

    return(list(x = x, y1 = y1, y2 = y2))
}

cross_point <- function(x1, x2, y11, y12, y21, y22){
    tmp <- (y11 - y21) / (y22 + y11 - y21 - y12)
    x <- x1 + (x2 - x1) * tmp
    y <- y11 + (y12 - y11) * tmp
    return(c(x, y))
}

find_cross_two_lines <- function(x1, x2, y1, y2){
    obj <- interpolate_two_lines(x1, x2, y1, y2)
    x <- obj$x
    y1 <- obj$y1
    y2 <- obj$y2
    n <- length(x)
    cross <- list()
    k <- 1
    for (i in 1:(n - 1)){
        if ((y1[i + 1] - y2[i + 1]) * (y1[i] - y2[i]) <= 0){
            cross[[k]] <- cross_point(x[i], x[i + 1], y1[i], y1[i + 1], y2[i], y2[i + 1])
            k <- k + 1
        }
    }
    return(unique(cross))
}
