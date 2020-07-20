create_fixed <- function(X, method = c("sdp", "equi"), sigma = NULL, y = NULL, randomize = F, intercept = T){
    method = match.arg(method)
    n = nrow(X)
    p = ncol(X)
    if (n <= p) 
        stop("Input X must have dimensions n > p")
    else if (n < 2 * p) {
        warning("Input X has dimensions p < n < 2p. ", "Augmenting the model with extra rows.", 
            immediate. = T)
        X.svd = svd(X, nu = n, nv = 0)
        u2 = X.svd$u[, (p + 1):n]
        X = rbind(X, matrix(0, 2 * p - n, p))
        if (is.null(sigma)) {
            if (is.null(y)) {
                stop("Either the noise level \"sigma\" or the response variables \"y\" must\n             be provided in order to augment the data with extra rows.")
            }
            else {
                sigma = sqrt(mean((t(u2) %*% y)^2))
            }
        }
        if (randomize) 
            y.extra = rnorm(2 * p - n, sd = sigma)
        else y.extra = with_seed(0, rnorm(2 * p - n, sd = sigma))
        y = c(y, y.extra)
    }
    X = knockoff:::normc(X, center = F)
    Xk = switch(match.arg(method), equi = create_equicorrelated(X, randomize, intercept), sdp = create_sdp(X, randomize, intercept))
    structure(list(X = X, Xk = Xk, y = y), class = "knockoff.variables")
}

create_equicorrelated <- function(X, randomize, intercept){
    X.svd = decompose(X, randomize, intercept)
    if (any(X.svd$d <= 1e-05 * max(X.svd$d))) 
        stop(paste("Data matrix is rank deficient.", "Equicorrelated knockoffs will have no power."))
    lambda_min = min(X.svd$d)^2
    s = min(2 * lambda_min, 1)
    s_diff = pmax(0, 2 * s - (s/X.svd$d)^2)
    X_ko = (X.svd$u %*diag% (X.svd$d - s/X.svd$d) + X.svd$u_perp %*diag% 
            sqrt(s_diff)) %*% t(X.svd$v)
    return(X_ko)
}

create_sdp <- function(X, randomize, intercept){
    X.svd = decompose(X, randomize, intercept)
    tol = 1e-05
    d = X.svd$d
    d_inv = 1/d
    d_zeros = d <= tol * max(d)
    if (any(d_zeros)) {
        warning(paste("Data matrix is rank deficient.", "Model is not identifiable, but proceeding with SDP knockoffs"), 
            immediate. = T)
        d_inv[d_zeros] = 0
    }
    G = (X.svd$v %*diag% d^2) %*% t(X.svd$v)
    G_inv = (X.svd$v %*diag% d_inv^2) %*% t(X.svd$v)
    s = knockoff:::create.solve_sdp(G)
    s[s <= tol] = 0
    C.svd = knockoff:::canonical_svd(2 * diag(s) - (s %diag*% G_inv %*diag% 
        s))
    X_ko = X - (X %*% G_inv %*diag% s) + (X.svd$u_perp %*diag% 
                                          sqrt(pmax(0, C.svd$d))) %*% t(C.svd$v)
    return(X_ko)    
}

decompose <- function(X, randomize, intercept){
    n = nrow(X)
    p = ncol(X)
    stopifnot(n >= 2 * p)
    result = knockoff:::canonical_svd(X)
    tmp = cbind(result$u, matrix(0, n, p))
    ids <- (p + 1):(2 * p)
    if (intercept){
        tmp = cbind(tmp, rep(1, n))
        ids <- ids + 1
    }
    Q = qr.Q(qr(tmp))
    u_perp = Q[, ids]
    if (randomize) {
        Q = qr.Q(qr(knockoff:::rnorm_matrix(p0, p0)))
        u_perp = u_perp %*% Q
    }
    result$u_perp = u_perp
    result
}

`%*diag%` <- function(X, d) t(t(X) * d)

`%diag*%` <- function(d, X) d * X
