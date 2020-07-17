source("dBH_mvgauss_qc.R")
source("dBH_mvgauss_qc_grid.R")

dBH_mvgauss <- function(zvals,
                        Sigma = NULL,
                        Sigmafun = NULL,
                        vars = NULL,
                        side = c("right", "left", "two"),
                        alpha = 0.05, gamma = NULL,
                        niter = 1,
                        tautype = "QC",
                        avals = NULL, 
                        avals_type = c("BH", "geom", "bonf", "manual"),
                        geom_fac = 2,
                        eps = 0.05,
                        qcap = 2,
                        gridsize = 20,
                        exptcap = 0.9,
                        is_safe = NULL){
    if (niter > 2){
        stop("\'niter\' can only be 1 or 2.")
    }

    side <- side[1]
    avals_type <- avals_type[1]
    tautype <- tautype[1]    
    n <- length(zvals)
    if (is.null(avals)){
        if (avals_type == "manual"){
            stop("avals must be inputted when avals_type = \"manual\"")
        } else if (avals_type == "geom" && geom_fac <= 1){
            stop("geom_fac must be larger than 1 when avals_type = \"geom\"")
        }
        avals <- switch(avals_type,
                        BH = 1:n,
                        geom = geom_avals(geom_fac, n),
                        bonf = 1)
    } else if (is.null(avals_type)) {
        if (avals[1] != 1){
            stop("The first element of avals must be 1.")
        }
        avals_type <- "manual"
        warning("avals is inputted and avals_type is set to be \"manual\" by default. This may slow down the code. Use the built-in avals_type (\"BH\", \"geom\" or \"bonf\") instead unless there is a good reason to use the inputted avals.")
    } else {
        stop("Set avals = NULL when avals_type is specified")
    }

    if (is.null(gamma)){
        gamma <- 1 / normalize(avals)
        is_safe <- TRUE
    } else if (is.null(is_safe)){
        is_safe <- (gamma <= 1 / normalize(avals))
    } else if (!is_safe){
        warning("Set is_safe = TRUE only if you can prove the procedure is safe (e.g. CPRDS case with gamma = 1)")
    }

    if (!is.null(Sigma) && any(diag(Sigma) != 1)){
        vars <- diag(Sigma)
        Sigma <- cov2cor(Sigma)        
    } else {
        if (is.null(vars)){
            vars <- rep(1, n)
        }
    }
    zvals <- zvals / sqrt(vars)    
    if (side == "left"){
        zvals <- -zvals
        side <- "one"
    } else if (side == "right"){
        side <- "one"
    }
    
    if (niter == 1){
        if (tautype == "QC"){
            dBH_mvgauss_qc(zvals = zvals,
                           Sigma = Sigma,
                           Sigmafun = Sigmafun,
                           side = side,
                           alpha = alpha,
                           gamma = gamma,
                           is_safe = is_safe,
                           avals = avals,
                           avals_type = avals_type,
                           geom_fac = geom_fac,
                           eps = eps,
                           qcap = qcap)
        }
    } else if (niter == 2){
        if (tautype == "QC"){
            dBH_mvgauss_qc_grid(zvals = zvals,
                                Sigma = Sigma,
                                Sigmafun = Sigmafun,
                                side = side,
                                alpha = alpha,
                                gamma = gamma,
                                is_safe = is_safe,
                                avals = avals,
                                avals_type = avals_type,
                                geom_fac = geom_fac,
                                eps = eps,
                                qcap = qcap,
                                gridsize = gridsize,
                                exptcap = exptcap)
        }
    }
}
