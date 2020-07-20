library("dbh")
source("utils.R")
source("expr_functions.R")

if (!file.exists("../data/HIV_data.RData")){
    source("HIV_preprocess.R")
}
load("../data/HIV_data.RData")

HIV_expr <- function(X, y, 
                     side,
                     alphas, 
                     alpha_fac = 0.9,
                     geom_fac = 2,
                     tautype = "QC",
                     skip_knockoff = TRUE,
                     skip_dBH2 = TRUE,
                     ...){
    ## Log-transform the drug resistance measurements.
    y <- log(y)

    ## Remove patients with missing measurements.
    missing <- is.na(y)
    y <- y[!missing]
    X <- X[!missing,]

    ## Remove predictors that appear less than 3 times.
    X <- X[,colSums(X) >= 3]

    ## Remove duplicate predictors.
    X <- X[,colSums(abs(cor(X)-1) < 1e-4) == 1]

    ## Get names
    genes <- colnames(X)

    ## Get stats
    nalphas <- length(alphas)
    n <- nrow(X)
    p <- ncol(X)
    stats <- lm_mvt(y, X, 1:p, intercept = FALSE)
    tvals <- stats$tvals
    df <- stats$df
    Sigma <- stats$Sigma
    pvals <- pvals_mvt(tvals, Sigma, df, side)    
    
    ## Get methods and experimental settings for dBH
    methods <- gen_methods(alpha_fac, geom_fac, tautype,
                           skip_knockoff, skip_dBH2)
    expr_params <- expand.grid(
        alpha_fac = alpha_fac,
        geom_fac = geom_fac,
        tautype = tautype
    )

    if (!skip_knockoff){
        if (n < 2 * p){
            obj <- lm(y ~ X)
            sigma <- summary(obj)$sigma
            nvars <- 2 * p
            Xnew <- matrix(0, nrow = nvars - n, ncol = p)
            ynew <- rnorm(nvars - n) * sigma
            X <- rbind(X, Xnew)
            y <- c(y, ynew)
            print("n < 2p. Approximated knockoff is used.")
        }
        Xk_equi <- knockoff::create.fixed(X, "equi")$Xk
        Xk_sdp <- knockoff::create.fixed(X, "sdp")$Xk
    }
    
    res <- list()
    pb <- txtProgressBar(style=3)
    for (k in 1:nalphas){
        alpha <- alphas[k]
        res[[k]] <- list(alpha = alpha, rejs = list())
        obj <- list()

        ## BH rejections
        for (x in union(NA, geom_fac)){
            if (is.na(x)){
                avals <- 1:p
            } else {
                avals <- geom_avals(x, p)
            }
            rejs_BH <- BH(pvals, alpha, avals, FALSE)
            rejs_BH$rejs <- genes[rejs_BH$rejs]
            rejs_BH_safe <- BH(pvals, alpha, avals, TRUE) 
            rejs_BH_safe$rejs <- genes[rejs_BH_safe$rejs]
            obj <- c(obj, list(rejs_BH, rejs_BH_safe))
        }

        ## BC rejections
        rejs_BC <- BC(pvals, alpha)
        rejs_BC$rejs <- genes[rejs_BC$rejs]        
        obj <- c(obj, list(rejs_BC))

        ## Knockoff rejections
        if (!skip_knockoff) {        
            rejs_knockoff_equi <- knockoff_lm(X, Xk_equi, y, alpha)
            rejs_knockoff_equi$rejs <- genes[rejs_knockoff_equi$rejs]
            rejs_knockoff_sdp <- knockoff_lm(X, Xk_sdp, y, alpha)
            rejs_knockoff_sdp$rejs <- genes[rejs_knockoff_sdp$rejs]
            obj <- c(obj, list(rejs_knockoff_equi, rejs_knockoff_sdp))
        }
        
        ## Number of methods so far
        nBHBCkn <- length(obj)
        
        ## dBH rejections
        for (j in 1:nrow(expr_params)){
            fac <- expr_params[j, 1]
            x <- expr_params[j, 2]
            type <- expr_params[j, 3]
            if (is.na(x)){
                avals_type <- "BH"
            } else {
                avals_type <- "geom"
            }
            if (is.na(fac)){
                gamma <- NULL
            } else {
                gamma <- fac
            }
            rejs_dBH <- dBH_mvt(
                tvals = tvals,
                df = df,
                Sigma = Sigma,
                side = side,
                alpha = alpha,
                gamma = gamma, 
                niter = 1,
                tautype = type,
                avals_type = avals_type,
                geom_fac = x, ...)
            rejs_dBH$rejs <- genes[rejs_dBH$rejs]
            rejs_dBH_init <- list(rejs = rejs_dBH$initrejs)
            rejs_dBH_init$rejs <- genes[rejs_dBH_init$rejs]
            obj <- c(obj, list(rejs_dBH, rejs_dBH_init))
        }

        if (!skip_dBH2){
            ## dBH2 rejections
            for (j in 1:nrow(expr_params)){
                fac <- expr_params[j, 1]
                x <- expr_params[j, 2]
                type <- expr_params[j, 3]
                if (is.na(x)){
                    avals_type <- "BH"
                } else {
                    avals_type <- "geom"
                }
                if (is.na(fac)){
                    alpha0 <- NULL
                } else {
                    alpha0 <- fac * alpha
                }
                rejs_dBH2 <- dBH_mvt(
                    tvals = tvals,
                    df = df,
                    Sigma = Sigma,
                    side = side,
                    alpha = alpha,
                    gamma = gamma,
                    niter = 2,
                    tautype = type,
                    avals_type = avals_type,
                    geom_fac = x, ...)
                rejs_dBH2$rejs <- genes[rejs_dBH2$rejs]
                rejs_dBH2_init <- list(rejs = rejs_dBH2$initrejs)
                rejs_dBH2_init$rejs <- genes[rejs_dBH2_init$rejs]
                obj <- c(obj, list(rejs_dBH2, rejs_dBH2_init))
                setTxtProgressBar(pb, ((k - 1) * nrow(expr_params) + j) / (nalphas * nrow(expr_params)))
            }
        }
        names(obj) <- methods
        res[[k]]$rejs <- obj
    }
    
    return(res)
}

set.seed(20200711)
alpha_fac <- c(0.9, NA)
gamma <- c(NA, 2)
alphas <- c(0.05, 0.2)
tautype <- "QC"
side <- "two"
skip_knockoff <- FALSE
skip_dBH2 <- FALSE
res <- list()

for (drug_class in names(data)){
    Y <- data[[drug_class]]$Y
    X <- data[[drug_class]]$X
    res[[drug_class]] <- lapply(1:length(Y), function(j){
        print(j)
        HIV_expr(X, Y[[j]],
                 side, alphas,
                 alpha_fac,
                 gamma,
                 tautype,
                 skip_knockoff,
                 skip_dBH2)
    })
}
save(file = "../data/HIV_res.RData", res)
