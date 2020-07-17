source("dBH_mvt.R")

dBH_lm <- function(y, X,
                   subset = 1:ncol(X),
                   side = c("right", "left", "two"),
                   alpha = 0.05, gamma = NULL,
                   tautype = "QC",
                   niter = 1,                   
                   avals = NULL,
                   avals_type = c("BH", "geom", "bonf", "manual"),
                   beta = 2,
                   eps = 0.05,
                   qcap = TRUE,
                   gridsize = 20,
                   exptcap = 0.9){
    stats <- lm_mvt(y, X)
    dBH_mvt(tvals = stats$tvals, 
            df = stats$df,
            Sigma = stats$Sigma,
            Sigmafun = NULL,
            vars = NULL,
            side = side,
            alpha = alpha,
            gamma = gamma,
            tautype = tautype,
            niter = niter,
            avals = avals,
            avals_type = avals_type,
            beta = beta,
            eps = eps,
            qcap = qcap,
            gridsize = gridsize,
            exptcap = exptcap)
}
