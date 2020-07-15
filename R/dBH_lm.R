source("dBH_mvt.R")

dBH_lm <- function(y, X,
                   side = c("right", "left", "two"),
                   alpha = 0.05, alpha0 = NULL,
                   tautype = c("GF", "QC"),
                   niter = 1,                   
                   avals = NULL,
                   avals_type = c("BH", "geom", "bonf", "manual"),
                   gamma = 2,
                   knotsfun_type = NULL,
                   eps = 0.05,
                   use_cap = TRUE,
                   acc_multiplier = 2,
                   acc_minrejs = 10,
                   gridfun = lingrid,
                   gridsize = 20,
                   if_thresh_expt = TRUE,
                   thresh_expt = min(0.9 * alpha, 1.1 * alpha0)){
    stats <- lm_mvt(y, X)
    dBH_mvt(stats$zvals, stats$Sigma,
            stats$sigmahat, stats$df,
            side,
            alpha, alpha0,
            tautype,
            niter,
            avals, avals_type, gamma,
            knotsfun_type,
            eps,
            use_cap,
            acc_multiplier, acc_minrejs,
            gridfun, gridsize,
            if_thresh_expt, thresh_expt)
}
