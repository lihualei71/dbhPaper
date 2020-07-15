source("dBH_mvt_gf.R")
source("dBH_mvt_qc.R")
source("dBH_mvt_gf_grid.R")
source("dBH_mvt_qc_grid.R")

dBH_mvt <- function(zvals, Sigma,
                    sigmahat, df,
                    side = c("right", "left", "two"),
                    alpha = 0.05, alpha0 = NULL,
                    tautype = c("GF", "QC"),
                    niter = 1,
                    avals = NULL, 
                    avals_type = c("BH", "geom", "bonf", "manual"),
                    gamma = 2,
                    eps = 0.05,
                    use_cap = TRUE,
                    acc_multiplier = 2,
                    acc_minrejs = 10,
                    gridfun = lingrid,
                    gridsize = 20,
                    if_thresh_expt = TRUE,
                    thresh_expt = min(0.9 * alpha, 1.1 * alpha0)
                    ){
    if (niter > 2){
        stop("\'niter\' can only be 1 or 2.")
    }
    tautype <- tautype[1]
    if (niter == 1){
        if (tautype == "GF"){
            dBH_mvt_gf(zvals, Sigma,
                       sigmahat, df,
                       side,
                       alpha, alpha0,
                       avals, avals_type, gamma,
                       eps,
                       use_cap,
                       acc_multiplier, acc_minrejs)
        } else if (tautype == "QC"){
            dBH_mvt_qc(zvals, Sigma,
                       sigmahat, df,
                       side,
                       alpha, alpha0,
                       avals, avals_type, gamma,
                       eps,
                       use_cap,
                       acc_multiplier, acc_minrejs)
        }
    } else if (niter == 2){
        if (tautype == "GF"){
            dBH_mvt_gf_grid(zvals, Sigma,
                            sigmahat, df,
                            side,
                            alpha, alpha0,
                            avals, avals_type, gamma,
                            eps,
                            gridfun = gridfun,
                            gridsize = gridsize,
                            if_thresh_expt = if_thresh_expt,
                            thresh_expt = thresh_expt,
                            use_cap = use_cap,
                            acc_multiplier = acc_multiplier,
                            acc_minrejs = acc_minrejs)
        } else if (tautype == "QC"){
            dBH_mvt_qc_grid(zvals, Sigma,
                            sigmahat, df,
                            side,
                            alpha, alpha0,
                            avals, avals_type, gamma,
                            eps,
                            gridfun = gridfun,
                            gridsize = gridsize,
                            if_thresh_expt = if_thresh_expt,
                            thresh_expt = thresh_expt,
                            use_cap = use_cap,
                            acc_multiplier = acc_multiplier,
                            acc_minrejs = acc_minrejs)
        }
    }
}
