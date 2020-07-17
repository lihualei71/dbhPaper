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
                        beta = 2,
                        eps = 0.05,
                        qcap = 2,
                        gridsize = 20,
                        exptcap = 0.9){
    if (niter > 2){
        stop("\'niter\' can only be 1 or 2.")
    }
    if (niter == 1){
        if (tautype == "QC"){
            dBH_mvgauss_qc(zvals = zvals,
                           Sigma = Sigma,
                           Sigmafun = Sigmafun,
                           vars = vars,
                           side = side,
                           alpha = alpha,
                           gamma = gamma,
                           avals = avals,
                           avals_type = avals_type,
                           beta = beta,
                           eps = eps,
                           qcap = qcap)
        }
    } else if (niter == 2){
        if (tautype == "QC"){
            dBH_mvgauss_qc_grid(zvals = zvals,
                                Sigma = Sigma,
                                Sigmafun = Sigmafun,
                                vars = vars,
                                side = side,
                                alpha = alpha,
                                gamma = gamma,
                                avals = avals,
                                avals_type = avals_type,
                                beta = beta,
                                eps = eps,
                                qcap = qcap,
                                gridsize = gridsize,
                                exptcap = exptcap)
        }
    }
}
