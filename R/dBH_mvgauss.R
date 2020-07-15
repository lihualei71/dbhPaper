source("dBH_mvgauss_qc.R")
source("dBH_mvgauss_qc_grid.R")

dBH_mvgauss <- function(zvals, Sigma,
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
            dBH_mvgauss_qc(zvals, Sigma, side,
                           alpha, gamma,
                           avals, avals_type, beta,
                           eps,
                           qcap)
        }
    } else if (niter == 2){
        if (tautype == "QC"){
            dBH_mvgauss_qc_grid(zvals, Sigma, side,
                                alpha, gamma,
                                avals, avals_type, beta,
                                eps,
                                qcap,
                                gridsize,
                                exptcap)
        }
    }
}
