source("Rcurve.R")
source("utils.R")

## Plot of the 1/R curve as a function of the rescaled p-value
## p*n/alpha and the reference curve p*n/alpha = 1/u(R) where
## u is the reshaping function.
##
## Inputs:
##    ...: objects produced by Rcurve_mvgauss/mvt/lm. See
##         Rcurve.R for details
##    n: number of hypotheses
##    side: "one" or "two"
##    alpha: target FDR level
##    avals: a values in the reshaping function
##    xlim, ylim: the limits of x/y-axis. Default values
##                will be computed if NULL
##    cols: a vector of colors of length the same as the
##          number of objects in ... for 1/R curves
##    refcols: a vector of colors of length the same as the
##          number of objects in ... for reference curves
##    ltys: a vector of linetypes of length the same as the
##          number of objects in ... for 1/R curves
##    refltys: a vector of linetypes of length the same as
##          the number of objects in ... for reference curves
##    labels: a vector of legend labels of length the same as
##          the number of objects in ... for 1/R curves
##    reflabels: a vector of legend labels of length the same
##          as the number of objects in ... for reference curves
##    title: plot title
plot_Rcurve <- function(...,
                        n, side, alpha, avals, 
                        xlim = NULL, ylim = NULL,
                        cols = NULL, refcols = NULL,
                        ltys = NULL, refltys = NULL,
                        labels = NULL, reflabels = NULL,
                        title = ""){
    rescale <- n / alpha    
    ntails <- ifelse(side == "two", 2, 1)
    par(mfrow = c(1, ntails))
    res_list <- list(...)
    nobjs <- length(res_list)
    if (is.null(cols)){
        cols <- 2:(nobjs + 1)
    }
    if (is.null(refcols)){
        refcols <- rep(1, nobjs)
    }
    if (is.null(ltys)){
        ltys <- rep(1, nobjs)
    }
    if (is.null(refltys)){
        refltys <- rep(1, nobjs)
    }

    low <- min(sapply(res_list, function(res){
        min(res[[1]]$knots)
    })) * rescale
    high <- max(sapply(res_list, function(res){
        max(res[[1]]$knots)
    })) * rescale

    for (j in 1:ntails){
        refx <- fill_int(n, avals) / ntails
        refy <- 1 / (1:n)
        if (j == 1){
            if (ntails == 2){
                main <- paste(title, "(right tail)")
            } else {
                main <- title
            }
        } else {
            main <- paste(title, "(left tail)")
        }
        if (is.null(xlim)){
            xlim <- c(low, high)
        }
        if (is.null(ylim)){
            ylim <- c(0, 1)
        }
        plot(xlim, ylim, type = "n",
             xlab = expression("Rescale p-value " (p %*% n / alpha)),
             ylab = "1 / R",
             main = main,
             xlim = xlim, ylim = ylim)

        for (i in 1:nobjs){
            refx_discount <- refx * res_list[[i]]$alpha0 / alpha
            lines(refx_discount, refy, type = "l", col = refcols[i], lty = refltys[i])
            res <- res_list[[i]][[j]]
            x <- rep(res$knots, each = 2)
            m <- length(x)
            x[seq(1, m, 2)] <- x[seq(1, m, 2)] - 1e-10
            x[seq(2, m, 2)] <- x[seq(1, m, 2)] + 1e-10
            x <- x[-c(1, m)] * rescale
            y <- rep(1 / res$nrejs, each = 2)
            lines(x, y, col = cols[i], lty = ltys[i])
        }

        if (!is.null(labels)){
            legend("topright", legend = labels, col = cols, lty = ltys, cex = 0.8, y.intersp = 1.5)
        }
    }
}
