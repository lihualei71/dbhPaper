source("plot_Rcurve.R")

load("../data/dBH_mvgauss_Rcurve.RData")

pdf("../figs/calibrate_gamma1.pdf", width = 5, height = 5)
plot_Rcurve(Rcurves$BH,
            Rcurves$dBH,
            Rcurves$dBH2,
            n = Rcurves$n,
            side = Rcurves$side,
            alpha = Rcurves$alpha,
            avals = Rcurves$avals,
            cols = c("red", "blue", "orange"),
            ltys = c(1, 2, 4),
            xlim = c(0, 40 * Rcurves$alpha / n),
            ylim = c(0.02, 0.05),
            labels = c(expression("BH"(alpha)),
                       expression("dBH"[1](alpha)),
                       expression("dBH"[1]^2*(alpha))))
dev.off()
