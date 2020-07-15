source("plot_Rcurve.R")

load("../data/dBH_mvgauss_Rcurve.RData")
pdf("../figs/calibrate_gamma1.pdf", width = 5, height = 5)
plot_Rcurve(Rcurves$Rcurve_BH,
            Rcurves$Rcurve_dBH,
            Rcurves$Rcurve_dBH2,
            n = n, side = side, alpha = alpha, avals = 1:n,
            cols = c("red", "blue", "orange"),
            ltys = c(1, 2, 3),
            lwd = 2,
            xlim = c(0, 40 * alpha / n), ylim = c(0.02, 0.05),
            labels = c(expression("BH"(alpha)),
                       expression("dBH"[1](alpha)),
                       expression("dBH"[1]^2*(alpha))))
dev.off()
