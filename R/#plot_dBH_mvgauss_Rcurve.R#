source("plot_Rcurve.R")

load("../data/dBH_mvgauss_Rcurve.RData")
load("../data/dBH_mvgauss_safe_Rcurve.RData")

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

pdf("../figs/calibrate_safe.pdf", width = 5, height = 5)
plot_Rcurve(Rcurves_safe$BH,
            Rcurves_safe$dBH,
            Rcurves_safe$dBH2,
            n = Rcurves_safe$n,
            side = Rcurves_safe$side,
            alpha = Rcurves_safe$alpha,
            avals = Rcurves_safe$avals,
            cols = c("red", "blue", "orange"),
            ltys = c(1, 2, 4),
            xlim = c(0, 40 * Rcurves_safe$alpha / n),
            ylim = c(0.02, 0.05),
            labels = c(expression("BH"(alpha)),
                       expression("dBH"[1](alpha)),
                       expression("dBH"[1]^2*(alpha))))
dev.off()
