source("Rcurve.R")

id <- 20
n <- 1000
alpha <- 0.05
side <- "one"
ntails <- ifelse(side == "two", 2, 1)

set.seed(20200212)
pi1 <- 0.1
mu1 <- 2.5
rho <- 0.8
Sigma <- genSigma(n, rho, type = "AR")
mu <- genmu(n, pi1, mu1, "fix")
zvals <- mvtnorm::rmvnorm(1, mu, Sigma)
zvals <- as.numeric(zvals)
low <- qnorm(alpha / ntails, lower.tail = FALSE)
high <- qnorm(alpha / n / ntails * 0.01, lower.tail = FALSE)

Rcurve_BH <- Rcurve_mvgauss(id, zvals, Sigma, side, alpha,
                            alpha0 = alpha,
                            avals_type = "BH",
                            low = low, high = high,
                            niter = 0)

Rcurve_dBH <- Rcurve_mvgauss(id, zvals, Sigma, side, alpha,
                             alpha0 = alpha,
                             avals_type = "BH",
                             low = low, high = high,
                             niter = 1)

Rcurve_dBH2 <- Rcurve_mvgauss(id, zvals, Sigma, side, alpha,
                              alpha0 = alpha,
                              avals_type = "BH",
                              low = low, high = high,
                              niter = 2)

plot_Rcurve(Rcurve_BH, Rcurve_dBH, Rcurve_dBH2,
            n = n, side = side, alpha = alpha, avals = 1:n,
            cols = c("red", "blue", "orange"),
            ltys = c(1, 2, 3),
            xlim = c(0, 40), ylim = c(0.02, 0.05),
            labels = c(expression("BH"(alpha)),
                       expression("dBH"[1](alpha)),
                       expression("dBH"[1]^2*(alpha))))
