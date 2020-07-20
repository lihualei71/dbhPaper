library("tidyverse")
library("ggplot2")
library("latex2exp")

load("../data/dBH_mvgauss_timing_aggregate.RData")

plot <- res %>% select(-model, -fac) %>%
    unite(method, c("method", "niter", "avals_type")) %>%
    group_by(n, method, side, nalt) %>%
    summarize(time_up = quantile(time, 0.95),
              time_lo = quantile(time, 0.05),
              time = median(time)) %>%
    ungroup() %>%
    mutate(nalt = factor(nalt, levels = c(10, 30),
                         labels = paste0("#non-nulls = ", c(10, 30))),
           side = factor(side, levels = c("right", "two"),
                         labels = c("One-sided test", "Two-sided test"))) %>%
    ggplot(aes(x = n, y = time,
               color = method,
               linetype = method)) +
    geom_line() +
    facet_wrap(side ~ nalt, nrow = 1) +
    xlab("Number of hypotheses") +
    ylab("Median running time (in seconds)") +
    scale_x_continuous(trans = "log10",
                       breaks = 10^(2:6),
                       labels = c(expression(10^2),
                                  expression(10^3),
                                  expression(10^4),
                                  expression(10^5),
                                  expression(10^6))) +
    scale_y_continuous(trans = "log10",
                       breaks = 10^(seq(-2, 3, 1)),
                       labels = c(expression(10^{-2}),
                                  expression(10^{-1}),
                                  expression(10^0),
                                  expression(10^1),
                                  expression(10^2),
                                  expression(10^3))) +
    scale_color_manual(
        name = "method",
        values = c("red", "red", "blue", "blue"),
        labels = c("dBH_1_BH" = expression(dBH[1]),
                   "dBH_2_BH" = expression(dBH[1]^{2}),
                   "dBH_1_geom" = expression(s-dBH[1]),
                   "dBH_2_geom" = expression(s-dBH[1]^{2}))
    ) +
    scale_linetype_manual(
        name = "method",
        guide = "none",
        values = c("solid", "longdash", "solid", "longdash"),
        labels = c("dBH_1_BH" = expression(dBH[1]),
                   "dBH_2_BH" = expression(dBH[1]^{2}),
                   "dBH_1_geom" = expression(s-dBH[1]),
                   "dBH_2_geom" = expression(s-dBH[1]^{2}))
    ) +
    guides(color = guide_legend(
               override.aes = list(
                   linetype = c("solid", "longdash", "solid", "longdash")))) + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(size = 15),
          axis.text = element_text(size = 8.5),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.position = "bottom")

ggsave(filename = "../figs/dBH_mvgauss_timing.pdf", 
       plot, width = 10, height = 3.8)
