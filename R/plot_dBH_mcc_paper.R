library("tidyverse")
library("ggplot2")
library("latex2exp")

load("../data/dBH_mcc_aggregate.RData")

methods_levels <- c("BH_full",    
                    "dBH_QC_full_0.9",
                    "dBH_QC_full_ 1",
                    "dBH2_QC_full_0.9",
                    "dBH2_QC_full_ 1",
                    "BH_full_safe",
                    "dBH_QC_full_safe",
                    "dBH2_QC_full_safe",
                    "Knockoff_sdp")
methods_labels <- c("BH",
                    parse(text = TeX("dBH$_{0.9}$")),
                    parse(text = TeX("dBH$_{1}$")),
                    parse(text = TeX("dBH$_{0.9}^2$")),
                    parse(text = TeX("dBH$_{1}^2$")),
                    "BY", "dBY",
                    parse(text = TeX("dBY$^2$")),
                    "KN")
names(methods_labels) <- methods_levels
methods_colors <- c("grey65",
                    "dodgerblue3", "dodgerblue3",
                    "indianred3", "indianred3",
                    "gray",
                    "dodgerblue1",
                    "indianred1",
                    "orange1")

color_df <- data.frame(method = methods_levels,
                       color = methods_colors)
    
settings <- res %>%
    filter(side == "two") %>%
    select(ng, nr, pi1, side) %>%
    unique
secBH <- data.frame()    

for (i in 1:nrow(settings)){
    data <- res %>%
        filter(ng == settings$ng[i],
               nr == settings$nr[i],
               pi1 == settings$pi1[i],
               side == settings$side[i],
               method %in% methods_levels) %>%
        left_join(color_df, by = "method") %>%
        mutate(method = factor(
                   method,
                   levels = methods_levels
               ), 
               alpha = factor(alpha,
                                  levels = c(0.05, 0.2),
                                  labels = c("alpha = 0.05", "alpha = 0.2")))
    
    file_root <- paste0("../figs/dBH-mcc",
                        "-ng", settings$ng[i],
                        "-nr", settings$nr[i],
                        "-pi1", settings$pi1[i])

    reference <- data.frame(type = c("FDR", "FDR", "power", "power"),
                            alpha = c("alpha = 0.05", "alpha = 0.2", "alpha = 0.05", "alpha = 0.2"),
                            value = c(0.05, 0.2, NA, NA))
    plot <- data %>%
        select(-`FDR (init)`, -`power (init)`, -secBH,
               -qmax, -q99, -q95,
               -ng, -nr, -pi1, -side) %>%
        gather("type", "value", -method, -alpha, -color) %>%
        ggplot(aes(x = method, y = value,
                   fill = color, color = color)) +
        geom_bar(stat = "identity", width = 0.8) +
        facet_grid(type ~ alpha, scales = "free_x", space = "free_x") +
        scale_colour_identity() +
        scale_fill_identity() +
        scale_x_discrete(labels = methods_labels) +
        scale_y_continuous(breaks = c(0.05, 0.2, 0.35, 0.5, 0.65)) + 
        geom_hline(data = reference, aes(yintercept = value), linetype = "longdash") +
        ylab("") + 
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.text.x = element_text(angle = 90),
              strip.text = element_text(size = 17.5),
              axis.text = element_text(size = 12.5),
              axis.title = element_text(size = 17.5),
              legend.position = "bottom")
    ggsave(filename = paste0(file_root, "-paper.pdf"),
           plot, width = 6.5, height = 5.5)
}
