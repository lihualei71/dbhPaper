library("tidyverse")
library("ggplot2")
library("latex2exp")

load("../data/dBH_mcc_aggregate.RData")

methods_levels <- c("BH_full",
                    "BH_sparse(2)",
                    "dBH_QC_full_0.9",
                    "dBH_QC_sparse(2)_0.9",
                    "dBH_QC_full_ 1",
                    "dBH_QC_sparse(2)_ 1",
                    "dBH2_QC_full_0.9",
                    "dBH2_QC_sparse(2)_0.9",
                    "dBH2_QC_full_ 1",
                    "dBH2_QC_sparse(2)_ 1",
                    "BH_full_safe",
                    "BH_sparse(2)_safe",
                    "dBH_QC_full_safe",
                    "dBH_QC_sparse(2)_safe",
                    "dBH2_QC_full_safe",
                    "dBH2_QC_sparse(2)_safe",
                    "Knockoff_equi",
                    "Knockoff_sdp")
safe_methods_levels <- methods_levels[11:16]
knockoff_methods_levels <- methods_levels[17:18]
methods_labels <- c("BH", "s-BH",
                    parse(text = TeX("dBH$_{0.9}$")),
                    parse(text = TeX("s-dBH$_{0.9}$")),
                    parse(text = TeX("dBH$_{1}$")),
                    parse(text = TeX("s-dBH$_{1}$")),
                    parse(text = TeX("dBH$_{0.9}^2$")),
                    parse(text = TeX("s-dBH$_{0.9}^2$")),
                    parse(text = TeX("dBH$_{1}^2$")),
                    parse(text = TeX("s-dBH$_{1}^2$")),
                    "BY", "s-BY",
                    "dBY", "s-dBY",
                    parse(text = TeX("dBY$^2$")),
                    parse(text = TeX("s-dBY$^2$")),
                    "KN-eq", "KN-sdp")

names(methods_labels) <- methods_levels
methods_colors <- c("grey65", "grey65", 
                    "dodgerblue3", "dodgerblue3", "dodgerblue3", "dodgerblue3",
                    "indianred3", "indianred3", "indianred3", "indianred3",
                    "gray", "gray",
                    "dodgerblue1", "dodgerblue1",
                    "indianred1", "indianred1",
                    "orange1", "orange3")

color_df <- data.frame(method = methods_levels,
                       color = methods_colors)
    
settings <- res %>%
    select(ng, nr, pi1, side, alpha) %>%
    unique
secBH <- data.frame()    

for (i in 1:nrow(settings)){
    data <- res %>%
        filter(ng == settings$ng[i],
               nr == settings$nr[i],
               pi1 == settings$pi1[i],
               side == settings$side[i],
               alpha == settings$alpha[i],
               method %in% methods_levels) %>%
        mutate(safe = method %in% safe_methods_levels,
               knockoff = method %in% knockoff_methods_levels,
               method_type = 2 * safe + knockoff) %>%
        select(-safe, -knockoff) %>%
        left_join(color_df, by = "method") %>%
        mutate(method = factor(
                   method,
                   levels = methods_levels)) %>%
            mutate(method_type =
                       factor(method_type,
                              levels = c(0, 2, 1),
                              labels = c("BH-type methods", "BY-type methods", "KN")))
    
    file_root <- paste0("../figs/dBH-mcc",
                        "-ng", settings$ng[i],
                        "-nr", settings$nr[i],
                        "-pi1", settings$pi1[i],
                        "-side", settings$side[i],
                        "-alpha", settings$alpha[i])

    reference <- data.frame(type = c("FDR", "power"),
                            value = c(settings$alpha[i], NA),
                            linetype = c("longdashed", "dotted"),
                            color = c("black", "grey45"))
    plot <- data %>%
        select(-`FDR (init)`, -`power (init)`, -secBH,
               -qmax, -q99, -q95,
               -ng, -nr, -pi1, -side, -alpha) %>%
        gather("type", "value", -method, -method_type, -color) %>%
        ggplot(aes(x = method, y = value,
                   fill = color, color = color)) +
        geom_bar(stat = "identity", width = 0.8) +
        facet_grid(type ~ method_type, scales = "free_x", space = "free_x") +
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
    ggsave(filename = paste0(file_root, ".pdf"),
           plot, width = 6.5, height = 5.5)

    secBH <- data %>% select(ng, nr, side, alpha,
                             method, secBH) %>%
        rbind(secBH, .)
}

filename <- paste0("../data/dBH_mcc_secBH.RData")
save(secBH, file = filename)
