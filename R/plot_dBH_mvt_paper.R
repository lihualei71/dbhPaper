library("tidyverse")
library("ggplot2")
library("latex2exp")

load("../data/dBH_mvt_aggregate.RData")

methods_levels <- c("BH_full",    
                    "dBH_QC_full_0.9",
                    "dBH_QC_full_ 1", 
                    "dBH2_QC_full_0.9",
                    "dBH2_QC_full_ 1",
                    "BH_full_safe",
                    "dBH_QC_full_safe",
                    "dBH2_QC_full_safe")
methods_labels <- c("BH",
                    parse(text = TeX("dBH$_{0.9}$")),
                    parse(text = TeX("dBH$_{1}$")),
                    parse(text = TeX("dBH$_{0.9}^2$")),
                    parse(text = TeX("dBH$_{1}^2$")),
                    "BY", "dBY",
                    parse(text = TeX("dBY$^2$")))

names(methods_labels) <- methods_levels
methods_colors <- c("grey65",
                    "dodgerblue3", "dodgerblue3",
                    "indianred3", "indianred3",
                    "gray",
                    "dodgerblue1",
                    "indianred1")

color_df <- data.frame(method = methods_levels,
                       color = methods_colors)

settings <- res %>%
    select(model, n, df, pi1, alpha) %>%
    unique
secBH <- data.frame()    

for (i in 1:nrow(settings)){
    data <- res %>%
        filter(model == settings$model[i],
               n == settings$n[i],
               df == settings$df[i],
               pi1 == settings$pi1[i],
               alpha == settings$alpha[i],
               method %in% methods_levels) %>%
        left_join(color_df, by = "method") %>%
        mutate(method = factor(method,
                               levels = methods_levels),
               side = factor(side,
                             levels = c("right", "two"),
                             labels = c("one-sided", "two-sided")))
    
    model_tmp <- switch(settings$model[i],
                        `AR (0.8)` = "AR0.8",
                        `iid` = "iid",
                        `block (0.5)` = "block0.5")
    file_root <- paste0("../figs/dBH-mvt",
                        "-", model_tmp,
                        "-n", settings$n[i],
                        "-df", settings$df[i],
                        "-pi1", settings$pi1[i],
                        "-alpha", settings$alpha[i])

    reference <- data.frame(type = c("FDR", "power"),
                            value = c(settings$alpha[i], NA))

    plot <- data %>%
        select(-`FDR (init)`, -`power (init)`, -secBH,
               -model, -n, -df, -pi1, -alpha) %>%
        gather("type", "value", -method, -color, -side) %>%
        ggplot(aes(x = method, y = value,
                   fill = color, color = color)) +
        geom_bar(stat = "identity", width = 0.8) +
        facet_grid(type ~ side, scales = "free_x") +
        scale_colour_identity() +
        scale_fill_identity() +
        scale_x_discrete(labels = methods_labels) +
        scale_y_continuous(breaks = c(0.05, 0.15, 0.25, 0.35)) + 
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
