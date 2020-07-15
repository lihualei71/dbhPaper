library("tidyverse")
library("ggplot2")
library("latex2exp")

load("../data/HIV_discoveries.RData")

methods_levels <- c("Knockoff_equi",
                    "BH_full",
                    "dBH2_QC_full_0.9",                    
                    "dBH2_QC_full_safe")
methods_labels <- c("KN", "BH",
                    parse(text = TeX("dBH$_{0.9}^2$")),
                    parse(text = TeX("dBY$^2$")))
names(methods_labels) <- methods_levels

for (al in unique(discoveries$alpha)){
    data <- discoveries %>%
        filter(alpha == al,
               method %in% methods_levels) %>%
        mutate(method = factor(method,
                               levels = methods_levels)) %>%
        mutate(drug = paste0(drug_name, " (", drug_class, ")"))
    drug_levels <- unique(data$drug)
    plot <- data %>%
        select(-alpha, -drug_name, -drug_class, -secBH) %>%
        gather("discoveries", "value", -drug, -method) %>%
        mutate(discoveries = factor(
                   discoveries,
                   levels = c("nfalse", "ntrue"),
                   labels = c("Not in TSM list", "In TSM list")),
               drug = factor(drug,
                             levels = drug_levels)) %>%
        ggplot(aes(x = method, y = value)) +
        geom_bar(stat = "identity", aes(fill = discoveries)) +
        facet_wrap(~ drug, nrow = 4) +
        scale_x_discrete(labels = methods_labels) +
        scale_fill_manual(values = c("orangered3", "navyblue")) +
        ylab("Number of discoveries") + 
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.text.x = element_text(angle = 90),
              strip.text = element_text(size = 10),
              axis.text = element_text(size = 10),
              axis.title = element_text(size = 12.5),
              legend.position = "bottom")
    ggsave(filename = paste0("../figs/dBH-HIV-", al, "-paper.pdf"),
           plot, width = 6.5, height = 7.5)
        
}
