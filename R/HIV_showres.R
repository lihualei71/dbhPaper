drug_class <- 'PI'
base_url <- 'http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006'
gene_url <- paste(base_url, 'DATA', paste0(drug_class, '_DATA.txt'), sep = '/')
tsm_url <- paste(base_url, 'MUTATIONLISTS', 'NP_TSM', drug_class, sep = '/')

gene_df <- read.delim(gene_url, na.string = c('NA', ''), stringsAsFactors = FALSE)
tsm_df <- read.delim(tsm_url, header = FALSE, stringsAsFactors = FALSE)
names(tsm_df) <- c('Position', 'Mutations')

get_position <- function(x)
  sapply(regmatches(x, regexpr('[[:digit:]]+', x)), as.numeric)

HIV_showres <- function(res){
    lapply(res, function(drug) {
        lapply(drug$rejs, function(method) {
            positions <- get_position(method) # remove possible duplicates
            discoveries <- length(positions)
            false_discoveries <- length(setdiff(positions, tsm_df$Position))
            list(true_discoveries = discoveries - false_discoveries,
                 false_discoveries = false_discoveries,
                 fdp = false_discoveries / max(1, discoveries))
        })
    })
}

methods <- c("Knockoff", "BHq", "dBH_0.8", "dBH2_0.8",
             "BY", "dBH_safe", "dBH2_safe")
methods_name <- c("Knockoff", "BHq", "dBH_0.8", "dBH2_0.8",
             "BY", "dBY", "dBY2")

load("../data/HIV_res.RData")
comparisons <- HIV_showres(res)
for (drug in names(comparisons)){
    comparisons[[drug]] <- comparisons[[drug]][methods]
    names(comparisons[[drug]]) <- methods_name
}
for (drug in names(comparisons)) {
    plot_data <- do.call(cbind, comparisons[[drug]])
    plot_data <- plot_data[c('true_discoveries','false_discoveries'),]
    pdf(paste0("../talks/Berkeley_Columbia_2020/figs/HIV_alpha0.2_", drug, ".pdf"), width = 8, height = 6)
    barplot(as.matrix(plot_data), main = paste('Resistance to', drug),
            col = c('navy','orange'))
    dev.off()
}

load("../data/HIV_alpha0.1.RData")
comparisons <- HIV_showres(res)
for (drug in names(comparisons)){
    comparisons[[drug]] <- comparisons[[drug]][methods]
    names(comparisons[[drug]]) <- methods_name
}
for (drug in names(comparisons)) {
    plot_data <- do.call(cbind, comparisons[[drug]])
    plot_data <- plot_data[c('true_discoveries','false_discoveries'),]
    pdf(paste0("../talks/Berkeley_Columbia_2020/figs/HIV_alpha0.1_", drug, ".pdf"), width = 8, height = 6)
    barplot(as.matrix(plot_data), main = paste('Resistance to', drug),
            col = c('navy','orange'))
    dev.off()
}

load("../data/HIV_alpha0.05.RData")
comparisons <- HIV_showres(res)
for (drug in names(comparisons)){
    comparisons[[drug]] <- comparisons[[drug]][methods]
    names(comparisons[[drug]]) <- methods_name
}
for (drug in names(comparisons)) {
    plot_data <- do.call(cbind, comparisons[[drug]])
    plot_data <- plot_data[c('true_discoveries','false_discoveries'),]
    pdf(paste0("../talks/Berkeley_Columbia_2020/figs/HIV_alpha0.05_", drug, ".pdf"), width = 8, height = 6)
    barplot(as.matrix(plot_data), main = paste('Resistance to', drug),
            col = c('navy','orange'))
    dev.off()
}

