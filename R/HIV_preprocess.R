#### This R script is a part of the tutorial in the "knockoff" package
#### Source: https://cran.r-project.org/web/packages/knockoff/vignettes/hiv.html

data <- list()
for (drug_class in c("PI", "NRTI", "NNRTI")){
    base_url <- 'http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006'
    gene_url <- paste(base_url, 'DATA', paste0(drug_class, '_DATA.txt'), sep = '/')
    tsm_url <- paste(base_url, 'MUTATIONLISTS', 'NP_TSM', drug_class, sep = '/')

    gene_df <- read.delim(gene_url, na.string = c('NA', ''), stringsAsFactors = FALSE)
    tsm_df <- read.delim(tsm_url, header = FALSE, stringsAsFactors = FALSE)
    names(tsm_df) <- c('Position', 'Mutations')

    grepl_rows <- function(pattern, df) {
        cell_matches = apply(df, c(1,2), function(x) grepl(pattern, x))
        apply(cell_matches, 1, all)
    }

    pos_start <- which(names(gene_df) == 'P1')
    pos_cols <- seq.int(pos_start, ncol(gene_df))
    valid_rows <- grepl_rows('^(\\.|-|[A-Zid]+)$', gene_df[,pos_cols])
    gene_df <- gene_df[valid_rows,]


    flatten_matrix <- function(M, sep = '.'){
        x <- c(M)
        names(x) <- c(outer(rownames(M), colnames(M),
                            function(...) paste(..., sep=sep)))
        x
    }

    ## Construct preliminary design matrix.
    muts <- c(LETTERS, 'i', 'd')
    X <- outer(muts, as.matrix(gene_df[,pos_cols]), Vectorize(grepl))
    X <- aperm(X, c(2,3,1))
    dimnames(X)[[3]] <- muts
    X <- t(apply(X, 1, flatten_matrix))
    mode(X) <- 'numeric'

    ## Remove any mutation/position pairs that never appear in the data.
    X <- X[,colSums(X) != 0]

    ## Extract response matrix.
    Y <- gene_df[,4:(pos_start-1)]

    ## Save data
    data[[drug_class]] <- list(Y = Y, X = X)
}
save(data, file = "../data/HIV_data.RData")
