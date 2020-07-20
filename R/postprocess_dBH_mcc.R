library("dplyr")
source("expr_functions.R")

params <- read.table("../jobs/dBH_mcc_params.txt")
names(params) <- c("seed", "ng", "nr", "pi1", "mutype", "side", "nreps", "dBH2", "knockoff")
settings <- params %>% select(-seed) %>% unique
alphas <- c(0.05, 0.2)

aggregate_res <- list()
for (i in 1:nrow(settings)){
    seedlist <- params %>%
        filter(ng == settings$ng[i],
               nr == settings$nr[i],
               pi1 == settings$pi1[i],
               mutype == settings$mutype[i],
               side == settings$side[i],
               nreps == settings$nreps[i],
               dBH2 == settings$dBH2[i],
               knockoff == settings$knockoff[i]) %>%
        .$seed
    file_root <- paste0("../cluster_raw_data/dBH_mcc",
                        "_ng", settings$ng[i],
                        "_nr", settings$nr[i],
                        "_pi1", settings$pi1[i],
                        "_mutype", settings$mutype[i],
                        "_side", settings$side[i],
                        "_nreps", settings$nreps[i],
                        "_dBH2", settings$dBH2[i],
                        "_knockoff", settings$knockoff[i])
    expr_list <- list()
    k <- 1
    
    for (seed in seedlist){
        tmp <- try(load(paste0(file_root, "_seed", seed, ".RData")))
        if (class(tmp) != "try-error"){        
            expr_list[[k]] <- res
        }
        k <- k + 1
    }

    expr_res <- aggregate_expr(expr_list, alphas) %>%
        postprocess
    aggregate_res[[i]] <- expr_res %>%
        mutate(ng = settings$ng[i],
               nr = settings$nr[i],
               pi1 = settings$pi1[i],
               side = settings$side[i])
}

res <- do.call(rbind, aggregate_res)
save(res, file = "../data/dBH_mcc_aggregate.RData")

print(summary(res$qmax[res$qmax > -Inf]))
print(summary(res$q99[res$q99 > -Inf]))
print(summary(res$q95[res$q95 > -Inf]))
