library("dplyr")
source("expr_functions.R")

params <- read.table("../jobs/dBH_mvgauss_timing_params.txt")
names(params) <- c("seed", "nalt", "rho", "fac", "nreps", "dBH2")
settings <- params %>% select(-seed) %>% unique

aggregate_res <- list()
for (i in 1:nrow(settings)){
    seedlist <- params %>%
        filter(nalt == settings$nalt[i],
               rho == settings$rho[i],
               fac == settings$fac[i],
               nreps == settings$nreps[i],
               dBH2 == settings$dBH2[i]) %>%
        .$seed
    file_root <- paste0("../cluster_raw_data/dBH_mvgauss_timing",
                        "_nalt", settings$nalt[i],
                        "_rho", settings$rho[i],
                        "_fac", settings$fac[i],
                        "_nreps", settings$nreps[i],
                        "_dBH2", settings$dBH2[i])

    expr_list <- list()
    k <- 1
    for (seed in seedlist){
        tmp <- try(load(paste0(file_root, "_seed", seed, ".RData")))
        if (class(tmp) != "try-error"){
            expr_list[[k]] <- res
            k <- k + 1
        }
    }
    aggregate_res[[i]] <- do.call(rbind, expr_list)
}

res <- do.call(rbind, aggregate_res)
save(res, file = "../data/dBH_mvgauss_timing_aggregate.RData")
