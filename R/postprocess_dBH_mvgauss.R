library("dplyr")
source("expr_functions.R")

params <- read.table("../jobs/dBH_mvgauss_params.txt")
names(params) <- c("seed", "n", "pi1", "mutype", "side", "nreps", "dBH2")
settings <- params %>% select(-seed) %>% unique
alphas <- 0.05

aggregate_res <- list()
for (i in 1:nrow(settings)){
    seedlist <- params %>%
        filter(n == settings$n[i],
               pi1 == settings$pi1[i],
               mutype == settings$mutype[i],
               side == settings$side[i],
               nreps == settings$nreps[i],
               dBH2 == settings$dBH2[i]) %>%
        .$seed
    file_root <- paste0("../cluster_raw_data/dBH_mvgauss",
                        "_n", settings$n[i],
                        "_pi1", settings$pi1[i],
                        "_mutype", settings$mutype[i],
                        "_side", settings$side[i],
                        "_nreps", settings$nreps[i],
                        "_dBH2", settings$dBH2[i])
    expr1_list <- list()
    expr2_list <- list()
    expr3_list <- list()
    k1 <- 1
    k2 <- 1
    k3 <- 1
    
    for (seed in seedlist){
        tmp <- try(load(paste0(file_root, "_seed", seed, "_AR(0.8)_fix.RData")))
        if (class(tmp) != "try-error"){        
            expr1_list[[k1]] <- res
            k1 <- k1 + 1
        }
        tmp <- try(load(paste0(file_root, "_seed", seed, "_AR(-0.8)_fix.RData")))
        if (class(tmp) != "try-error"){        
            expr2_list[[k2]] <- res
            k2 <- k2 + 1
        }
        tmp <- try(load(paste0(file_root, "_seed", seed, "_block(0.5)_fix.RData")))
        if (class(tmp) != "try-error"){        
            expr3_list[[k3]] <- res
            k3 <- k3 + 1
        }
    }

    expr1_res <- aggregate_expr(expr1_list, alphas) %>%
        postprocess %>%
        mutate(model = "AR (0.8)")
    expr2_res <- aggregate_expr(expr2_list, alphas) %>%
        postprocess %>%
        mutate(model = "AR (-0.8)")
    expr3_res <- aggregate_expr(expr3_list, alphas) %>%
        postprocess %>%
        mutate(model = "block (0.5)")
    aggregate_res[[i]] <- rbind(expr1_res,
                                expr2_res,
                                expr3_res) %>%
        mutate(n = settings$n[i],
               pi1 = settings$pi1[i],
               side = settings$side[i])
}

res <- do.call(rbind, aggregate_res)
save(res, file = "../data/dBH_mvgauss_aggregate.RData")

print(summary(res$qmax[res$qmax > -Inf]))
print(summary(res$q99[res$q99 > -Inf]))
print(summary(res$q95[res$q95 > -Inf]))
