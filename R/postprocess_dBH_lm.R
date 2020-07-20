library("dplyr")
source("expr_functions.R")

params <- read.table("../jobs/dBH_lm_params.txt")
names(params) <- c("seed", "Xseed", "n", "p", "pi1", "mutype", "side", "nreps", "dBH2", "knockoff")
params <- params %>% filter(seed < 100)
settings <- params %>% select(-seed) %>% unique
alphas <- c(0.05, 0.2)

aggregate_res <- list()
for (i in 1:nrow(settings)){
    seedlist <- params %>%
        filter(Xseed == settings$Xseed[i],
               n == settings$n[i],
               p == settings$p[i],
               pi1 == settings$pi1[i],
               mutype == settings$mutype[i],
               side == settings$side[i],
               nreps == settings$nreps[i],
               dBH2 == settings$dBH2[i],
               knockoff == settings$knockoff[i]) %>%
        .$seed
    file_root <- paste0("../cluster_raw_data/dBH_lm",
                        "_n", settings$n[i],
                        "_p", settings$p[i],
                        "_pi1", settings$pi1[i],
                        "_mutype", settings$mutype[i],
                        "_side", settings$side[i],
                        "_nreps", settings$nreps[i],
                        "_dBH2", settings$dBH2[i],
                        "_knockoff", settings$knockoff[i])
    expr1_list <- list()
    expr2_list <- list()
    expr3_list <- list()
    k1 <- 1
    k2 <- 1
    k3 <- 1
    
    for (seed in seedlist){
        tmp <- try(load(paste0(file_root, "_seed", seed, "_Xseed", settings$Xseed[i], "_AR(0.8)_fix.RData")))
        if (class(tmp) != "try-error"){        
            expr1_list[[k1]] <- res
            k1 <- k1 + 1
        }
        tmp <- try(load(paste0(file_root, "_seed", seed, "_Xseed", settings$Xseed[i], "_iid(0)_fix.RData")))
        if (class(tmp) != "try-error"){        
            expr2_list[[k2]] <- res
            k2 <- k2 + 1
        }
        tmp <- try(load(paste0(file_root, "_seed", seed, "_Xseed", settings$Xseed[i], "_block(0.5)_fix.RData")))
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
        mutate(model = "iid")
    expr3_res <- aggregate_expr(expr3_list, alphas) %>%
        postprocess %>%
        mutate(model = "block (0.5)")
    aggregate_res[[i]] <- rbind(expr1_res,
                                expr2_res,
                                expr3_res) %>%
        mutate(n = settings$n[i],
               p = settings$p[i],
               pi1 = settings$pi1[i],
               side = settings$side[i])
}

res <- do.call(rbind, aggregate_res)
save(res, file = "../data/dBH_lm_aggregate.RData")

print(summary(res$qmax[res$qmax > -Inf]))
print(summary(res$q99[res$q99 > -Inf]))
print(summary(res$q95[res$q95 > -Inf]))
