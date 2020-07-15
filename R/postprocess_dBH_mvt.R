library("dplyr")
source("expr_functions.R")

aggregate_expr <- function(objlist){
    alphas <- sapply(objlist[[1]], function(y){
        y$alpha
    })
    res <- list()
    for (k in 1:length(alphas)){
        FDP <- do.call(cbind, lapply(objlist, function(x){
            x[[k]]$FDP
        }))
        power <- do.call(cbind, lapply(objlist, function(x){
            x[[k]]$power
        }))
        secBH <- do.call(cbind, lapply(objlist, function(x){
            x[[k]]$secBH
        }))
        res[[k]] <- list(alpha = alphas[k], FDP = FDP, power = power, secBH = secBH)
    }
    return(res)
}

params <- read.table("../jobs/dBH_mvt_params.txt")
names(params) <- c("seed", "n", "df", "pi1", "mutype", "side", "nreps", "dBH2")
settings <- params %>% select(-seed) %>% unique

aggregate_res <- list()
for (i in 1:nrow(settings)){
    seedlist <- params %>%
        filter(n == settings$n[i],
               df == settings$df[i],
               pi1 == settings$pi1[i],
               mutype == settings$mutype[i],
               side == settings$side[i],
               nreps == settings$nreps[i],
               dBH2 == settings$dBH2[i]) %>%
        .$seed
    file_root <- paste0("../cluster_raw_data/dBH_mvt",
                        "_n", settings$n[i],
                        "_df", settings$df[i],
                        "_pi1", settings$pi1[i],
                        "_mutype", settings$mutype[i],
                        "_side", settings$side[i],
                        "_nreps", settings$nreps[i],
                        "_dBH2", settings$dBH2[i])
    expr1_list <- list()
    expr2_list <- list()
    expr3_list <- list()
    k <- 1
    
    for (seed in seedlist){
        tmp <- try(load(paste0(file_root, "_seed", seed, "_AR(0.8)_fix.RData")))
        if (class(tmp) != "try-error"){        
            expr1_list[[k]] <- res
        }
        tmp <- try(load(paste0(file_root, "_seed", seed, "_iid(0)_fix.RData")))
        if (class(tmp) != "try-error"){        
            expr2_list[[k]] <- res
        }
        tmp <- try(load(paste0(file_root, "_seed", seed, "_block(0.5)_fix.RData")))
        if (class(tmp) != "try-error"){        
            expr3_list[[k]] <- res
        }
        k <- k + 1
    }

    expr1_res <- aggregate_expr(expr1_list) %>%
        postprocess %>%
        mutate(model = "AR (0.8)")
    expr2_res <- aggregate_expr(expr2_list) %>%
        postprocess %>%
        mutate(model = "iid")
    expr3_res <- aggregate_expr(expr3_list) %>%
        postprocess %>%
        mutate(model = "block (0.5)")
    aggregate_res[[i]] <- rbind(expr1_res,
                                expr2_res,
                                expr3_res) %>%
        mutate(n = settings$n[i],
               df = settings$df[i],
               pi1 = settings$pi1[i],
               side = settings$side[i])
}

res <- do.call(rbind, aggregate_res)
save(res, file = "../data/dBH_mvt_aggregate.RData")
