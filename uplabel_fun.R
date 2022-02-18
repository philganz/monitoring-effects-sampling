uplabel_fun <- function(data, units, metrics, groups, yn_var, n_rep, bonferroni = FALSE, n_pop, uplabel_prop, boxcox = FALSE, lambda = seq(-10, 10, 0.01)){

# Start time ----------------------------------------------------------------------------
start_time <- Sys.time()
    
# Save function arguments for output-----------------------------------------------------
args <- data.table(n_rep, n_pop, uplabel_prop)
  
# Initial permutation test---------------------------------------------------------------
  
## Separate the permutation data from the input data and isolate the columns of interest
pdat <- data[, c(..units, ..groups, ..metrics, ..yn_var)]

## Print status report
message("Running initial permutation test to confirm evidence of observer effects")

## Run initial permutation test to confirm evidence of observer effects
res <- perm_fun(data = pdat, units, metrics, groups, yn_var, n_rep, bonferroni, boxcox, lambda)

## Identify strata that show any evidence of an observer effect and are therefore candidates for increased sampling
upstrata <- unique(melt(res$`p-value`, measure.vars = metrics, value.name = "p_val")[p_val < 0.05, ..groups])

## Get starting sample sizes for candidate strata
upstrata <- res$N_table[upstrata, on = c(groups)]

## Order groups by how the were entered into the function
setorderv(upstrata, groups)

## Create empty lists to hold permutation results
N_table <- list()
obs     <- list()
exp     <- list()
p_val   <- list()

## Create empty lists to hold permutation results for each stratum
N_table_strat <- list()
obs_strat     <- list()
exp_strat     <- list()
p_val_strat   <- list()

# Group-level permutation tests-------------------------------------------------------- 
for(h in 1:nrow(upstrata)){

  ## Start an increment count
  i <- 0
  
  ## Start counting additional trips
  t <- 0

  ## Isolate one stratum at a time
  sdat <- pdat[upstrata[h], on = c(groups)]

  ## Isolate observed and unobserved trips
  observed   <- as.vector(unlist(unique(sdat[get(yn_var) == "Y", c(..units)])))
  unobserved <- as.vector(unlist(unique(sdat[get(yn_var) == "N", c(..units)])))

  ## Start test value at initial permutation p-values
  test <- res$`p-value`[upstrata[h, ..groups], on = c(groups)]

  ## Create empty lists to hold permutation results for each increment
  N_table_inc <- list()
  obs_inc     <- list()
  exp_inc     <- list()
  p_val_inc   <- list()

# Increment-level permutation tests------------------------------------------------------
while((sum(melt(test[, ..metrics], measure.vars = metrics)$value < 0.05) != 0) & (t < length(unobserved))){
  
    ## Increase the increment count
    i <- i + 1
  
    ## Increase the number of additional trips to be observed
    t <- round(uplabel_prop * upstrata[h, n]) * i 
  
    ## Create empty lists to hold permutation results for each population
    N_table_pop <- list()
    obs_pop     <- list()
    exp_pop     <- list()
    p_val_pop   <- list()
    
    ## Print status report
    message(paste0(" Adding ", t, " trips to the ", paste(as.vector(unlist(upstrata[h, c(..groups)])), collapse = " "), " group "))
    
# Population-level permutation tests-----------------------------------------------------  
  for(p in 1:n_pop){
    
      ## Choose t unobserved trip(s) to relabel as observed
      updat <- copy(sdat)[get(units) %in% sample(unobserved, min(t, length(unobserved))), (yn_var) := "Y"]
    
      ## Run permutation test on this population of trips  
      res_pop <- perm_fun(data = updat, units, metrics, groups, yn_var, n_rep, bonferroni, boxcox, lambda)
    
# Save population-level results----------------------------------------------------------
      N_table_pop[[p]] <- res_pop$N_table[, POP := p]
      obs_pop[[p]]     <- res_pop$obs[, POP := p]
      exp_pop[[p]]     <- res_pop$exp[, POP := p]
      p_val_pop[[p]]   <- res_pop$`p-value`[, POP := p]
    
      }
  
# Save increment-level results-----------------------------------------------------------
    N_table_inc[[i]] <- rbindlist(N_table_pop)[, TRIPS_ADDED := t]
    obs_inc[[i]]     <- rbindlist(obs_pop)[, TRIPS_ADDED := t]
    exp_inc[[i]]     <- rbindlist(exp_pop)[, TRIPS_ADDED := t]
    p_val_inc[[i]]   <- rbindlist(p_val_pop)[, TRIPS_ADDED := t]
  
    ### Create object to test whether observer effects are still present (according to p-value)
    test <- unique(p_val_inc[[i]][, ..metrics][, lapply(.SD, function(x) median(x, na.rm = TRUE)), .SDcols = metrics])
  
    }

# Save stratum-level permutation results-------------------------------------------------
  N_table_strat[[h]] <- rbindlist(N_table_inc)
  obs_strat[[h]]     <- rbindlist(obs_inc)
  exp_strat[[h]]     <- rbindlist(exp_inc)
  p_val_strat[[h]]   <- rbindlist(p_val_inc)

  }

### Save all results---------------------------------------------------------------------
N_table <- rbindlist(N_table_strat)
obs     <- rbindlist(obs_strat)
exp     <- rbindlist(exp_strat)
p_val   <- rbindlist(p_val_strat)

# End time ------------------------------------------------------------------------------
end_time <- Sys.time()

# Return results ------------------------------------------------------------------------
return(list("run_time" = end_time - start_time, "args" = args, "N_table" = N_table, "obs" = obs, "exp" = exp, "p-value" = p_val))

}