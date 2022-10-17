# Set random number seed --------------------------------------------------
set.seed(49)

# Get packages ------------------------------------------------------------
library(dplyr)
library(data.table)
library(flextable)
library(getPass)
library(ggplot2)
library(ggridges)
library(officer)
library(ROracle)

# Get data ----------------------------------------------------------------
if(!file.exists("Observer_Effects/data/work.data.rds")){
  source("Observer_Effects/get_data.r")
  
} else{
  work.data <- readRDS("Observer_Effects/data/work.data.rds")
  setDT(work.data)
}

# Get user-defined functions ----------------------------------------------
source("Observer_Effects/perm_fun.R")
source("Observer_Effects/uplabel_fun.R")

# Wrangle data for permutation test ---------------------------------------

# Select rows of interest
perm.data <- work.data[!(STRATA %in% c("FULL", "ZERO", "EM_POT"))]

# Remove jig trips
jig       <- unique(perm.data[AGENCY_GEAR_CODE == "JIG", TRIP_ID])
perm.data <- perm.data[!(TRIP_ID %in% jig)]

# Reformat EM_HAL stratum name
perm.data[STRATA == "EM_HAL", STRATA := "EM HAL"]

# Create post-strata that include tender
perm.data[, P_STRATA := ifelse(STRATA == "POT" & TENDER == "Y", "TenP", STRATA)]
perm.data[, P_STRATA := ifelse(STRATA == "TRW" & TENDER == "Y", "TenTR", P_STRATA)]

# Create gear type column that also reflects EM status
perm.data[, EM_GEAR := ifelse(STRATA == "EM HAL", "EM HAL", AGENCY_GEAR_CODE)]

# Calculate permutation test metrics --------------------------------------

# Number of areas fished
perm.data[, N_AREAS := as.numeric(length(unique(REPORTING_AREA_CODE))), by=.(TRIP_ID, EM_GEAR)]

# Trip duration in days (removes 6 trips that had landing date before trip start)
perm.data <- perm.data[, MIN_DATE := min(as.Date(TRIP_TARGET_DATE)), by = .(TRIP_ID, EM_GEAR)
                       ][, MAX_DATE := max(as.Date(LANDING_DATE)), by = .(TRIP_ID, EM_GEAR)
                         ][, DAYS := as.numeric(difftime(MAX_DATE, MIN_DATE, units = "days") + 1), by = .(TRIP_ID, EM_GEAR)
                           ][, c("TRIP_TARGET_DATE", "LANDING_DATE", "MIN_DATE", "MAX_DATE") := NULL
                             ][DAYS > 0]
# Length overall  
perm.data[, LENGTH_OVERALL := as.numeric(LENGTH_OVERALL)]

# Number of species
perm.data[, N_SPECIES := as.numeric(length(unique(AGENCY_SPECIES_CODE[SOURCE_TABLE == "Y"]))), by = .(TRIP_ID, EM_GEAR)]

# Landed catch
perm.data[, LANDED_CATCH := sum(as.numeric(WEIGHT_POSTED[SOURCE_TABLE == "Y"])), by = .(TRIP_ID, EM_GEAR)]

# Proportion of catch that is predominant species (pMax)
perm.data[, SPECIES_WEIGHT := sum(WEIGHT_POSTED[SOURCE_TABLE == "Y"], na.rm = TRUE), by = .(TRIP_ID, EM_GEAR, AGENCY_SPECIES_CODE)
          ][, PMAX := max(SPECIES_WEIGHT/LANDED_CATCH), by = .(TRIP_ID, EM_GEAR)
            ][, SPECIES_WEIGHT := NULL]

# Isolate columns of interest and make rows unique
perm.data <- unique(perm.data[, .(ADP, VESSEL_ID, TRIP_ID, COVERAGE_TYPE, EM_GEAR, TENDER, OBSERVED_FLAG, N_AREAS, DAYS, LENGTH_OVERALL, N_SPECIES, PMAX, LANDED_CATCH)])

# Define variables of interest --------------------------------------------
units        <- "TRIP_ID"
metrics      <- c("N_AREAS", "DAYS", "LENGTH_OVERALL", "N_SPECIES", "PMAX", "LANDED_CATCH")
groups       <- c("ADP", "EM_GEAR")
yn_var       <- "OBSERVED_FLAG"
n_rep        <- 1000
n_pop        <- 100
uplabel_prop <- 0.2

# Permutation tests -------------------------------------------------------

# With transformed data 
if(file.exists("Observer_Effects/data/perm_results_t.rdata")){
  load("Observer_Effects/data/perm_results_t.rdata")
  
} else{
  
  perm_results_t <- perm_fun(data = perm.data, units, metrics, groups, yn_var, n_rep, boxcox = TRUE)
  
  # Save results
  save(perm_results_t, file = "Observer_Effects/data/perm_results_t.rdata")
}

# With untransformed data
if(file.exists("Observer_Effects/data/perm_results_u.rdata")){
  load("Observer_Effects/data/perm_results_u.rdata")
  
} else{
  
perm_results_u <- perm_fun(data = perm.data, units, metrics, groups, yn_var, n_rep, boxcox = FALSE)
                         
# Save results
save(perm_results_u, file = "Observer_Effects/data/perm_results_u.rdata")
}

# Increased monitoring exercise -----------------------------------------------------

# Sample to 100% with uplabel_prop = 0.2

# With transformed data, uplabeling (flip) method
if(file.exists("Observer_Effects/data/uplabel_results_tf.rdata")){
  load("Observer_Effects/data/uplabel_results_tf.rdata")
  
} else{
  
  uplabel_results_tf <- uplabel_fun(perm.data, units, metrics, groups, yn_var, n_rep, n_pop, uplabel_prop, method = "flip", boxcox = TRUE)
  
  # Save results
  save(uplabel_results_tf, file = "Observer_Effects/data/uplabel_results_tf.rdata")
}

# With transformed data, resampling (swap) method  
if(file.exists("Observer_Effects/data/uplabel_results_ts.rdata")){
  load("Observer_Effects/data/uplabel_results_ts.rdata")
  
} else{
  
  uplabel_results_ts <- uplabel_fun(perm.data, units, metrics, groups, yn_var, n_rep, n_pop, uplabel_prop, method  = "swap", boxcox = TRUE)
  
  # Save results
  save(uplabel_results_ts, file = "Observer_Effects/data/uplabel_results_ts.rdata")
}

# With untransformed data, uplabeling (flip) method  
if(file.exists("Observer_Effects/data/uplabel_results_uf.rdata")){
  load("Observer_Effects/data/uplabel_results_uf.rdata")
  
} else{
  
  uplabel_results_uf <- uplabel_fun(perm.data, units, metrics, groups, yn_var, n_rep, n_pop, uplabel_prop, method = "flip", boxcox = FALSE)

  # Save results
  save(uplabel_results_uf, file = "Observer_Effects/data/uplabel_results_uf.rdata")
}

# With untransformed data, resampling (swap) method  
if(file.exists("Observer_Effects/data/uplabel_results_us.rdata")){
  load("Observer_Effects/data/uplabel_results_us.rdata")
  
} else{
  
  uplabel_results_us <- uplabel_fun(perm.data, units, metrics, groups, yn_var, n_rep, n_pop, uplabel_prop, method  = "swap", boxcox = FALSE)
  
  # Save results
  save(uplabel_results_us, file = "Observer_Effects/data/uplabel_results_us.rdata")
}

# Sample to median p-value > 0.05/length(metrics) with uplabel_prop = 0.05
set.seed(49)
uplabel_prop <- 0.05

# With transformed data, uplabeling (flip) method
if(file.exists("Observer_Effects/data/uplabel_results_tf_p_stop.rdata")){
  load("Observer_Effects/data/uplabel_results_tf_p_stop.rdata")
  
} else{
  
  uplabel_results_tf_p_stop <- uplabel_fun(perm.data, units, metrics, groups, yn_var, n_rep, n_pop, uplabel_prop, method = "flip", boxcox = TRUE, p_stop = TRUE)
  
  # Save results
  save(uplabel_results_tf_p_stop, file = "Observer_Effects/data/uplabel_results_tf_p_stop.rdata")
}

# With transformed data, resampling (swap) method  
if(file.exists("Observer_Effects/data/uplabel_results_ts_p_stop.rdata")){
  load("Observer_Effects/data/uplabel_results_ts_p_stop.rdata")
  
} else{
  
  uplabel_results_ts_p_stop <- uplabel_fun(perm.data, units, metrics, groups, yn_var, n_rep, n_pop, uplabel_prop, method  = "swap", boxcox = TRUE, p_stop = TRUE)
  
  # Save results
  save(uplabel_results_ts_p_stop, file = "Observer_Effects/data/uplabel_results_ts_p_stop.rdata")
}

# With untransformed data, uplabeling (flip) method  
if(file.exists("Observer_Effects/data/uplabel_results_uf_p_stop.rdata")){
  load("Observer_Effects/data/uplabel_results_uf_p_stop.rdata")
  
} else{
  
  uplabel_results_uf_p_stop <- uplabel_fun(perm.data, units, metrics, groups, yn_var, n_rep, n_pop, uplabel_prop, method = "flip", boxcox = FALSE, p_stop = TRUE)
  
  # Save results
  save(uplabel_results_uf_p_stop, file = "Observer_Effects/data/uplabel_results_uf_p_stop.rdata")
}

# With untransformed data, resampling (swap) method  
if(file.exists("Observer_Effects/data/uplabel_results_us_p_stop.rdata")){
  load("Observer_Effects/data/uplabel_results_us_p_stop.rdata")
  
} else{
  
  uplabel_results_us_p_stop <- uplabel_fun(perm.data, units, metrics, groups, yn_var, n_rep, n_pop, uplabel_prop, method  = "swap", boxcox = FALSE, p_stop = TRUE)
  
  # Save results
  save(uplabel_results_us_p_stop, file = "Observer_Effects/data/uplabel_results_us_p_stop.rdata")
}


# Tables ------------------------------------------------------------------

# Table 1: percent differences at the simulated monitoring rates needed to bring the median p-value (across populations) above 0.05/length(metrics) for all metrics

# Get the median p-values that resulted for each combination of year, gear type, metric, increased monitoring method, transformation status, and sampling rate
dat <- rbind(copy(uplabel_results_uf_p_stop$`p-value`)[, ':=' (BOXCOX = uplabel_results_uf_p_stop$args$boxcox, METHOD = uplabel_results_uf_p_stop$args$method, VARIABLE = "PVAL")], 
             copy(uplabel_results_tf_p_stop$`p-value`)[, ':=' (BOXCOX = uplabel_results_tf_p_stop$args$boxcox, METHOD = uplabel_results_tf_p_stop$args$method, VARIABLE = "PVAL")],
             copy(uplabel_results_us_p_stop$`p-value`)[, ':=' (BOXCOX = uplabel_results_us_p_stop$args$boxcox, METHOD = uplabel_results_us_p_stop$args$method, VARIABLE = "PVAL")],
             copy(uplabel_results_ts_p_stop$`p-value`)[, ':=' (BOXCOX = uplabel_results_ts_p_stop$args$boxcox, METHOD = uplabel_results_ts_p_stop$args$method, VARIABLE = "PVAL")],
             copy(uplabel_results_uf_p_stop$obs[SCALE == "PCT"])[, ':=' (SCALE = NULL, BOXCOX = uplabel_results_uf_p_stop$args$boxcox, METHOD = uplabel_results_uf_p_stop$args$method, VARIABLE = "PCT_DIFF")], 
             copy(uplabel_results_tf_p_stop$obs[SCALE == "PCT"])[, ':=' (SCALE = NULL, BOXCOX = uplabel_results_tf_p_stop$args$boxcox, METHOD = uplabel_results_tf_p_stop$args$method, VARIABLE = "PCT_DIFF")],
             copy(uplabel_results_us_p_stop$obs[SCALE == "PCT"])[, ':=' (SCALE = NULL, BOXCOX = uplabel_results_us_p_stop$args$boxcox, METHOD = uplabel_results_us_p_stop$args$method, VARIABLE = "PCT_DIFF")],
             copy(uplabel_results_ts_p_stop$obs[SCALE == "PCT"])[, ':=' (SCALE = NULL, BOXCOX = uplabel_results_ts_p_stop$args$boxcox, METHOD = uplabel_results_ts_p_stop$args$method, VARIABLE = "PCT_DIFF")])
dat <- dat[, lapply(.SD, function(x) median(x)), .SDcols = metrics, keyby = c(groups, "BOXCOX", "METHOD", "TRIPS_ADDED", "VARIABLE")]
dat <- melt(dat, measure.vars = metrics, variable.name = "METRICS", value.name = "MEDIAN")

# Add sample size and rates achieved at each given number of trips added  
dat <- unique(
  rbind(
    uplabel_results_uf_p_stop$N_table[, .(ADP, EM_GEAR, BOXCOX = uplabel_results_uf_p_stop$args$boxcox, METHOD = uplabel_results_uf_p_stop$args$method, TRIPS_ADDED, N, n, r = rate)],
    uplabel_results_tf_p_stop$N_table[, .(ADP, EM_GEAR, BOXCOX = uplabel_results_tf_p_stop$args$boxcox, METHOD = uplabel_results_tf_p_stop$args$method, TRIPS_ADDED, N, n, r = rate)],
    uplabel_results_us_p_stop$N_table[, .(ADP, EM_GEAR, BOXCOX = uplabel_results_us_p_stop$args$boxcox, METHOD = uplabel_results_us_p_stop$args$method, TRIPS_ADDED, N, n, r = rate)],
    uplabel_results_ts_p_stop$N_table[, .(ADP, EM_GEAR, BOXCOX = uplabel_results_ts_p_stop$args$boxcox, METHOD = uplabel_results_ts_p_stop$args$method, TRIPS_ADDED, N, n, r = rate)]
  )
)[dat, on = .(ADP, EM_GEAR, BOXCOX, METHOD, TRIPS_ADDED)]

# Replace NAs at 100% sampling
dat[N == n & VARIABLE == "PVAL", MEDIAN := 1]
dat[N == n & VARIABLE == "PCT_DIFF", MEDIAN := 0]

# Cast to wide format in order to get p-value and percent difference on the same row
dat <- dcast(dat, ADP + EM_GEAR + BOXCOX + METHOD + TRIPS_ADDED + N + n + r + METRICS ~ paste0("MEDIAN_", VARIABLE), value.var = "MEDIAN")

# Add column for whether p <= 0.05/length(metrics)
dat[, SIG := MEDIAN_PVAL <= 0.05/length(metrics)]

# Isolate the lowest sampling rate at which all metrics have p > 0.05/length(metrics)
pass <- unique(dat[, .SD[sum(SIG) == 0], by = .(ADP, EM_GEAR, BOXCOX, METHOD, TRIPS_ADDED)
][, .SD[TRIPS_ADDED == min(TRIPS_ADDED)], by = .(ADP, EM_GEAR, BOXCOX, METHOD)]
[, .(ADP, EM_GEAR, BOXCOX, METHOD, TRIPS_ADDED)])
dat <- dat[pass, on = .(ADP, EM_GEAR, BOXCOX, METHOD, TRIPS_ADDED)]

# Isolate the metric with the lowest p-value of those that passed
dat <- dat[, .SD[MEDIAN_PVAL == min(MEDIAN_PVAL)], by = .(ADP, EM_GEAR, BOXCOX, METHOD)][, ':='(MEDIAN_PVAL = NULL, SIG = NULL)]

# Save these "lagging" metrics for later use with figures
lagging_metrics <- unique(dat[, .(ADP, EM_GEAR, BOXCOX, METHOD, METRICS)])

# Add original sampling rates and percent differences
dat <- rbind(
  merge(perm_results_u$N_table[, .(EM_GEAR, ADP, N_o = N, n_o = n, r_o = round(rate * 100, 2))], on = .(ADP, EM_GEAR, METRICS),
        melt(perm_results_u$obs[SCALE == "PCT"], measure.vars = metrics, variable.name = "METRICS", value.name = "d_o")
  )[, ":=" (SCALE = NULL, BOXCOX = FALSE)],
  merge(perm_results_t$N_table[, .(EM_GEAR, ADP, N_o = N, n_o = n, r_o = round(rate * 100, 2))], on = .(ADP, EM_GEAR, METRICS),
        melt(perm_results_u$obs[SCALE == "PCT"], measure.vars = metrics, variable.name = "METRICS", value.name = "d_o")
  )[, ":=" (SCALE = NULL, BOXCOX = TRUE)]
)[dat, on = .(ADP, EM_GEAR, BOXCOX, METRICS)]

# Express rates as percentages
dat[, r := r * 100]

# Round values
dat[, c("d_o", "r", "MEDIAN_PCT_DIFF") := lapply(.(d_o, r, MEDIAN_PCT_DIFF), round, 2)]

# Recode metrics
dat[, METRICS := recode(METRICS, "N_AREAS" = "Areas", "DAYS" = "Days", "LENGTH_OVERALL" = "Length (ft)", "N_SPECIES" = "Species", "PMAX" = "pMax", "LANDED_CATCH" = "Landings (t)")]

# Start creating wide table
dat[, TREATMENT := paste(BOXCOX, METHOD, sep = "_")]
dat[, TREATMENT := recode(TREATMENT, "FALSE_flip" = "uu", "FALSE_swap" = "ur", "TRUE_flip" = "tu", "TRUE_swap" = "tr")]
dat <- dcast(dat, EM_GEAR + ADP + d_o + METRICS ~ paste0("n_", TREATMENT), value.var = "n")[
  dcast(dat, EM_GEAR + ADP + d_o + METRICS ~ paste0("r_", TREATMENT), value.var = "r")][
    dcast(dat, EM_GEAR + ADP + d_o + METRICS ~ paste0("d_", TREATMENT), value.var = "MEDIAN_PCT_DIFF")]

# Add original sample rates
dat <- dat[perm_results_u$N_table[, .(EM_GEAR, ADP, N, n, r = round(rate * 100, 2))], on = .(EM_GEAR, ADP)][order(EM_GEAR, ADP)]

# Where no increased sampling was needed, replace sample sizes and rates with originals
dat[is.na(METRICS), ':=' (n_uu = n, r_uu = r)]
dat[is.na(METRICS), ':=' (n_ur = n, r_ur = r)]
dat[is.na(METRICS), ':=' (n_tu = n, r_tu = r)]
dat[is.na(METRICS), ':=' (n_tr = n, r_tr = r)]

dat[, c("n_uu", "n_ur", "n_tu", "n_tr") := lapply(.(n_uu, n_ur, n_tu, n_tr), function(x) replace(x, is.na(x), 0))]

# Add totals
totals <- dat[, lapply(.SD, function(x) sum(x, na.rm = TRUE)), .SDcols = c("n_uu", "n_ur", "n_tu", "n_tr")]
totals[, ':=' (EM_GEAR = "Total", ADP = "", N = sum(unique(dat[, .(EM_GEAR, ADP, N)])$N), n = sum(unique(dat[, .(EM_GEAR, ADP, n)])$n))]
totals[, ':=' (r = round(n / N * 100, 2), d_o = "", METRICS = "",
               d_uu = "", r_uu = round(n_uu / N * 100, 2), 
               d_ur = "", r_ur = round(n_ur / N * 100, 2), 
               d_tu = "", r_tu = round(n_tu / N * 100, 2), 
               d_tr = "", r_tr = round(n_tr / N * 100, 2))]

dat <- rbind(dat, totals)

# Reorder columns and rows
dat <- dat[, .(Gear = EM_GEAR, Year = ADP, N, r, d = d_o, m = METRICS, r_uu, d_uu, r_ur, d_ur, r_tu, d_tu, r_tr, d_tr)][order(Gear, Year)]

# Format columns so that digits to the right of the decimal are retained
dat <- dat[, lapply(.SD, function(x) format(x, nsmall = 2, big.mark = ","))]

# Replace NAs with dashes
dat <- dat[, lapply(.SD, function(x) replace(x, grepl("NA", x), "-"))]

# Create Table 1 header
t1_header <- data.frame(col_keys  = c("Gear", "Year", "N", "r", "m", "r_uu", "d_uu", "r_ur", "d_ur", "r_tu", "d_tu", "r_tr", "d_tr"),
                        line1 = c(rep("", 2), rep("Original", 2), rep("Increased monitoring results", 9)),
                        line2 = c(rep("", 5), rep("Using untransformed data", 4), rep("Using transformed data", 4)),
                        line3  = c("Gear", "Year", "N", "r", "m", "ru", "du", "rr", "dr", "ru", "du", "rr", "dr"),
                        stringsAsFactors = FALSE)

# Create Table 1
t1 <- dat %>% 
      flextable(col_keys = t1_header$col_keys) %>%
      set_header_df(mapping = t1_header, key = "col_keys") %>% 
      theme_booktabs() %>% 
      merge_h(part = "header") %>%
      merge_v(j = "Gear") %>% 
      hline(i = c(2, 5, 8, 11, 15), border = fp_border(color = "black", style = "solid", width = 1))

# Save Table 1
# save_as_docx(t1, path = "Observer_Effects/tables/t1.docx")

# Figures -----------------------------------------------------------------

## Figure 2: median p-values against observation rate for lagging metrics
dat <- uplabel_results_uf$N_table[uplabel_results_uf$`p-value`, on = .(ADP, EM_GEAR, POP, TRIPS_ADDED)][, POP := NULL][, ':=' (BOXCOX = uplabel_results_uf$args$boxcox, METHOD = uplabel_results_uf$args$method)]
dat <- rbind(dat, uplabel_results_tf$N_table[uplabel_results_tf$`p-value`, on = .(ADP, EM_GEAR, POP, TRIPS_ADDED)][, POP := NULL][, ':=' (BOXCOX = uplabel_results_tf$args$boxcox, METHOD = uplabel_results_tf$args$method)])
dat <- rbind(dat, uplabel_results_us$N_table[uplabel_results_us$`p-value`, on = .(ADP, EM_GEAR, POP, TRIPS_ADDED)][, POP := NULL][, ':=' (BOXCOX = uplabel_results_us$args$boxcox, METHOD = uplabel_results_us$args$method)])
dat <- rbind(dat, uplabel_results_ts$N_table[uplabel_results_ts$`p-value`, on = .(ADP, EM_GEAR, POP, TRIPS_ADDED)][, POP := NULL][, ':=' (BOXCOX = uplabel_results_ts$args$boxcox, METHOD = uplabel_results_ts$args$method)])
dat <- rbind(dat, perm_results_u$N_table[perm_results_u$`p-value`, on = .(ADP, EM_GEAR)][, ":=" (TRIPS_ADDED = 0, BOXCOX = FALSE, METHOD = "flip")])
dat <- rbind(dat, perm_results_u$N_table[perm_results_u$`p-value`, on = .(ADP, EM_GEAR)][, ":=" (TRIPS_ADDED = 0, BOXCOX = FALSE, METHOD = "swap")])
dat <- rbind(dat, perm_results_t$N_table[perm_results_t$`p-value`, on = .(ADP, EM_GEAR)][, ":=" (TRIPS_ADDED = 0, BOXCOX = TRUE, METHOD = "flip")])
dat <- rbind(dat, perm_results_t$N_table[perm_results_t$`p-value`, on = .(ADP, EM_GEAR)][, ":=" (TRIPS_ADDED = 0, BOXCOX = TRUE, METHOD = "swap")])
dat <- melt(dat, measure.vars = metrics, variable.name = "METRICS")
dat <- dat[lagging_metrics, on = .(ADP, EM_GEAR, BOXCOX, METHOD, METRICS)]
dat[, EM_GEAR := factor(EM_GEAR, levels = c("EM HAL", "HAL", "POT", "NPT", "PTR"))]
dat[, BOXCOX := ifelse(BOXCOX, "Yes", "No")]
dat[, BOXCOX := factor(BOXCOX, levels = c("Yes", "No"))]
dat[, METHOD := recode(METHOD, "flip" = "Uplabel", "swap" = "Resample")]
dat[, METHOD := factor(METHOD, levels = c("Uplabel", "Resample"))]
dat[, METRICS := as.character(recode(METRICS, "N_AREAS" = "Areas", "DAYS" = "Days", "LENGTH_OVERALL" = "Length (ft)", "N_SPECIES" = "Species", "PMAX" = "pMax", "LANDED_CATCH" = "Landings (t)"))]

# Replace NAs at 100% sampling with 1, since 100% of permutations are at least as extreme as the observed difference (the two are equal)
dat[N == n, value := 1]

# Calculate median and confidence interval  
dat <- dat[, .(lci = quantile(value, 0.25), median = median(value), uci = quantile(value, 0.75)), keyby = .(ADP, EM_GEAR, N, n, rate, TRIPS_ADDED, BOXCOX, METHOD, METRICS)]

# Create Figure 2
f2 <- ggplot(dat) +
      geom_line(aes(x = rate * 100, y = median, shape = METHOD, color = BOXCOX)) +
      geom_point(aes(x = rate * 100, y = median, shape = METHOD, color = BOXCOX), size = 2) + 
      geom_ribbon(aes(x = rate * 100, ymin = lci, ymax = uci, shape = METHOD, fill = BOXCOX), alpha = 0.2) +
      geom_hline(yintercept = 0.05/length(metrics), linetype = 2, color = "red") +
      facet_wrap(c("EM_GEAR", "ADP", "METRICS"), scales = "free") +
      labs(x = "Simulated monitoring rate (%)", y = "p-value", fill = "Transformed", shape = "Method") +
      theme_bw() +
      scale_fill_manual(values = c("red", "black")) +
      scale_color_manual(values = c("red", "black")) +
      theme(text = element_text(size = 16),
            legend.position = "bottom") +
      guides(color = "none")

# Save Figure 2
# png("Observer_Effects/figures/f2.png", width = 11, height = 8.5, units = 'in', res = 300)
# f2
# dev.off()

## Figure 3: median difference (%) against observation for lagging metrics 
dat <- uplabel_results_uf$N_table[uplabel_results_uf$obs[SCALE == "PCT"], on = .(ADP, EM_GEAR, POP, TRIPS_ADDED)][, POP := NULL][, ':=' (BOXCOX = uplabel_results_uf$args$boxcox, METHOD = uplabel_results_uf$args$method)]
dat <- rbind(dat, uplabel_results_tf$N_table[uplabel_results_tf$obs[SCALE == "PCT"], on = .(ADP, EM_GEAR, POP, TRIPS_ADDED)][, POP := NULL][, ':=' (BOXCOX = uplabel_results_tf$args$boxcox, METHOD = uplabel_results_tf$args$method)])
dat <- rbind(dat, uplabel_results_us$N_table[uplabel_results_us$obs[SCALE == "PCT"], on = .(ADP, EM_GEAR, POP, TRIPS_ADDED)][, POP := NULL][, ':=' (BOXCOX = uplabel_results_us$args$boxcox, METHOD = uplabel_results_us$args$method)])
dat <- rbind(dat, uplabel_results_ts$N_table[uplabel_results_ts$obs[SCALE == "PCT"], on = .(ADP, EM_GEAR, POP, TRIPS_ADDED)][, POP := NULL][, ':=' (BOXCOX = uplabel_results_ts$args$boxcox, METHOD = uplabel_results_ts$args$method)])
dat <- rbind(dat, perm_results_u$N_table[perm_results_u$obs[SCALE == "PCT"], on = .(ADP, EM_GEAR)][, ":=" (TRIPS_ADDED = 0, BOXCOX = FALSE, METHOD = "flip")])
dat <- rbind(dat, perm_results_u$N_table[perm_results_u$obs[SCALE == "PCT"], on = .(ADP, EM_GEAR)][, ":=" (TRIPS_ADDED = 0, BOXCOX = FALSE, METHOD = "swap")])
dat <- rbind(dat, perm_results_t$N_table[perm_results_t$obs[SCALE == "PCT"], on = .(ADP, EM_GEAR)][, ":=" (TRIPS_ADDED = 0, BOXCOX = TRUE, METHOD = "flip")])
dat <- rbind(dat, perm_results_t$N_table[perm_results_t$obs[SCALE == "PCT"], on = .(ADP, EM_GEAR)][, ":=" (TRIPS_ADDED = 0, BOXCOX = TRUE, METHOD = "swap")])
dat <- melt(dat, measure.vars = metrics, variable.name = "METRICS")
dat <- dat[lagging_metrics, on = .(ADP, EM_GEAR, BOXCOX, METHOD, METRICS)]
dat[, EM_GEAR := factor(EM_GEAR, levels = c("EM HAL", "HAL", "POT", "NPT", "PTR"))]
dat[, BOXCOX := ifelse(BOXCOX, "Yes", "No")]
dat[, BOXCOX := factor(BOXCOX, levels = c("Yes", "No"))]
dat[, METHOD := recode(METHOD, "flip" = "Uplabel", "swap" = "Resample")]
dat[, METHOD := factor(METHOD, levels = c("Uplabel", "Resample"))]
dat[, METRICS := as.character(recode(METRICS, "N_AREAS" = "Areas", "DAYS" = "Days", "LENGTH_OVERALL" = "Length (ft)", "N_SPECIES" = "Species", "PMAX" = "pMax", "LANDED_CATCH" = "Landings (t)"))]

# Replace NAs at 100% sampling with 0, since there is no difference between monitored and unmonitored trips (all trips are monitored)
dat[N == n, value := 0]

# Calculate median and confidence interval  
dat <- dat[, .(lci = quantile(value, 0.25), median = median(value), uci = quantile(value, 0.75)), keyby = .(ADP, EM_GEAR, N, n, rate, TRIPS_ADDED, BOXCOX, METHOD, METRICS)]

# Create Figure 3
f3 <- ggplot(dat) +
      geom_line(aes(x = rate * 100, y = median, shape = METHOD, color = BOXCOX)) +
      geom_point(aes(x = rate * 100, y = median, shape = METHOD, color = BOXCOX), size = 2) + 
      geom_ribbon(aes(x = rate * 100, ymin = lci, ymax = uci, shape = METHOD, fill = BOXCOX), alpha = 0.2) +
      geom_hline(yintercept = 0.05/length(metrics), linetype = 2, color = "red") +
      facet_wrap(c("EM_GEAR", "ADP", "METRICS"), scales = "free") +
      labs(x = "Simulated monitoring rate (%)", y = "Observed difference (%)", fill = "Transformed", shape = "Method") +
      theme_bw() +
      scale_fill_manual(values = c("red", "black")) +
      scale_color_manual(values = c("red", "black")) +
      theme(text = element_text(size = 16),
            legend.position = "bottom") +
      guides(color = "none")

# Save Figure 3
# png("Observer_Effects/figures/f3.png", width = 11, height = 8.5, units = 'in', res = 300)
# f3
# dev.off()
