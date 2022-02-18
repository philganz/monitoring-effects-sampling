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
uplabel_prop <- 0.05

# Permutation tests -------------------------------------------------------

# Without data transformation  
if(file.exists("Observer_Effects/data/perm_results.rdata")){
  load("Observer_Effects/data/perm_results.rdata")
  
} else{
  
perm_results <- perm_fun(data = perm.data, units, metrics, groups, yn_var, n_rep, bonferroni = TRUE)
                         
# Save results
save(perm_results, file = "Observer_Effects/data/perm_results.rdata")
}

# With data transformation  
if(file.exists("Observer_Effects/data/perm_results_trans.rdata")){
  load("Observer_Effects/data/perm_results_trans.rdata")
  
} else{
  
perm_results_trans <- perm_fun(data = perm.data, units, metrics, groups, yn_var, n_rep, bonferroni = TRUE, boxcox = TRUE)
  
# Save results
save(perm_results_trans, file = "Observer_Effects/data/perm_results_trans.rdata")
}

# Uplabeling exercise -----------------------------------------------------

# Without data transformation  
if(file.exists("Observer_Effects/data/uplabel_results.rdata")){
  load("Observer_Effects/data/uplabel_results.rdata")
  
} else{
  
uplabel_results <- uplabel_fun(perm.data, units, metrics, groups, yn_var, n_rep, bonferroni = TRUE, n_pop, uplabel_prop)

# Save results
save(uplabel_results, file = "Observer_Effects/data/uplabel_results.rdata")
}

# With data transformation  
if(file.exists("Observer_Effects/data/uplabel_results_trans.rdata")){
  load("Observer_Effects/data/uplabel_results_trans.rdata")
  
} else{
  
uplabel_results_trans <- uplabel_fun(perm.data, units, metrics, groups, yn_var, n_rep, bonferroni = TRUE, n_pop, uplabel_prop, boxcox = TRUE)
  
# Save results
save(uplabel_results_trans, file = "Observer_Effects/data/uplabel_results_trans.rdata")
}

# Tables ------------------------------------------------------------------

# Table 1: simulated monitoring rates needed to reduce monitoring effects according to uplabeling exercise

# Get original sampling rates
dat <- perm_results$N_table[, .(EM_GEAR, ADP, N, n, r = round(rate * 100, 2))][, .SD, keyby = .(EM_GEAR, ADP)]

# Append uplabeling rates
dat <- unique(uplabel_results$N_table[, .SD[n == max(n)], keyby = .(EM_GEAR, ADP)][, .(EM_GEAR, ADP, n_u = n, r_u = round(rate * 100, 2))])[dat]
dat <- unique(uplabel_results_trans$N_table[, .SD[n == max(n)], keyby = .(EM_GEAR, ADP)][, .(EM_GEAR, ADP, n_t = n, r_t = round(rate * 100, 2))])[dat]

# Use the original sample size and rate for groups that didn't require uplabeling
dat[is.na(n_u), ':=' (n_u = n, r_u = r)]
dat[is.na(n_t), ':=' (n_t = n, r_t = r)]

# Get the median p-values that resulted for each combination of year, gear type, metric, transformation status, and sampling rate
all_metrics <- rbind(copy(uplabel_results$`p-value`)[, BOXCOX := "N_BOXCOX"], copy(uplabel_results_trans$`p-value`)[, BOXCOX := "Y_BOXCOX"])
all_metrics <- all_metrics[, lapply(.SD, function(x) median(x)), .SDcols = metrics, keyby = c(groups, "BOXCOX", "TRIPS_ADDED")]
all_metrics <- melt(all_metrics, measure.vars = metrics, variable.name = "METRICS", value.name = "MEDIAN_P_VALUE")

# Isolate the highest two sampling rates for each combination of year, gear type, metric, and transformation status
lagging_metrics <- all_metrics[order(-TRIPS_ADDED), .SD[1:2], keyby = c(groups, "BOXCOX", "METRICS")]
lagging_metrics <- lagging_metrics[!is.na(MEDIAN_P_VALUE)]

# Remove combinations that achieved a p-value above 0.05 before the highest sampling rate
passing_metrics <- unique(lagging_metrics[, .SD[MEDIAN_P_VALUE >= 0.05 & TRIPS_ADDED < max(TRIPS_ADDED)], keyby = c(groups, "BOXCOX")][, .(ADP, EM_GEAR, BOXCOX, METRICS)])
lagging_metrics <- lagging_metrics[!passing_metrics, on = .(ADP, EM_GEAR, BOXCOX, METRICS)][order(ADP, EM_GEAR, BOXCOX, METRICS)]

# From the remaining metrics, isolate those that had the lowest p-value at the highest sampling rate  
lagging_metrics <- lagging_metrics[, .SD[TRIPS_ADDED == max(TRIPS_ADDED)], keyby = c(groups, "BOXCOX")]
lagging_metrics <- lagging_metrics[, .SD[MEDIAN_P_VALUE == min(MEDIAN_P_VALUE)], keyby = c(groups, "BOXCOX")]

# Reformat the lagging metrics in preparation for joining  
lagging_metrics <- unique(lagging_metrics[, .(ADP, EM_GEAR, BOXCOX, METRICS = recode(METRICS, "N_AREAS" = "NMFS areas", "DAYS" = "Days fished", "LENGTH_OVERALL" = "Vessel length (ft)", "N_SPECIES" = "Species landed", "PMAX" = "pMax species", "LANDED_CATCH" = "Landed catch (t)"))])
lagging_metrics <- dcast(lagging_metrics, ... ~ BOXCOX, value.var = "METRICS")

# Attach those lagging metrics to the table  
dat <- lagging_metrics[dat, on = .(ADP, EM_GEAR)]

# Consolidate rows where the lagging metric was the same  
dat <- merge(dat[, .(EM_GEAR, ADP, N, n, r, N_BOXCOX, n_u, r_u)], dat[, .(EM_GEAR, ADP, N, n, r, Y_BOXCOX, n_t, r_t)], by.x = c("EM_GEAR", "ADP", "N", "n", "r", "N_BOXCOX"), by.y = c("EM_GEAR", "ADP", "N", "n", "r", "Y_BOXCOX"), all = TRUE)

# Reorder and rename columns
dat <- dat[, .(Gear = EM_GEAR, Year = ADP, N, n, r, m = N_BOXCOX, n_u, n_t, r_u, r_t)]

# Add totals
totals <- dat[, lapply(.SD, function(x) sum(x, na.rm = TRUE)), .SDcols = c("N", "n", "n_u", "n_t")]
totals <- totals[, ':=' (Gear = "Total", Year = "", r = round(n / N * 100, 2), m = "", r_u = round(n_u / N * 100, 2), r_t = round(n_t / N * 100, 2))]
dat    <- rbind(dat, totals)

# Replace NAs with a dash
dat <- dat[, lapply(.SD, function(x) ifelse(is.na(x), "-", as.character(x)))]

# Create Table 1
t1 <- flextable(dat) %>% 
      merge_v(j = "Gear") %>% 
      compose(j = "n_u", part = "header", value = as_paragraph("n", as_chunk("u", props = fp_text(vertical.align = "subscript")))) %>% 
      compose(j = "r_u", part = "header", value = as_paragraph("r", as_chunk("u", props = fp_text(vertical.align = "subscript")))) %>%  
      compose(j = "n_t", part = "header", value = as_paragraph("n", as_chunk("t", props = fp_text(vertical.align = "subscript")))) %>% 
      compose(j = "r_t", part = "header", value = as_paragraph("r", as_chunk("t", props = fp_text(vertical.align = "subscript")))) %>%  
      border_remove() %>% 
      hline_top(border = fp_border(color = "black", style = "solid", width = 1)) %>% 
      hline(i = c(2, 5, 8, 11), border = fp_border(color = "black", style = "solid", width = 1)) %>% 
      hline(i = 12, j = 1, border = fp_border(color = "black", style = "solid", width = 1)) %>% 
      hline(i = 14, border = fp_border(color = "black", style = "solid", width = 1))

# Save Table 1
# save_as_docx(t1, path = "Observer_Effects/tables/t1.docx")

# Table 2: initial lambda estimates on raw data, before any uplabeling

# Put data in long format
dat <- melt.data.table(perm.data, measure.vars = metrics, variable.name = "METRICS", value.name = "vals")

# Center values about the mean within group and yes / no variable
dat[, vals_centered := vals - mean(vals), keyby = c(groups, "METRICS", yn_var)]

# Shift values to be positive within each group
dat[, vals_shifted := vals_centered - min(vals_centered) + 0.01, keyby = c(groups, "METRICS")]

# Calculate the likelihood of many values of Box-Cox lambda for each group and metric
bc <- dat[, .(lambda = seq(-10, 10, 0.01), logLik = MASS::boxcox(vals_shifted ~ get(yn_var), lambda = seq(-10, 10, 0.01), plotit = FALSE)$y), keyby = c(groups, "METRICS")]

# Find the value of lambda that maximizes the likelihood for each group and metric
bc <- bc[, .(lambda = lambda[logLik == max(logLik)]), keyby = c(groups, "METRICS")]

# Put estimates into wide format
bc <- dcast(bc, ... ~ METRICS, value.var = "lambda")

# Reorder rows
bc <- bc[order(EM_GEAR, ADP)]

# Rename columns
bc <- bc[, .(Gear = EM_GEAR, Year = as.character(ADP), `NMFS areas` = N_AREAS, `Days fished` = DAYS, `Vessel length (ft)` = LENGTH_OVERALL, `Species landed` = N_SPECIES, pMax = PMAX, `Landed catch (t)` = LANDED_CATCH)]

# Create Table 2
t2 <- flextable(bc) %>% 
      merge_v(j = "Gear") %>%
      border_remove() %>% 
      hline_top(border = fp_border(color = "black", style = "solid", width = 1)) %>% 
      hline(i = c(2, 5, 8, 11), border = fp_border(color = "black", style = "solid", width = 1)) %>% 
      hline(i = 12, j = 1, border = fp_border(color = "black", style = "solid", width = 1)) %>% 
      hline(i = 14, border = fp_border(color = "black", style = "solid", width = 1))

# Save Table 2
# save_as_docx(t2, path = "Observer_Effects/tables/t2.docx")

# Figures -----------------------------------------------------------------

## Figure 2: median p-values against observation rate for each metric, gear type, and transformation status
dat <- uplabel_results$N_table[uplabel_results$`p-value`, on = .(ADP, EM_GEAR, POP, TRIPS_ADDED)][, POP := NULL][, BOXCOX := "No"]
dat <- rbind(dat, uplabel_results_trans$N_table[uplabel_results_trans$`p-value`, on = .(ADP, EM_GEAR, POP, TRIPS_ADDED)][, POP := NULL][, BOXCOX := "Yes"])
dat <- rbind(dat, perm_results$N_table[perm_results$`p-value`, on = .(ADP, EM_GEAR)][, ":=" (TRIPS_ADDED = 0, BOXCOX = "No")])
dat <- rbind(dat, perm_results_trans$N_table[perm_results_trans$`p-value`, on = .(ADP, EM_GEAR)][, ":=" (TRIPS_ADDED = 0, BOXCOX = "Yes")])
dat <- dat[, lapply(.SD, function(x) median(x, na.rm = TRUE)), .SDcols = metrics, keyby = .(ADP, EM_GEAR, N, n, rate, TRIPS_ADDED, BOXCOX)]
dat <- unique(dat[, .(ADP, EM_GEAR, N, n, rate, TRIPS_ADDED, N_AREAS, DAYS, LENGTH_OVERALL, N_SPECIES, PMAX, LANDED_CATCH, BOXCOX)])
dat <- melt(dat, measure.vars = metrics, variable.name = "METRIC")
dat[, EM_GEAR := factor(EM_GEAR, levels = c("EM HAL", "HAL", "POT", "NPT", "PTR"))]
dat[, BOXCOX := factor(BOXCOX, levels = c("Yes", "No"))]
dat[, METRIC := recode(METRIC, "N_AREAS" = "NMFS\nareas", "DAYS" = "Days\nfished", "LENGTH_OVERALL" = "Vessel\nlength (ft)", "N_SPECIES" = "Species\nlanded", "PMAX" = "pMax\nspecies", "LANDED_CATCH" = "Landed\ncatch (t)")]
dat <- dat[ADP == 2019]

# Create Figure 2
f2 <- ggplot(dat, aes(x = rate * 100, y = value, shape = BOXCOX, color = BOXCOX)) +
      geom_line() +
      geom_point(size = 2) + 
      geom_hline(yintercept = 0.05, linetype = 2, color = "red") +
      facet_grid(METRIC ~ EM_GEAR) +
      labs(x = "Simulated monitoring rate (%)", y = "p-value", shape = "Transformed", color = "Transformed") +
      theme_bw() +
      scale_color_manual(values = c("red", "black")) +
      theme(text = element_text(size = 16),
            legend.position = "bottom")

# Save Figure 2
# png("Observer_Effects/figures/f2.png", width = 11, height = 8.5, units = 'in', res = 300)
# f2
# dev.off()

## Figure 3: median difference (%) against observation rate for each metric, gear type, and transformation status  
dat <- uplabel_results$N_table[uplabel_results$obs[SCALE == "PCT"], on = .(ADP, EM_GEAR, POP, TRIPS_ADDED)][, POP := NULL][, BOXCOX := "No"]
dat <- rbind(dat, uplabel_results_trans$N_table[uplabel_results_trans$obs[SCALE == "PCT"], on = .(ADP, EM_GEAR, POP, TRIPS_ADDED)][, POP := NULL][, BOXCOX := "Yes"])
dat <- rbind(dat, perm_results$N_table[perm_results$obs[SCALE == "PCT"], on = .(ADP, EM_GEAR)][, ":=" (TRIPS_ADDED = 0, BOXCOX = "No")])
dat <- rbind(dat, perm_results_trans$N_table[perm_results_trans$obs[SCALE == "PCT"], on = .(ADP, EM_GEAR)][, ":=" (TRIPS_ADDED = 0, BOXCOX = "Yes")])
dat <- dat[, lapply(.SD, function(x) median(x, na.rm = TRUE)), .SDcols = metrics, keyby = .(ADP, EM_GEAR, N, n, rate, TRIPS_ADDED, BOXCOX)]
dat <- unique(dat[, .(ADP, EM_GEAR, N, n, rate, TRIPS_ADDED, N_AREAS, DAYS, LENGTH_OVERALL, N_SPECIES, PMAX, LANDED_CATCH, BOXCOX)])
dat <- melt(dat, measure.vars = metrics, variable.name = "METRIC")
dat[, EM_GEAR := factor(EM_GEAR, levels = c("EM HAL", "HAL", "POT", "NPT", "PTR"))]
dat[, BOXCOX := factor(BOXCOX, levels = c("Yes", "No"))]
dat[, METRIC := recode(METRIC, "N_AREAS" = "NMFS\nareas", "DAYS" = "Days\nfished", "LENGTH_OVERALL" = "Vessel\nlength (ft)", "N_SPECIES" = "Species\nlanded", "PMAX" = "pMax\nspecies", "LANDED_CATCH" = "Landed\ncatch (t)")]
dat <- dat[ADP == 2019]

# Create Figure 3
f3 <- ggplot(dat, aes(x = rate * 100, y = value, shape = BOXCOX, color = BOXCOX)) +
      geom_line() +
      geom_point(size = 2) + 
      geom_hline(yintercept = 0.05, linetype = 2, color = "red") +
      facet_grid(METRIC ~ EM_GEAR) +
      labs(x = "Simulated monitoring rate (%)", y = "Observed difference (%)", shape = "Transformed", color = "Transformed") +
      theme_bw() +
      scale_color_manual(values = c("red", "black")) +
      theme(text = element_text(size = 16),
            legend.position = "bottom")

# Save Figure 3
# png("Observer_Effects/figures/f3.png", width = 11, height = 8.5, units = 'in', res = 300)
# f3
# dev.off()
