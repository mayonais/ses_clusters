library(tidyverse)
library(lubridate)
library(INLA)
library(dplyr)
library(tidyr)

file_name <- "Los Angeles"
folder_name <- "Los_Angeles_discrete_none_HEAT_by_cluster"

clustered_zctas <- readRDS(paste0(
  "create_cluster outputs/", folder_name, "/", file_name, "_ACS_zcta_clustered_for_ED.rds"))

daily_ed <- readRDS(paste0(
  "create_cluster outputs/", folder_name, "/", file_name, "_daily_ed_no_heat.rds"))

cluster_labels <- clustered_zctas %>% distinct(cluster, cluster_label) %>%
  mutate(cluster = as.character(cluster))

vars <- c("pct_white", "pct_black", "pct_asian", "pct_hispanic",
          "poverty_rate", "renter_burden_rate", "housing_overcrowding",
          "median_income", "education_low", "pct_age17", "no_phone_rate",
          "pct_age65", "uninsured_rate", "language_isolation", "disability_rate",
          "outdoor_worker_rate", "unemployment_rate", "elderly_alone_rate",
          "alt_housing_rate", "old_housing_rate")

cluster_profiles <- clustered_zctas %>% group_by(cluster) %>%
  summarise(across(all_of(vars), ~ median(.x, na.rm = TRUE)), n = n())

daily_ed <- daily_ed %>%
  mutate(cluster = clustered_zctas$cluster[match(zip, clustered_zctas$GEOID)])

sum(is.na(daily_ed$cluster))

# set lowest vulnerability score as reference cluster
source("vulnerability_score.R")

cluster_vuln <- clustered_zctas %>%
  compute_vulnerability_score() %>%
  group_by(cluster) %>%
  summarise(vulnerability_score = mean(vulnerability_score, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  arrange(vulnerability_score)

reference_cluster <- as.character(cluster_vuln$cluster[1])

# ----- CHRONIC AND ACUTE HEAT VARIABLES (CHS-style) -------------------------

heat_mean <- mean(daily_ed$Max_HI_Value, na.rm = TRUE)
heat_sd <- sd(daily_ed$Max_HI_Value, na.rm = TRUE)

daily_ed <- daily_ed %>%
  mutate(heat_z = (Max_HI_Value - heat_mean) / heat_sd) %>% # daily standardized heat
  group_by(zip) %>% mutate(
    chronic_heat = mean(heat_z, na.rm = TRUE), # zip's chronic mean heat 
    # how unusually hot/cold that day was for that ZIP compared to its long-term
    acute_heat = heat_z - chronic_heat) %>% ungroup() 

#acute_heat_group <- inla.group(daily_ed$acute_heat, n = 50)
#acute_heat_group <- match(acute_heat_group, unique(acute_heat_group))
#daily_ed$acute_heat_group <- acute_heat_group

# -----------------------------------------------------------------------------

daily_ed <- daily_ed %>%
  mutate(day_of_week = factor(wday(date, label = TRUE), ordered = FALSE),
    doy = as.integer(doy), month = factor(month(date)),
    cluster = fct_relevel(
      factor(as.character(cluster)), as.character(reference_cluster)))

daily_ed <- daily_ed %>%
  select(-`...1`, -D1Dx1, -D3Dx1, -D4Dx1, -Percentile.95,
         -heat_day, -D2_suppressed, -n_days,
         -n_suppressed_days, -daily_sum_raw, -D2Dx1)

head(daily_ed)
gc()

# ------ POISSON ----------------------------------------------------------

capture.output({
  cat("====================================================================\n")
  cat("POISSON INLA - SUPPRESSED 2018 ED DATA\n")
  cat("INLA version:\n")
  print(packageVersion("INLA"))
  cat("\nR version:\n")
  print(R.version.string)
  cat("\nReference cluster (lowest vulnerability):", reference_cluster, "\n")
  cat("====================================================================\n\n")
  
  poisson_inla <- inla(
    D2 ~ cluster + chronic_heat + acute_heat:cluster + month + day_of_week +
      f(doy, model = "rw2", scale.model = TRUE, constr = TRUE),
    family = "poisson", data = daily_ed, E = Population,
    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
    control.predictor = list(compute = TRUE))

    # f(acute_heat_group, model = "rw2", scale.model = TRUE, constr = TRUE)
  
  print(summary(poisson_inla))
  cat("\n\n====================================================================\n")
  cat("\nDIC:\n")
  print(poisson_inla$dic$dic)
  cat("\nWAIC:\n")
  print(poisson_inla$waic$waic)
  
  poisson_p_table <- as.data.frame(poisson_inla$summary.fixed)
  poisson_p_table$term <- rownames(poisson_p_table)
  
  saveRDS(poisson_p_table, paste0("create_cluster outputs/", folder_name, "/",
                                  folder_name, "_suppressed_INLA_poisson_coef_table.rds"))  
}, file = paste0("create_cluster outputs/", folder_name, "/",
                 folder_name, "_suppressed_ED_INLA_poisson_summary.txt"))
saveRDS(poisson_inla, paste0("create_cluster outputs/", folder_name, "/",
                             folder_name, "_suppressed_INLA_poisson_model.rds"))
rm(poisson_inla)
gc()

# ------- NEGATIVE BINOMIAL ---------------------------------------------------

capture.output({
  cat("====================================================================\n")
  cat("NEGATIVE BINOMIAL INLA - SUPPRESSED 2018 ED DATA\n")
  cat("INLA version:\n")
  print(packageVersion("INLA"))
  cat("\nR version:\n")
  print(R.version.string)
  cat("\nReference cluster (lowest vulnerability):", reference_cluster, "\n")
  cat("====================================================================\n\n")
  
  nb_inla <- inla(
    D2 ~ cluster + chronic_heat + acute_heat:cluster + month + day_of_week +
      f(doy, model = "rw2", scale.model = TRUE, constr = TRUE),
    family = "nbinomial", data = daily_ed, E = Population,
    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
    control.predictor = list(compute = TRUE))
  
  print(summary(nb_inla))
  cat("\n====================================================================\n")
  cat("\nDIC:\n")
  print(nb_inla$dic$dic)
  cat("\nWAIC:\n")
  print(nb_inla$waic$waic)
  cat("\n====================================================================\n")
  
  cat("\nZERO-INFLATION DIAGNOSTIC (observed vs. NB-predicted)\n")
  observed_zero_prop <- mean(daily_ed$D2 == 0)
  expected_counts <- nb_inla$summary.fitted.values$mean * daily_ed$Population
  size_est <- nb_inla$summary.hyperpar[
    grep("size", rownames(nb_inla$summary.hyperpar), ignore.case = TRUE), "mean"][1]
  predicted_zero_prop <- mean(dnbinom(0, mu = expected_counts, size = size_est))
  cat("Observed zero proportion:", observed_zero_prop, "\n")
  cat("NB-predicted zero proportion:", predicted_zero_prop, "\n")
  cat("Excess zeros (observed - predicted):", observed_zero_prop - predicted_zero_prop, "\n")
  
  nb_p_table <- as.data.frame(nb_inla$summary.fixed)
  nb_p_table$term <- rownames(nb_p_table)
  
  saveRDS(nb_p_table, paste0("create_cluster outputs/", folder_name, "/",
                             file_name, "_suppressed_INLA_nb_coef_table.rds"))
}, file = paste0("create_cluster outputs/", folder_name, "/",
                 file_name, "_suppressed_ED_INLA_nb_summary.txt"))
saveRDS(nb_inla, paste0("create_cluster outputs/", folder_name, "/",
                             folder_name, "_suppressed_INLA_nb_model.rds"))
rm(nb_inla, expected_counts)
gc()

# ===========================================================================
# extract INLA coefficients and calculate IRRs
# ===========================================================================

nb_p_table <- readRDS(paste0("create_cluster outputs/", folder_name, "/",
                             file_name, "_suppressed_INLA_nb_coef_table.rds")) %>%
  mutate(p_IRR_gt_1 = 1 - pnorm(0, mean = mean, sd = sd))

coef_table <- nb_p_table %>%
  rename(estimate = mean, std.error = sd, ci_low_log = `0.025quant`,
         ci_high_log = `0.975quant`) %>%
  filter(term != "(Intercept)") %>% mutate(IRR = exp(estimate),
         ci_low = exp(ci_low_log), ci_high = exp(ci_high_log),
         predictor = case_when(
           grepl("^cluster", term) ~ "cluster",
           grepl("^day_of_week", term) ~ "day_of_week",
           grepl("^month", term) ~ "month",
           grepl("^doy", term) ~ "season", TRUE ~ term),
         level = gsub("[^0-9]+", "", term)) %>%
  select(predictor, level, IRR, ci_low, ci_high) %>%
  mutate(cluster_label = if_else(predictor == "cluster",
    cluster_labels$cluster_label[match(level, cluster_labels$cluster)],
    NA_character_))

reference_row <- data.frame(
  predictor = "cluster", level = as.character(reference_cluster),
  IRR = 1, ci_low = 1, ci_high = 1,
  cluster_label = cluster_labels$cluster_label[
    match(reference_cluster, cluster_labels$cluster)])

coef_table <- bind_rows(reference_row, coef_table) %>%
  mutate(predictor = factor(predictor, levels = c("cluster", "day_of_week", "season")),
         level_numeric = as.numeric(level)) %>%
  arrange(predictor, level_numeric) %>% select(-level_numeric)

final_table <- coef_table %>% filter(predictor == "cluster") %>%
  mutate(cluster = as.integer(level)) %>% left_join(nb_p_table %>%
      filter(grepl("^cluster", term)) %>%
      mutate(cluster = as.integer(gsub("[^0-9]", "", term))) %>%
      select(cluster, p_IRR_gt_1), by = "cluster") %>%
  mutate(p_IRR_gt_1 = if_else(cluster == as.integer(reference_cluster),
                              NA_real_, p_IRR_gt_1)) %>%
  left_join(cluster_vuln, by = "cluster") %>%
  left_join(cluster_profiles %>%
              mutate(cluster = as.integer(as.character(cluster))),
            by = "cluster") %>% arrange(desc(IRR)) %>% select(-level)

write.csv(final_table, paste0("create_cluster outputs/", folder_name, "/",
                 file_name, "_cluster_results.csv"), row.names = FALSE)

# ===========================================================================
# CHS-STYLE CLUSTER THRESHOLDS USING CLUSTER-SPECIFIC LINEAR HEAT EFFECTS
# ===========================================================================

inla <- readRDS(paste0("create_cluster outputs/", folder_name, "/",
                       folder_name, "_suppressed_INLA_poisson_model.rds"))

# ----- model coefficients: cluster-specific intercepts and slopes ---------

p_table <- as.data.frame(inla$summary.fixed)
p_table$term <- rownames(p_table)

reference_cluster <- "cluster6"
alpha_ref <- p_table %>% filter(term == "(Intercept)") %>% pull(mean)

cluster_intercepts <- tibble(cluster = paste0(
  "cluster",sort(unique(daily_ed$cluster))), alpha = alpha_ref) %>%
  left_join(p_table %>% filter(grepl("^cluster", term)) %>%
              transmute(cluster = term, cluster_effect = mean), by = "cluster") %>%
  mutate(alpha = ifelse(cluster == reference_cluster, alpha,
                        alpha + coalesce(cluster_effect, 0))) %>%
  select(cluster, alpha)

heat_slopes <- p_table %>%
  filter(grepl("^cluster.*:acute_heat$", term)) %>%
  transmute(cluster = sub(":acute_heat$", "", term), beta = mean)

cluster_model_results <- cluster_intercepts %>%
  left_join(heat_slopes, by = "cluster") %>%
  arrange(as.numeric(sub("cluster", "", cluster)))

# ----- cluster-specific chronic heat --------------------------------------

cluster_heat <- daily_ed %>%
  mutate(cluster_name = paste0("cluster", as.character(cluster))) %>%
  group_by(cluster_name) %>%
  summarise(H_bar = mean(chronic_heat, na.rm = TRUE), .groups = "drop")

threshold_data <- heat_slopes %>%
  left_join(cluster_heat, by = c("cluster" = "cluster_name"))

# ----- daily counterfactual rate and attributable risk --------------------

daily_thresholds <- daily_ed %>%
  mutate(cluster_name = paste0("cluster", as.character(cluster)),
         beta = threshold_data$beta[match(cluster_name, threshold_data$cluster)],
         H_bar = threshold_data$H_bar[match(cluster_name, threshold_data$cluster)],
         lp_obs = inla$summary.linear.predictor$mean,
         acute_heat = heat_z - H_bar,
         lp_cf = lp_obs - beta * acute_heat,
         rate_obs = exp(lp_obs),
         rate_cf = exp(lp_cf),
         AR = pmax((rate_obs - rate_cf) * 10000, 0)) # excess rate per 10,000 population

# ----- cluster-specific intercepts, slopes, and P(AR > 0) -------------------

cluster_slopes <- p_table %>%
  filter(grepl("^cluster.*:acute_heat$", term)) %>%
  transmute(cluster = sub(":acute_heat$", "", term), beta = mean)

cluster_exceedance <- daily_thresholds %>%
  group_by(cluster_name) %>%
  summarise(P_AR_gt_0 = mean(AR > 0, na.rm = TRUE), .groups = "drop") %>%
  rename(cluster = cluster_name)

cluster_results <- cluster_intercepts %>%
  left_join(cluster_slopes, by = "cluster") %>%
  left_join(cluster_exceedance, by = "cluster") %>%
  arrange(as.numeric(sub("cluster", "", cluster)))

print(cluster_results)
gc()
