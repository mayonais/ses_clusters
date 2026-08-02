library(tidyverse)
library(lubridate)
library(INLA)
library(dplyr)
library(tidyr)

area <- "Los Angeles"
folder_name <- "Los_Angeles_discrete_binarycluster"

clustered_zctas <- readRDS(paste0(
  "create_cluster outputs/", folder_name, "/", area, "_ACS_zcta_clustered_for_ED.rds"))

daily_ed <- readRDS(paste0(
  "create_cluster outputs/", folder_name, "/", area, "_daily_ed_no_heat.rds"))

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

# standardize daily heat index across the full study period
heat_mean <- mean(daily_ed$Max_HI_Value, na.rm = TRUE)
heat_sd <- sd(daily_ed$Max_HI_Value, na.rm = TRUE)

daily_ed <- daily_ed %>%
  mutate(heat_z = (Max_HI_Value - heat_mean) / heat_sd) %>%
  group_by(zip) %>% mutate(
    chronic_heat = mean(heat_z, na.rm = TRUE),
    acute_heat = heat_z - chronic_heat) %>% ungroup()

# -----------------------------------------------------------------------------

daily_ed <- daily_ed %>%
  mutate(day_of_week = factor(wday(date, label = TRUE), ordered = FALSE),
    doy = as.integer(doy), cluster = fct_relevel(factor(as.character(cluster)),
      as.character(reference_cluster)))

daily_ed <- daily_ed %>%
  select(-`...1`, -D1Dx1, -D3Dx1, -D4Dx1, -Percentile.95, -month,
         -heat_day, -D2_suppressed, -n_days,
         -n_suppressed_days, -daily_sum_raw, -D2Dx1)

head(daily_ed)
gc()

# ------ POISSON (only with heat) ----------------------------------------------------------

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
    D2 ~ cluster + day_of_week + chronic_heat +
      f(acute_heat, model = "rw2", scale.model = TRUE) +
      day_of_week + f(doy, model = "rw2", scale.model = TRUE),
    family = "poisson", data = daily_ed, E = Population,
    control.compute = list(dic = TRUE, waic = TRUE),
    control.predictor = list(compute = TRUE))
  
  print(summary(poisson_inla))
  cat("\n\n====================================================================\n")
  cat("\nDIC:\n")
  print(poisson_inla$dic$dic)
  cat("\nWAIC:\n")
  print(poisson_inla$waic$waic)
  
  poisson_p_table <- as.data.frame(poisson_inla$summary.fixed)
  poisson_p_table$term <- rownames(poisson_p_table)
  
  saveRDS(poisson_p_table, paste0("create_cluster outputs/", folder_name, "/",
                             area, "_suppressed_INLA_poisson_coef_table.rds"))  
}, file = paste0("create_cluster outputs/", folder_name, "/",
                 area, "_suppressed_ED_INLA_poisson_summary.txt"))
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
    D2 ~ cluster + day_of_week +  f(doy, model = "rw2", scale.model = TRUE)#,
    + chronic_heat + f(acute_heat, model = "rw2", scale.model = TRUE),
    family = "nbinomial", data = daily_ed, E = Population,
  
    control.compute = list(dic = TRUE, waic = TRUE),
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
                             area, "_suppressed_INLA_nb_coef_table.rds"))
}, file = paste0("create_cluster outputs/", folder_name, "/",
                 area, "_suppressed_ED_INLA_nb_summary.txt"))
rm(nb_inla, expected_counts)
gc()

# ----- extract INLA coefficients and calculate IRRs -------------------------

nb_p_table <- readRDS(paste0("create_cluster outputs/", folder_name, "/",
                             area, "_suppressed_INLA_nb_coef_table.rds")) %>%
  mutate(p_IRR_gt_1 = 1 - pnorm(0, mean = mean, sd = sd))

coef_table <- nb_p_table %>%
  rename(estimate = mean, std.error = sd, ci_low_log = `0.025quant`,
         ci_high_log = `0.975quant`) %>%
  filter(term != "(Intercept)") %>% mutate(IRR = exp(estimate),
         ci_low = exp(ci_low_log), ci_high = exp(ci_high_log),
         predictor = case_when(
           grepl("^cluster", term) ~ "cluster",
           grepl("^day_of_week", term) ~ "day_of_week",
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

composite_test <- cor.test(final_table$vulnerability_score,
                           final_table$IRR, method = "spearman")

irr_corrs <- lapply(vars, function(v) {
  test <- suppressWarnings(cor.test(final_table[[v]], final_table$IRR,
             method = "spearman"))
  data.frame(variable = v, rho = unname(test$estimate))
}) %>% bind_rows() %>% arrange(desc(abs(rho)))

write.csv(final_table, paste0("create_cluster outputs/", folder_name, "/",
                 area, "_cluster_results.csv"), row.names = FALSE)

write.csv(irr_corrs, paste0("create_cluster outputs/", folder_name, "/",
                 area, "_vulnerability_IRR_correlation.csv"), row.names = FALSE)
