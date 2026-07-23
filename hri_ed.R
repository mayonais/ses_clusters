library(tidyverse)
library(lubridate)
library(INLA)
library(dplyr)
library(tidyr)
library(car)

state <- "CA"

clustered_zctas <- readRDS(paste0(
  "ed_analysis_ready/", gsub(" ", "_", state), "_ACS_zcta_clustered_for_ED_MIN.rds"))

#daily_ed <- readRDS("rds/CA_daily_ed_full.rds")
daily_ed <- readRDS("rds/CA_daily_ed_no_heat.rds")

cluster_labels <- clustered_zctas %>% distinct(cluster, cluster_label) %>%
  mutate(cluster = as.character(cluster))

#vars <- c("pct_black", "pct_asian", "pct_hispanic",
#          "pct_age17", "pct_age65", "poverty_rate", "median_income",
#          "unemployment_rate", "education_low", "uninsured_rate", "vehicle_rate",
#          "language_isolation", "disability_rate", "housing_overcrowding",
#          "owner_burden_rate", "renter_burden_rate", "outdoor_worker_rate", "renter_rate",
#          "no_phone_rate","elderly_alone_rate", "alt_housing_rate", "old_housing_rate")

clustered_zctas$dummy_outcome <- rnorm(nrow(clustered_zctas))
model <- lm(dummy_outcome ~ ., data = clustered_zctas[, c("dummy_outcome", vars)])
vif_results <- data.frame(variable = names(vif(model)), VIF = vif(model)) %>% arrange(desc(VIF))
print(vif_results)

vars <- c(
  "poverty_rate", "pct_white", "pct_black", "pct_asian", "pct_hispanic",
  "median_income", "education_low", "pct_age17",
  "pct_age65", "uninsured_rate", "language_isolation", "disability_rate",
  "outdoor_worker_rate", "renter_burden_rate", "housing_overcrowding",
  "alt_housing_rate")

cluster_profiles <- clustered_zctas %>% group_by(cluster) %>%
  summarise(across(all_of(vars), ~ median(.x, na.rm = TRUE)), n = n())

daily_ed <- daily_ed %>%
  mutate(cluster = clustered_zctas$cluster[match(zip, clustered_zctas$GEOID)])

sum(is.na(daily_ed$cluster))

# ============================================================
# set lowest vulnerability score as reference cluster
# ============================================================

# zcta-level vulnerability score: SVI-informed, equal-weighted domains
cluster_vuln <- clustered_zctas %>%  mutate(
    socioeconomic = rowMeans(cbind(
      percent_rank(poverty_rate),
      1 - percent_rank(median_income),
      #percent_rank(unemployment_rate),
      percent_rank(education_low),
      percent_rank(uninsured_rate)), na.rm = TRUE),
    
    household = rowMeans(cbind(
      percent_rank(pct_age65),
      percent_rank(pct_age17),
      percent_rank(disability_rate),
      percent_rank(language_isolation)),
      #percent_rank(vehicle_rate),
      #percent_rank(no_phone_rate),
      #percent_rank(elderly_alone_rate)),
      na.rm = TRUE),
    
    minority_status = percent_rank(1 - pct_white),
    
    housing = rowMeans(cbind(
      percent_rank(housing_overcrowding),
      #percent_rank(owner_burden_rate),
      percent_rank(renter_burden_rate),
      percent_rank(alt_housing_rate)),
      #percent_rank(old_housing_rate)),
      na.rm = TRUE),
    
    occupational_exposure = percent_rank(outdoor_worker_rate),
    
    vulnerability_score = rowMeans(cbind(
      socioeconomic, household, minority_status,
      housing, occupational_exposure
    ), na.rm = TRUE)) %>%
  select(GEOID, cluster, vulnerability_score) %>%
  arrange(vulnerability_score)

cluster_vuln %>% group_by(cluster) %>% summarise(n = n(),
    mean_vuln = mean(vulnerability_score, na.rm = TRUE),
    median_vuln = median(vulnerability_score, na.rm = TRUE),
    sd_vuln = sd(vulnerability_score, na.rm = TRUE),
    min_vuln = min(vulnerability_score, na.rm = TRUE),
    max_vuln = max(vulnerability_score, na.rm = TRUE)) %>% print(n=Inf)

cluster_vuln <- cluster_vuln %>% group_by(cluster) %>%
  summarise(vulnerability_score = mean(vulnerability_score, na.rm = TRUE),
            n = n()) %>% arrange(vulnerability_score)

reference_cluster <- as.character(cluster_vuln$cluster[1])

# ========================================================================
# FIT POISSON AND NEGATIVE BINOMIAL MODELS WITH INLA
# ========================================================================

daily_ed <- daily_ed %>%
  mutate(day_of_week = factor(wday(date, label = TRUE), ordered = FALSE),
    doy = as.integer(doy),
    cluster = fct_relevel(factor(as.character(cluster)),
      as.character(reference_cluster)))

daily_ed <- daily_ed %>%
  select(-`...1`, -D1Dx1, -D3Dx1, -D4Dx1, -Percentile.95, -month,
         -heat_day, -D2_suppressed, -n_days, -Max_HI_Value,
         -n_suppressed_days, -daily_sum_raw, -D2Dx1)
#daily_ed <- daily_ed %>% select(-HI_scaled)

daily_ed <- daily_ed %>%
  left_join(clustered_zctas %>%
      distinct(GEOID, poverty_rate, pct_hispanic, pct_black,
               outdoor_worker_rate, pct_age65, uninsured_rate),
    by = c("zip" = "GEOID"))

head(daily_ed)
gc()

# --------------------------------------------------------------

capture.output({
  
  cat("====================================================================\n")
  cat("POISSON INLA - SUPPRESSED 2018 ED DATA (uniform no-heat imputation + minimal set of SES vars)\n")
  cat("INLA version:\n")
  print(packageVersion("INLA"))
  cat("\nR version:\n")
  print(R.version.string)
  cat("\nReference cluster (lowest vulnerability):", reference_cluster, "\n")
  cat("====================================================================\n\n")

  poisson_inla <- inla(
    D2 ~ cluster + day_of_week +
      f(doy, model = "rw2", scale.model = TRUE),
    family = "poisson",
    data = daily_ed,
    E = Population,
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
  
  saveRDS(poisson_p_table, "rds/CA_suppressed_inla_poisson_coef_table_NO_HEAT_MIN.rds")
}, file = "outputs/CA_suppressed_ED_INLA_poisson_summary_NO_HEAT_MIN.txt")

rm(poisson_inla)
gc()

capture.output({
  
  cat("====================================================================\n")
  cat("NEGATIVE BINOMIAL INLA - SUPPRESSED 2018 ED DATA (uniform no-heat imputation + minimal set of SES vars)\n")
  cat("INLA version:\n")
  print(packageVersion("INLA"))
  cat("\nR version:\n")
  print(R.version.string)
  cat("\nReference cluster (lowest vulnerability):", reference_cluster, "\n")
  cat("====================================================================\n\n")
  
  nb_inla <- inla(
    D2 ~ cluster + day_of_week +
      f(doy, model = "rw2", scale.model = TRUE),
    family = "nbinomial",
    data = daily_ed,
    E = Population,
    control.compute = list(
      dic = TRUE,
      waic = TRUE),
    control.predictor = list(
      compute = TRUE))
  
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
  
  saveRDS(nb_p_table, "rds/CA_suppressed_inla_nb_coef_table_NO_HEAT_MIN.rds")
}, file = "outputs/CA_suppressed_ED_INLA_nb_summary_NO_HEAT_MIN.txt")
rm(nb_inla, expected_counts)
gc()

# ============================================================
# extract INLA coefficients and calculate IRRs
# ============================================================

nb_p_table <- readRDS("rds/CA_suppressed_inla_nb_coef_table_INDV.rds")

nb_p_table <- nb_p_table %>%
  mutate(p_IRR_gt_1 = 1 - pnorm(0, mean = mean, sd = sd))

coef_table <- nb_p_table %>%
  rename(estimate = mean, std.error = sd, ci_low_log = `0.025quant`, 
         ci_high_log = `0.975quant`) %>% filter(term != "(Intercept)") %>%
  mutate(IRR = exp(estimate), ci_low = exp(ci_low_log), 
         ci_high = exp(ci_high_log), 
         predictor = case_when(grepl("^cluster", term) ~ "cluster",
                               grepl("^day_of_week", term) ~ "day_of_week",
                               grepl("^doy", term) ~ "season", TRUE ~ term),
         level = gsub("[^0-9]+", "", term)) %>%
  select(predictor, level, IRR, ci_low, ci_high)

coef_table <- coef_table %>%  mutate(cluster_label = if_else(
    predictor == "cluster", cluster_labels$cluster_label[
      match(level, cluster_labels$cluster)], NA_character_))

reference_row <- data.frame(predictor = "cluster",
  level = as.character(reference_cluster),
  IRR = 1, ci_low = 1, ci_high = 1,
  cluster_label = cluster_labels$cluster_label[
    match(reference_cluster, cluster_labels$cluster)])

coef_table <- bind_rows(reference_row, coef_table) %>%
  mutate(predictor = factor(predictor, levels = c("cluster", "day_of_week", "season")),
         level_numeric = as.numeric(level)) %>%
  arrange(predictor, level_numeric) %>% select(-level_numeric)

rownames(coef_table) <- NULL

# ============================================================
# cluster IRR/CI + vulnerability score + SES variables
# ============================================================

svi <- read.csv("datasets/svi2022.csv") %>% mutate(FIPS = as.character(FIPS))

cluster_svi <- clustered_zctas %>% left_join(svi %>% 
              select(FIPS, RPL_THEMES), by = c("GEOID" = "FIPS")) %>%
  group_by(cluster) %>% summarise(svi_score = median(RPL_THEMES, na.rm = TRUE))

cluster_race <- clustered_zctas %>% group_by(cluster) %>%
  summarise(pct_white = median(pct_white, na.rm = TRUE),
            pct_black = median(pct_black, na.rm = TRUE),
            pct_asian = median(pct_asian, na.rm = TRUE),
            pct_hispanic = median(pct_hispanic, na.rm = TRUE))

final_table <- coef_table %>% filter(predictor == "cluster") %>%
  mutate(cluster = as.integer(level)) %>% left_join(
    nb_p_table %>% filter(grepl("^cluster", term)) %>%
      mutate(cluster = as.integer(gsub("[^0-9]", "", term))) %>%
      select(cluster, p_IRR_gt_1), by = "cluster") %>%
  mutate(p_IRR_gt_1 = if_else(cluster == reference_cluster, NA_real_, p_IRR_gt_1)) %>%
  left_join(cluster_vuln, by = "cluster") %>% 
  left_join(cluster_svi, by = "cluster") %>%
  left_join(cluster_profiles %>% 
              mutate(cluster = as.integer(as.character(cluster))), 
            by = "cluster") %>% 
  #left_join(cluster_race, by = "cluster") %>%
  arrange(desc(IRR)) %>% select(-level)

race_vars <- c("pct_white", "pct_black", "pct_asian", "pct_hispanic")
race_irr_corrs <- lapply(race_vars, function(v) {
  test <- cor.test(final_table[[v]], final_table$IRR, method = "spearman")
  data.frame(variable = v, rho = unname(test$estimate),
             p_value = test$p.value)}) %>% bind_rows() %>%
  mutate(p_bonferroni = p.adjust(p_value, "bonferroni"))

# correlations for each SES variable 
irr_corrs <- lapply(vars, function(v) {
  test <- suppressWarnings(cor.test(
    final_table[[v]], final_table$IRR, method = "spearman"))
  data.frame(variable = v, rho = unname(test$estimate), p_value = test$p.value)
}) %>% bind_rows() %>%
  mutate(p_bonferroni = p.adjust(p_value, method = "bonferroni")) %>%
  arrange(desc(abs(rho)))

composite_test <- cor.test(
  final_table$vulnerability_score, final_table$IRR, method = "spearman")


irr_corrs <- bind_rows(
  irr_corrs, data.frame(variable = "vulnerability_score_composite",
                      rho = unname(composite_test$estimate),
                      p_value = composite_test$p.value, p_bonferroni = NA_real_))

write.csv(final_table,
          "outputs/CA_suppressed_NB_cluster_results_NO_HEAT_MIN.csv", row.names = FALSE)
write.csv(irr_corrs, 
          "outputs/CA_suppressed_NB_cluster_IRR_SES_correlations_NO_HEAT.csv",
          row.names = FALSE)
#write.csv(race_irr_corrs, 
#          "outputs/CA_suppressed_NB_cluster_IRR_RACE_correlations_NO_HEAT_MIN_no_race.csv",
#          row.names = FALSE)