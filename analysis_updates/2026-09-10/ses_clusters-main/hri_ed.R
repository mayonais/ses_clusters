library(tidyverse)
library(lubridate)
library(INLA)
library(dplyr)
library(tidyr)
source("cluster_ar.R")

file_name <- "Los Angeles"
folder_name <- "Los_Angeles_discrete_none_HEAT_by_cluster"

clustered_zctas <- readRDS(paste0(
  folder_name, "/", file_name, "_ACS_zcta_clustered_for_ED.rds"))

daily_ed <- readRDS(paste0(
  folder_name, "/", file_name, "_daily_ed_no_heat.rds"))

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

# Joint posterior marginals for shared heat slope + cluster deviation.
# INLA linear combinations retain covariance between the two coefficients.
cluster_terms <- paste0("cluster", levels(daily_ed$cluster))
heat_lincombs <- setNames(lapply(levels(daily_ed$cluster), function(cluster) {
  weights <- list(acute_heat = 1)
  if (cluster != reference_cluster) {
    weights[[paste0("cluster", cluster, ":acute_heat")]] <- 1
  }
  do.call(inla.make.lincomb, weights)[[1]]
}), cluster_terms)

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
    D2 ~ cluster + chronic_heat + acute_heat + acute_heat:cluster + month + day_of_week +
      f(doy, model = "rw2", scale.model = TRUE, constr = TRUE),
    family = "poisson", data = daily_ed, E = Population,
    lincomb = heat_lincombs,
    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
    control.predictor = list(compute = TRUE))
  
  # f(acute_heat_group, model = "rw2", scale.model = TRUE, constr = TRUE)
  
  print(summary(poisson_inla))
  cat("\n\n====================================================================\n")
  cat("\nDIC:\n")
  print(poisson_inla$dic$dic)
  cat("\nWAIC:\n")
  print(poisson_inla$waic$waic)
  
  cat("\n====================================================================\n")
  
  cat("\nZERO-INFLATION DIAGNOSTIC (observed vs. Poisson-predicted)\n")
  observed_zero_prop <- mean(daily_ed$D2 == 0)
  expected_counts <- poisson_inla$summary.fitted.values$mean * daily_ed$Population
  predicted_zero_prop <- mean(dpois(0, lambda = expected_counts))
  cat("Observed zero proportion:", observed_zero_prop, "\n")
  cat("Poisson-predicted zero proportion:", predicted_zero_prop, "\n")
  cat("Excess zeros (observed - predicted):", observed_zero_prop - predicted_zero_prop, "\n")
  
  poisson_p_table <- as.data.frame(poisson_inla$summary.fixed)
  poisson_p_table$term <- rownames(poisson_p_table)
  
  saveRDS(poisson_p_table, paste0(folder_name, "/",
                                  folder_name, "_suppressed_INLA_poisson_coef_table.rds"))  
}, file = paste0(folder_name, "/",
                 folder_name, "_suppressed_ED_INLA_poisson_summary.txt"))
saveRDS(poisson_inla, paste0(folder_name, "/",
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
      D2 ~ cluster + chronic_heat + acute_heat + acute_heat:cluster + month + day_of_week +
        f(doy, model = "rw2", scale.model = TRUE, constr = TRUE),
      family = "nbinomial", data = daily_ed, E = Population,
      lincomb = heat_lincombs,
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
    
    saveRDS(nb_p_table, paste0(folder_name, "/",
                               file_name, "_suppressed_INLA_nb_coef_table.rds"))
  }, file = paste0(folder_name, "/",
                   file_name, "_suppressed_ED_INLA_nb_summary.txt"))
saveRDS(nb_inla, paste0(folder_name, "/",
                        folder_name, "_suppressed_INLA_nb_model.rds"))
rm(nb_inla, expected_counts)
gc()

# ===========================================================================
# extract INLA coefficients and calculate IRRs
# ===========================================================================

poisson_inla <- readRDS(paste0(folder_name, "/",
                              folder_name, "_suppressed_INLA_poisson_model.rds"))
p_table <- as.data.frame(poisson_inla$summary.fixed)
p_table$term <- rownames(p_table)
p_table$p_IRR_gt_1 <- vapply(poisson_inla$marginals.fixed[p_table$term],
  function(marginal) 1 - inla.pmarginal(0, marginal), numeric(1))

# Cluster main effects compare rates at acute_heat = 0, not heat slopes.
# Exact term matching excludes cluster-by-heat interactions from this table.
coef_table <- p_table %>%
  filter(term %in% cluster_terms) %>%
  transmute(predictor = "cluster", cluster = as.integer(sub("^cluster", "", term)),
            IRR = exp(mean), ci_low = exp(`0.025quant`),
            ci_high = exp(`0.975quant`), p_IRR_gt_1)

reference_row <- tibble(
  predictor = "cluster", cluster = as.integer(reference_cluster),
  IRR = 1, ci_low = 1, ci_high = 1, p_IRR_gt_1 = NA_real_)

final_table <- bind_rows(reference_row, coef_table) %>%
  left_join(cluster_labels %>% mutate(cluster = as.integer(cluster)),
            by = "cluster") %>%
  left_join(cluster_vuln, by = "cluster") %>%
  left_join(cluster_profiles %>%
              mutate(cluster = as.integer(as.character(cluster))),
            by = "cluster") %>% arrange(desc(IRR))

write.csv(final_table, paste0(folder_name, "/",
                 file_name, "_poisson_cluster_results.csv"), row.names = FALSE)

# ===========================================================================
# POSTERIOR CLUSTER-AVERAGE ATTRIBUTABLE RISK
# ===========================================================================

# ----- model coefficients: cluster-specific intercepts and slopes ---------

p_table <- as.data.frame(poisson_inla$summary.fixed)
p_table$term <- rownames(p_table)

reference_term <- paste0("cluster", reference_cluster)
alpha_ref <- p_table %>% filter(term == "(Intercept)") %>% pull(mean)

cluster_intercepts <- tibble(cluster = paste0(
  "cluster",sort(unique(daily_ed$cluster))), alpha = alpha_ref) %>%
  left_join(p_table %>% filter(term %in% cluster_terms) %>%
              transmute(cluster = term, cluster_effect = mean), by = "cluster") %>%
  mutate(alpha = ifelse(cluster == reference_term, alpha,
                        alpha + coalesce(cluster_effect, 0))) %>%
  select(cluster, alpha)

# ----------------------------------------------------------------------

slope_marginals <- poisson_inla$marginals.lincomb.derived[cluster_terms]
heat_slopes <- tibble(
  cluster = cluster_terms,
  beta = poisson_inla$summary.lincomb.derived[cluster_terms, "mean"],
  P_beta_gt_0 = vapply(slope_marginals,
    function(marginal) 1 - inla.pmarginal(0, marginal), numeric(1)))

# Equal-weight average of signed ZIP-day rate contrasts within each draw.
cluster_posterior <- cluster_ar_exceedance(
  poisson_inla, daily_ed, reference_cluster,
  n_samples = 4000L, batch_size = 200L, seed = 20260910L)
saveRDS(cluster_posterior, paste0(folder_name, "/",
  file_name, "_poisson_cluster_AR_posterior.rds"))

cluster_results <- cluster_intercepts %>%
  left_join(heat_slopes, by = "cluster") %>%
  left_join(cluster_posterior$summary, by = "cluster") %>%
  arrange(as.numeric(sub("cluster", "", cluster)))

write.csv(cluster_results, paste0(folder_name, "/",
  file_name, "_poisson_cluster_AR_exceedance.csv"), row.names = FALSE)

print(cluster_results)
gc()
