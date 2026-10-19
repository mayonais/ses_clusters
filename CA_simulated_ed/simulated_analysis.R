library(tidyverse)
library(lubridate)
library(INLA)
library(dplyr)
library(tidyr)

state <- "CA"


file_name <- "Los Angeles"
folder_name <- "Los_Angeles_discrete_none_HEAT_by_cluster"

clustered_zctas <- readRDS(paste0(
  folder_name, "/", file_name, "_ACS_zcta_clustered_for_ED.rds"))

ed_data <- read.csv("CA_simulated_ed/daily_counts_2017_2018.csv",
                    colClasses = c(GEOID = "character"))

ses_clusters <- clustered_zctas %>% dplyr::select(GEOID, cluster, cluster_label)
ed_combined <- ed_data %>% inner_join(ses_clusters, by = "GEOID")

daily_ed <- ed_combined %>% mutate(cluster = factor(cluster)) %>%
  filter(!is.na(cluster)) %>% mutate(date = as.Date(date),
                                     year = factor(year(date)),
                                     day_of_week = factor(wday(date, label = TRUE), ordered = FALSE), # matches CHS's delta_DOW(d)
                                     day_of_year = yday(date)) %>%
  dplyr::select(hri_ed, cluster, pop, year, day_of_week, day_of_year)

head(daily_ed)

rm(ed_data)
rm(ed_combined)
rm(ses_clusters)
gc()

# ============================================================
# set lowest vulnerability score as reference cluster
# ============================================================

cluster_labels <- clustered_zctas %>% distinct(cluster, cluster_label) %>%
  mutate(cluster = as.character(cluster))

vars <- c("pct_white", "pct_black", "pct_asian", "pct_hispanic",
          "pct_age17", "pct_age65", "poverty_rate", "median_income",
          "unemployment_rate", "education_low", "uninsured_rate", "vehicle_rate",
          "language_isolation", "disability_rate", "housing_overcrowding",
          "owner_burden_rate", "renter_burden_rate", "outdoor_worker_rate", "renter_rate",
          "no_phone_rate","elderly_alone_rate", "alt_housing_rate", "old_housing_rate")

cluster_profiles <- clustered_zctas %>% group_by(cluster) %>%
  summarise(across(all_of(vars), ~ median(.x, na.rm = TRUE)), n = n())

source("shared_functions.R")

cluster_vuln <- clustered_zctas %>%
  compute_vulnerability_score() %>%
  group_by(cluster) %>%
  summarise(vulnerability_score = mean(vulnerability_score, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  arrange(vulnerability_score)

reference_cluster <- as.character(cluster_vuln$cluster[1])

# set reference cluster for INLA
daily_ed <- daily_ed %>%
  mutate(cluster = forcats::fct_relevel(cluster, reference_cluster))
head(daily_ed)

cluster_terms <- paste0("cluster", setdiff(levels(daily_ed$cluster), reference_cluster))

# cluster_labels, cluster_vuln, cluster_profiles kept for the results table below
rm(clustered_zctas)
gc()

# ========================================================================
# FIT POISSON AND NEGATIVE BINOMIAL MODELS WITH INLA
# ========================================================================

# ------ POISSON ----------------------------------------------------------

capture.output({
  
  cat("====================================================================\n")
  cat("POISSON INLA - SIMULATED 2018 ED DATA\n")
  cat("INLA version:\n")
  print(packageVersion("INLA"))
  cat("\nR version:\n")
  print(R.version.string)
  cat("\nReference cluster (lowest vulnerability):", reference_cluster, "\n")
  cat("====================================================================\n\n")
  
  poisson_inla <- inla(
    hri_ed ~ cluster + year + day_of_week +
      f(day_of_year, model = "rw2", scale.model = TRUE, constr = TRUE),
    family = "poisson",
    data = daily_ed,
    E = pop,
    control.compute = list(
      dic = TRUE,
      waic = TRUE,
      config = TRUE),
    control.predictor = list(
      compute = TRUE))
  
  print(summary(poisson_inla))
  cat("\n\n====================================================================\n")
  
  cat("\nDIC:\n")
  print(poisson_inla$dic$dic)
  cat("\nWAIC:\n")
  print(poisson_inla$waic$waic)
  
  cat("\n====================================================================\n")
  
  cat("\nZERO-INFLATION DIAGNOSTIC (observed vs. Poisson-predicted)\n")
  observed_zero_prop <- mean(daily_ed$hri_ed == 0)
  expected_counts <- poisson_inla$summary.fitted.values$mean * daily_ed$pop
  predicted_zero_prop <- mean(dpois(0, lambda = expected_counts))
  cat("Observed zero proportion:", observed_zero_prop, "\n")
  cat("Poisson-predicted zero proportion:", predicted_zero_prop, "\n")
  cat("Excess zeros (observed - predicted):", observed_zero_prop - predicted_zero_prop, "\n")
  
  poisson_p_table <- as.data.frame(poisson_inla$summary.fixed)
  poisson_p_table$term <- rownames(poisson_p_table)
  
  saveRDS(poisson_p_table, paste0(folder_name, "/",
                                  folder_name, "_simulated_INLA_poisson_coef_table.rds"))
}, file = paste0(folder_name, "/",
                 folder_name, "_simulated_ED_INLA_poisson_summary.txt"))
saveRDS(poisson_inla, paste0(folder_name, "/",
                             folder_name, "_simulated_INLA_poisson_model.rds"))
rm(poisson_inla, expected_counts)
gc()

# ------- NEGATIVE BINOMIAL ---------------------------------------------------

capture.output({
  
  cat("====================================================================\n")
  cat("NEGATIVE BINOMIAL INLA - SIMULATED 2018 ED DATA\n")
  cat("INLA version:\n")
  print(packageVersion("INLA"))
  cat("\nR version:\n")
  print(R.version.string)
  cat("\nReference cluster (lowest vulnerability):", reference_cluster, "\n")
  cat("====================================================================\n\n")
  
  nb_inla <- inla(
    hri_ed ~ cluster + year + day_of_week +
      f(day_of_year, model = "rw2", scale.model = TRUE, constr = TRUE),
    family = "nbinomial",
    data = daily_ed,
    E = pop,
    control.compute = list(
      dic = TRUE,
      waic = TRUE,
      config = TRUE),
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
  observed_zero_prop <- mean(daily_ed$hri_ed == 0)
  expected_counts <- nb_inla$summary.fitted.values$mean * daily_ed$pop
  size_est <- nb_inla$summary.hyperpar[
    grep("size", rownames(nb_inla$summary.hyperpar), ignore.case = TRUE), "mean"][1]
  predicted_zero_prop <- mean(dnbinom(0, mu = expected_counts, size = size_est))
  cat("Observed zero proportion:", observed_zero_prop, "\n")
  cat("NB-predicted zero proportion:", predicted_zero_prop, "\n")
  cat("Excess zeros (observed - predicted):", observed_zero_prop - predicted_zero_prop, "\n")
  
  nb_p_table <- as.data.frame(nb_inla$summary.fixed)
  nb_p_table$term <- rownames(nb_p_table)
  
  saveRDS(nb_p_table, paste0(folder_name, "/",
                             folder_name, "_simulated_INLA_nb_coef_table.rds"))
  
}, file = paste0(folder_name, "/",
                 folder_name, "_simulated_ED_INLA_nb_summary.txt"))
saveRDS(nb_inla, paste0(folder_name, "/",
                        folder_name, "_simulated_INLA_nb_model.rds"))
rm(nb_inla, expected_counts)
gc()


# ============================================================
# extract INLA coefficients and calculate IRRs
# ============================================================

nb_inla <- readRDS(paste0(folder_name, "/",
                          folder_name, "_simulated_INLA_nb_model.rds"))

p_table <- as.data.frame(nb_inla$summary.fixed)
p_table$term <- rownames(p_table)
p_table$p_IRR_gt_1 <- vapply(nb_inla$marginals.fixed[p_table$term],
                             function(marginal) 1 - inla.pmarginal(0, marginal), numeric(1))

coef_table <- p_table %>%
  rename(
    estimate = mean,
    std.error = sd,
    ci_low_log = `0.025quant`,
    ci_high_log = `0.975quant`
  ) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    IRR = exp(estimate),
    ci_low = exp(ci_low_log),
    ci_high = exp(ci_high_log),
    pct_change = (IRR - 1) * 100,
    predictor = case_when(
      grepl("^cluster", term) ~ "cluster",
      grepl("^year", term) ~ "year",
      grepl("^day_of_week", term) ~ "day_of_week",
      TRUE ~ term
    ),
    level = sub("^(cluster|year|day_of_week)", "", term)
  ) %>%
  select(predictor, level, IRR, ci_low, ci_high, pct_change, p_IRR_gt_1)

season_effect <- nb_inla$summary.random$day_of_year %>%
  rename(day_of_year = ID, estimate = mean,
         ci_low_log = `0.025quant`, ci_high_log = `0.975quant`) %>%
  mutate(
    IRR = exp(estimate),
    ci_low = exp(ci_low_log),
    ci_high = exp(ci_high_log)
  ) %>%
  select(day_of_year, IRR, ci_low, ci_high)

# ============================================================
# attach cluster labels
# ============================================================

coef_table <- coef_table %>%
  mutate(
    cluster_label = if_else(
      predictor == "cluster",
      cluster_labels$cluster_label[
        match(level, cluster_labels$cluster)
      ],
      NA_character_
    )
  )

# ============================================================
# add reference cluster
# ============================================================

reference_row <- data.frame(
  predictor = "cluster",
  level = reference_cluster,
  IRR = 1,
  ci_low = 1,
  ci_high = 1,
  pct_change = 0,
  p_IRR_gt_1 = NA_real_,
  cluster_label = cluster_labels$cluster_label[
    match(reference_cluster, cluster_labels$cluster)
  ]
)

coef_table <- bind_rows(reference_row, coef_table) %>%
  mutate(
    predictor = factor(
      predictor,
      levels = c("cluster", "year", "day_of_week")
    ),
    level_numeric = suppressWarnings(as.numeric(level))
  ) %>%
  arrange(predictor, level_numeric) %>%
  select(-level_numeric)

rownames(coef_table) <- NULL

print(coef_table, row.names = FALSE)

write.csv(
  coef_table,
  paste0(folder_name, "/", file_name, "_simulated_NB_incidence_rate_ratios.csv"),
  row.names = FALSE
)


# ============================================================
# validate recovered cluster IRRs against injected truth
# ============================================================

cluster_truth <- read.csv(
  "CA_simulated_ed/cluster_multiplier.csv",
  colClasses = c(cluster = "character")
)

reference_multiplier <- cluster_truth %>%
  filter(cluster == as.character(reference_cluster)) %>%
  pull(true_multiplier)

validation <- coef_table %>%
  filter(predictor == "cluster") %>%
  left_join(
    cluster_truth %>%
      select(cluster, true_multiplier) %>%
      mutate(true_IRR = true_multiplier / reference_multiplier),
    by = c("level" = "cluster")
  ) %>%
  left_join(cluster_vuln %>% mutate(cluster = as.character(cluster)),
            by = c("level" = "cluster")) %>%
  left_join(cluster_profiles %>% mutate(cluster = as.character(cluster)),
            by = c("level" = "cluster")) %>%
  arrange(desc(IRR))


validation_table <- validation %>%
  select(
    cluster = level,
    cluster_label,
    estimated_IRR = IRR,
    assigned_IRR = true_IRR
  )


irr_correlation <- cor(
  validation$IRR,
  validation$true_IRR,
  method = "spearman"
)

print(validation_table)
print(paste("Spearman correlation:", round(irr_correlation, 3)))

write.csv(validation, paste0(folder_name, "/",
                             file_name, "_simulated_NB_cluster_results.csv"),
          row.names = FALSE)


# ============================================================
# plot estimated vs true IRR
# ============================================================

irr_plot <- ggplot(validation, aes(x = true_IRR, y = IRR)) +
  geom_point() +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  labs(
    x = "Assigned IRR from simulated ED data",
    y = "Estimated IRR from NB INLA",
    title = "Recovery of simulated cluster IRRs",
    subtitle = paste0(
      "Spearman correlation = ",
      round(irr_correlation, 3)
    )
  ) +
  theme_minimal()


pdf(
  "CA_simulated_ed/cluster_IRR_validation.pdf",
  width = 8,
  height = 8
)

print(irr_plot)

dev.off()
