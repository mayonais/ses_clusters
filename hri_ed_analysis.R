# Full HRI ED analysis: Poisson/NB diagnostics, cluster comparisons and excess rates.
# Shared calculations and compatible cached fits come from analysis_functions.R.
library(dplyr)
library(tidyr)
library(INLA)

# Locate inputs relative to this document, without changing the working directory.
script_files <- as.character(unlist(lapply(sys.frames(), function(frame) frame$ofile)))
script_files <- script_files[basename(script_files) == "hri_ed_analysis.R"]
command_file <- grep("^--file=", commandArgs(), value = TRUE)
project_dir <- if (length(script_files)) {
  dirname(normalizePath(tail(script_files, 1), mustWork = TRUE))
} else if (length(command_file) && sub("^--file=", "", command_file[1]) != "-") {
  dirname(normalizePath(sub("^--file=", "", command_file[1]), mustWork = TRUE))
} else {
  getwd() # Running selected lines: use the project's working directory.
}

# ---- Settings ---------------------------------------------------------------
file_name <- "Los Angeles"
folder_name <- "Los_Angeles_discrete_none_HEAT_by_cluster"
input_dir <- file.path(project_dir, folder_name)
model_cache_dir <- file.path(input_dir, "model_cache")
force_refit <- FALSE
ar_n_samples <- 4000L
ar_batch_size <- 200L
ar_seed <- 20260910L
source(file.path(project_dir, "analysis_functions.R"), local = TRUE)
source(file.path(project_dir, "cluster_ar.R"), local = TRUE)

prepared <- load_heat_data(input_dir, file_name)
daily_ed <- prepared$daily_ed
clustered_zctas <- prepared$clustered_zctas
cluster_labels <- prepared$cluster_labels
cluster_vuln <- prepared$cluster_vuln %>% mutate(cluster = as.integer(cluster))
reference_cluster <- prepared$reference_cluster
cluster_terms <- paste0("cluster", levels(daily_ed$cluster))

vars <- c("pct_white", "pct_black", "pct_asian", "pct_hispanic",
          "poverty_rate", "renter_burden_rate", "housing_overcrowding",
          "median_income", "education_low", "pct_age17", "no_phone_rate",
          "pct_age65", "uninsured_rate", "language_isolation", "disability_rate",
          "outdoor_worker_rate", "unemployment_rate", "elderly_alone_rate",
          "alt_housing_rate", "old_housing_rate")
cluster_profiles <- clustered_zctas %>% group_by(cluster) %>%
  summarise(across(all_of(vars), ~ median(.x, na.rm = TRUE)), n = n())

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
  
  poisson_inla <- fit_heat_model(prepared, model_cache_dir, family = "poisson", force_refit = force_refit)
  
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
  
  saveRDS(poisson_p_table, paste0(input_dir, "/",
                                  folder_name, "_suppressed_INLA_poisson_coef_table.rds"))  
}, file = paste0(input_dir, "/",
                 folder_name, "_suppressed_ED_INLA_poisson_summary.txt"))
saveRDS(poisson_inla, paste0(input_dir, "/",
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
  
  nb_inla <- fit_heat_model(prepared, model_cache_dir, family = "nbinomial", force_refit = force_refit)
  
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
  
  saveRDS(nb_p_table, paste0(input_dir, "/",
                             file_name, "_suppressed_INLA_nb_coef_table.rds"))
}, file = paste0(input_dir, "/",
                 file_name, "_suppressed_ED_INLA_nb_summary.txt"))
saveRDS(nb_inla, paste0(input_dir, "/",
                        folder_name, "_suppressed_INLA_nb_model.rds"))
rm(nb_inla, expected_counts)
gc()

# ===========================================================================
# extract INLA coefficients and calculate IRRs
# ===========================================================================

poisson_inla <- readRDS(paste0(input_dir, "/",
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

write.csv(final_table, paste0(input_dir, "/",
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

heat_slopes <- summarise_heat_slopes(poisson_inla, levels(daily_ed$cluster)) %>%
  transmute(cluster = paste0("cluster", cluster), beta = beta_inla_mean, P_beta_gt_0)
# Equal-weight average of signed ZIP-day rate contrasts within each draw.
cluster_posterior <- cluster_ar_exceedance(
  poisson_inla, daily_ed, reference_cluster,
  n_samples = ar_n_samples, batch_size = ar_batch_size, seed = ar_seed)
saveRDS(cluster_posterior, paste0(input_dir, "/",
                                  file_name, "_poisson_cluster_AR_posterior.rds"))

cluster_results <- cluster_intercepts %>%
  left_join(heat_slopes, by = "cluster") %>%
  left_join(cluster_posterior$summary, by = "cluster") %>%
  arrange(as.numeric(sub("cluster", "", cluster)))

write.csv(cluster_results, paste0(input_dir, "/",
                                  file_name, "_poisson_cluster_AR_exceedance.csv"), row.names = FALSE)

print(cluster_results)
gc()
