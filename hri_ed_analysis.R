# Full HRI ED analysis: Poisson/NB diagnostics, cluster comparisons and excess rates.
# Shared calculations and compatible cached fits come from shared_functions.R.
library(dplyr)
library(tidyr)
library(INLA)

script_files <- as.character(unlist(lapply(sys.frames(), function(frame) frame$ofile)))
script_files <- script_files[basename(script_files) == "hri_ed_analysis.R"]
command_file <- grep("^--file=", commandArgs(), value = TRUE)
project_dir <- if (length(script_files)) {
  dirname(normalizePath(tail(script_files, 1), mustWork = TRUE))
} else if (length(command_file) && sub("^--file=", "", command_file[1]) != "-") {
  dirname(normalizePath(sub("^--file=", "", command_file[1]), mustWork = TRUE))
} else {
  getwd() 
}

# ---- Settings ---------------------------------------------------------------
# Input files and fitted-model cache
file_name <- "Los Angeles"
model_prefix <- "Los_Angeles_discrete_none_HEAT_by_cluster"
input_dir <- file.path(project_dir, "Los_Angeles_discrete_none")
results_dir <- file.path(input_dir, "climate_zone")
model_cache_dir <- file.path(input_dir, "model_cache")
force_refit <- FALSE
fit_negative_binomial <- FALSE

# Nonlinear heat resolution and posterior AR sampling
acute_rw2_n <- 50L
ar_n_samples <- 2000L
ar_batch_size <- 50L
ar_seed <- 20260910L
source(file.path(project_dir, "shared_functions.R"), local = TRUE)
source(file.path(project_dir, "cluster_ar.R"), local = TRUE)

# ---- Data preparation -------------------------------------------------------
prepared <- load_heat_data(input_dir, file_name)
clustered_zctas <- prepared$clustered_zctas
cluster_labels <- prepared$cluster_labels
cluster_vuln <- prepared$cluster_vuln %>% mutate(cluster = as.integer(cluster))
reference_cluster <- prepared$reference_cluster
cluster_terms <- paste0("cluster", levels(prepared$daily_ed$cluster))

# ---- Linear model component -------------------------------------------------
linear_formula <- D2 ~ cluster + chronic_heat +
  # acute_heat:cluster +  # SES-cluster heat slopes disabled for this fit.
  month + day_of_week

# ---- Existing nonlinear seasonal component ----------------------------------
seasonal_formula <- stats::update(linear_formula,
  . ~ . + f(doy, model = "rw2", scale.model = TRUE, constr = TRUE))

# ---- Nonlinear heat: f_heat,cz(x) = f_shared(x) + f_dev,cz(x) -------------------
# x is the existing centered acute_heat exposure. All ZIPs use f_shared;
# ZIPs without climate-zone assignments have no climate-zone deviation.
# Socioeconomic clusters span multiple climate zones, so the deviation uses
# each ZIP's climate zone; socioeconomic cluster terms remain in the linear part.
cz_lookup <- read.csv(file.path(project_dir, "datasets", "cz_zip.csv"),
                      colClasses = "character")
stopifnot(all(c("zip", "cz") %in% names(cz_lookup)), !anyDuplicated(cz_lookup$zip))
prepared$daily_ed$climate_zone <- cz_lookup$cz[
  match(as.character(prepared$daily_ed$zip), cz_lookup$zip)]
missing_cz <- is.na(prepared$daily_ed$climate_zone) |
  !nzchar(prepared$daily_ed$climate_zone)
prepared$daily_ed$climate_zone[missing_cz] <- NA_character_
message("ED ZIPs missing climate-zone assignments (retained with the shared RW2 curve only): ",
        paste(sort(unique(prepared$daily_ed$zip[missing_cz])), collapse = ", "))

# SES-cluster-month mean heat, expressed on each ZIP's overall centered heat scale.
acute_cf <- cluster_heat_reference(prepared$daily_ed)
acute_rw2_grid <- seq(min(c(prepared$daily_ed$acute_heat, acute_cf)),
  max(c(prepared$daily_ed$acute_heat, acute_cf)), length.out = acute_rw2_n)
heat_grid_index <- function(x) {
  midpoints <- (head(acute_rw2_grid, -1L) + tail(acute_rw2_grid, -1L)) / 2
  findInterval(x, midpoints) + 1L
}

# Shared RW2 curve: an index is present for every observation.
shared_rw2_field <- "acute_heat_rw2_shared"
prepared$daily_ed[[shared_rw2_field]] <- heat_grid_index(prepared$daily_ed$acute_heat)

# Climate-zone RW2 deviations: NA indices contribute zero for missing zones.
climate_zones <- sort(unique(prepared$daily_ed$climate_zone[!missing_cz]))
deviation_rw2_field <- "acute_heat_rw2_dev"
prepared$daily_ed[[deviation_rw2_field]] <- ifelse(missing_cz,
  NA_integer_, prepared$daily_ed[[shared_rw2_field]])
prepared$daily_ed$climate_zone_rep <- match(prepared$daily_ed$climate_zone, climate_zones)
# The placeholder replicate is unused when the deviation index is NA.
prepared$daily_ed$climate_zone_rep[missing_cz] <- 1L
prepared$daily_ed$acute_heat_rw2_cf <- heat_grid_index(acute_cf)

# All zone deviations share one precision through replicate; the shared curve
# has a separate precision. Retain only the sum-to-zero constraint per curve.
# Linear trends are allowed and are not penalized by the RW2 curvature prior.
# Full-response contrasts include both the fixed slopes and the entire curves.
rw2_term <- function(field) sprintf(
  paste0("f(%s, model = 'rw2', n = %dL, values = seq_len(%dL), ",
    "scale.model = TRUE, constr = TRUE)"),
  field, acute_rw2_n, acute_rw2_n)
deviation_term <- sub("model = 'rw2',", sprintf(
  "model = 'rw2', replicate = climate_zone_rep, nrep = %dL,",
  length(climate_zones)), rw2_term(deviation_rw2_field), fixed = TRUE)
rw2_terms <- c(rw2_term(shared_rw2_field), deviation_term)

# ---- Full formula used by both likelihoods -----------------------------------
nonlinear_formula <- stats::update(seasonal_formula,
  paste(". ~ . +", paste(rw2_terms, collapse = " + ")))
nonlinear_heat <- list(
  shared_field = shared_rw2_field, deviation_field = deviation_rw2_field,
  zone_levels = climate_zones, pooled = TRUE,
  heat_constraints = "sum-to-zero only", ar_reference = "SES-cluster-month mean heat",
  ses_heat_slopes = FALSE,
  common_heat_slope = FALSE,
  grid = acute_rw2_grid, zone = prepared$daily_ed$climate_zone,
  observed_index = prepared$daily_ed[[shared_rw2_field]],
  counterfactual_index = prepared$daily_ed$acute_heat_rw2_cf)
daily_ed <- prepared$daily_ed

# ---- Cluster profiles for reporting ------------------------------------------

vars <- c("pct_white", "pct_black", "pct_asian", "pct_hispanic",
          "poverty_rate", "renter_burden_rate", "housing_overcrowding",
          "median_income", "education_low", "pct_age17", "no_phone_rate",
          "pct_age65", "uninsured_rate", "language_isolation", "disability_rate",
          "outdoor_worker_rate", "unemployment_rate", "elderly_alone_rate",
          "alt_housing_rate", "old_housing_rate")
cluster_profiles <- clustered_zctas %>% group_by(cluster) %>%
  summarise(across(all_of(vars), ~ median(.x, na.rm = TRUE)), n = n())

# ------ POISSON ----------------------------------------------------------

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
capture.output({
  cat("====================================================================\n")
  cat("POISSON INLA - SUPPRESSED 2018 ED DATA\n")
  cat("INLA version:\n")
  print(packageVersion("INLA"))
  cat("\nR version:\n")
  print(R.version.string)
  cat("\nReference cluster (lowest vulnerability):", reference_cluster, "\n")
  cat("====================================================================\n\n")
  
  poisson_inla <- fit_heat_model(prepared, model_cache_dir, family = "poisson",
    force_refit = force_refit, model_formula = nonlinear_formula,
    cache_tag = "_shared_cz_rw2_pooled_no_fixed_heat", nonlinear_heat = nonlinear_heat)
  poisson_inla$nonlinear_heat <- nonlinear_heat
  
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
  
  saveRDS(poisson_p_table, paste0(results_dir, "/",
                                  model_prefix, "_suppressed_INLA_poisson_coef_table.rds"))  
}, file = paste0(results_dir, "/",
                 model_prefix, "_suppressed_ED_INLA_poisson_summary.txt"))
saveRDS(poisson_inla, paste0(results_dir, "/",
                             model_prefix, "_suppressed_INLA_poisson_model.rds"))
rm(poisson_inla)
gc()


# ------- NEGATIVE BINOMIAL ---------------------------------------------------

if (fit_negative_binomial) {
capture.output({
  cat("====================================================================\n")
  cat("NEGATIVE BINOMIAL INLA - SUPPRESSED 2018 ED DATA\n")
  cat("INLA version:\n")
  print(packageVersion("INLA"))
  cat("\nR version:\n")
  print(R.version.string)
  cat("\nReference cluster (lowest vulnerability):", reference_cluster, "\n")
  cat("====================================================================\n\n")
  
  nb_inla <- fit_heat_model(prepared, model_cache_dir, family = "nbinomial",
    force_refit = force_refit, model_formula = nonlinear_formula,
    cache_tag = "_shared_cz_rw2_pooled_no_fixed_heat", nonlinear_heat = nonlinear_heat)
  nb_inla$nonlinear_heat <- nonlinear_heat
  
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
  
  saveRDS(nb_p_table, paste0(results_dir, "/",
                             file_name, "_suppressed_INLA_nb_coef_table.rds"))
}, file = paste0(results_dir, "/",
                 file_name, "_suppressed_ED_INLA_nb_summary.txt"))
saveRDS(nb_inla, paste0(results_dir, "/",
                        model_prefix, "_suppressed_INLA_nb_model.rds"))
rm(nb_inla, expected_counts)
}
gc()

# ===========================================================================
# extract INLA coefficients and calculate IRRs
# ===========================================================================

poisson_inla <- readRDS(paste0(results_dir, "/",
                               model_prefix, "_suppressed_INLA_poisson_model.rds"))
p_table <- as.data.frame(poisson_inla$summary.fixed)
p_table$term <- rownames(p_table)
p_table$p_IRR_gt_1 <- vapply(poisson_inla$marginals.fixed[p_table$term],
                             function(marginal) 1 - inla.pmarginal(0, marginal), numeric(1))

# Cluster coefficients describe the linear intercept component; full comparisons
# also include the shared curve and any climate-zone deviation.
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

write.csv(final_table, paste0(results_dir, "/",
                              file_name, "_poisson_cluster_results.csv"), row.names = FALSE)

# ===========================================================================
# POSTERIOR CLUSTER-AVERAGE ATTRIBUTABLE RISK
# ===========================================================================

# Model coefficients: cluster-specific intercepts and slopes
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

heat_slopes <- if ("acute_heat" %in% poisson_inla$names.fixed) {
  summarise_heat_slopes(poisson_inla, levels(daily_ed$cluster)) %>%
    transmute(cluster = paste0("cluster", cluster), beta = beta_inla_mean, P_beta_gt_0)
} else tibble(cluster = paste0("cluster", levels(daily_ed$cluster)))
# These beta summaries cover only the linear component, not the full nonlinear heat response.
# Equal-weight average of signed ZIP-day rate contrasts within each draw.
# Sampling workers read the saved fit. Release the large parent copy first.
ar_sampling_info <- list(n_observations = nrow(poisson_inla$summary.linear.predictor),
  names.fixed = poisson_inla$names.fixed)
ar_nonlinear_heat <- poisson_inla$nonlinear_heat
stopifnot(isTRUE(all.equal(as.numeric(poisson_inla$.args$E), as.numeric(daily_ed$Population))),
  all(is.finite(poisson_inla$summary.fixed$mean)),
  all(is.finite(poisson_inla$summary.fixed$sd)),
  all(is.finite(poisson_inla$summary.hyperpar$mean)))
message("Fitted model mode status: ", poisson_inla$mode$mode.status)
rm(poisson_inla)
gc(verbose = FALSE)
cluster_posterior <- cluster_ar_exceedance(
  ar_sampling_info, daily_ed, reference_cluster,
  n_samples = ar_n_samples, batch_size = ar_batch_size, seed = ar_seed,
  nonlinear_heat = ar_nonlinear_heat,
  sampling_model_file = file.path(results_dir, paste0(model_prefix,
    "_suppressed_INLA_poisson_model.rds")))
saveRDS(cluster_posterior, paste0(results_dir, "/",
                                  file_name, "_poisson_cluster_AR_posterior.rds"))

cluster_results <- cluster_intercepts %>%
  left_join(heat_slopes, by = "cluster") %>%
  left_join(cluster_posterior$summary, by = "cluster") %>%
  arrange(as.numeric(sub("cluster", "", cluster)))

write.csv(cluster_results, paste0(results_dir, "/",
                                  file_name, "_poisson_cluster_AR_exceedance.csv"), row.names = FALSE)

print(cluster_results)
gc()
