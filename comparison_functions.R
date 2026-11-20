# Optional comparisons invoked through hri_ed_analysis.R.

protected_output_hashes <- function(output_dir) {
  paths <- c(list.files("Los_Angeles_discrete_none/climate_zone_added_nb_pc", full.names = TRUE),
    list.files("maps/nb_pc_prior", full.names = TRUE, recursive = TRUE),
    list.files("methods", full.names = TRUE, recursive = TRUE))
  paths <- sort(unique(normalizePath(paths[!dir.exists(paths)], winslash = "/", mustWork = TRUE)))
  hashes <- tools::md5sum(paths)
  file <- file.path(output_dir, "protected_files.rds")
  if (file.exists(file)) stopifnot(identical(hashes, readRDS(file))) else saveRDS(hashes, file)
  hashes
}

run_nb_prior_sensitivity <- function(output_dir = "prior_sensitivity_nb", restart_scenarios = character()) {
  # Vary the four precision priors jointly; preserve every other model choice.
  # Fit-only comparison: no excess-rate posterior sampling.
  # The saved current fit is the exact baseline and is reused without changing it.
  scenarios <- list(
    pc_baseline = list(prior = "pc.prec", param = c(1, 0.01)),
    pc_stronger = list(prior = "pc.prec", param = c(0.5, 0.01)),
    pc_weaker = list(prior = "pc.prec", param = c(2, 0.01)),
    gamma_comparison = list(prior = "loggamma", param = c(1, 0.00005)))
  stopifnot(all(restart_scenarios %in% names(scenarios)))
  baseline_dir <- "Los_Angeles_discrete_none/climate_zone_added_nb_pc"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  baseline_output <- file.path(output_dir, "pc_baseline")
  dir.create(baseline_output, showWarnings = FALSE)
  for (name in c("nbinomial_model.rds", "signature.txt")) {
    from <- file.path(baseline_dir, name)
    to <- file.path(baseline_output, name)
    stopifnot(file.exists(from))
    if (!file.exists(to)) stopifnot(file.copy(from, to))
    stopifnot(identical(unname(tools::md5sum(from)), unname(tools::md5sum(to))))
  }
  summaries <- invariants <- list()
  for (scenario in names(scenarios)) {
    message("CORRECTED MODEL PRIOR SENSITIVITY: ", scenario)
    scenario_dir <- file.path(output_dir, scenario)
    summaries[[scenario]] <- fit_preserved_nb(output_dir = scenario_dir, fit_only = TRUE,
      save_config = identical(scenario, "pc_baseline"),
      restart_fit = scenario %in% restart_scenarios,
      precision_prior = scenarios[[scenario]])
    validate_nb_pc_model(results_dir = scenario_dir, with_maps = FALSE,
      with_excess_rates = FALSE, check_sampling = FALSE, precision_prior = scenarios[[scenario]])
    fit <- readRDS(file.path(scenario_dir, "nbinomial_model.rds"))
    size <- fit$all.hyper$family[[1L]]$hyper[[1L]]
    invariants[[scenario]] <- list(data = heat_object_md5(fit$.args$data),
      formula = deparse(fit$.args$formula), offset = heat_object_md5(fit$.args$E),
      fixed_design = fit$preserved_nb_specification$fixed_design_md5,
      grid = fit$nonlinear_heat$grid, correlation = fit$preserved_nb_specification$correlation,
      fixed_priors = fit$.args$control.fixed,
      # INLA's expanded control.family contains closures with process-specific
      # environments. Compare the actual fitted prior, parameters and fixed flag.
      size_prior = list(prior = as.character(size$prior), param = as.numeric(size$param),
        fixed = as.logical(size$fixed)),
      full_contrasts = fit$full_heat_irr$combinations)
    stopifnot(identical(invariants[[scenario]], invariants[["pc_baseline"]]))
    rm(fit); gc(verbose = FALSE)
  }
  bind_scenarios <- function(values) dplyr::bind_rows(lapply(names(values), function(scenario)
    data.frame(scenario = scenario, values[[scenario]], row.names = NULL)))
  metrics <- bind_scenarios(lapply(summaries, `[[`, "metrics"))
  heat <- bind_scenarios(lapply(summaries, `[[`, "heat"))
  clusters <- bind_scenarios(lapply(summaries, `[[`, "clusters"))
  baseline <- summaries$pc_baseline
  metrics$delta_WAIC_vs_baseline <- metrics$WAIC - baseline$metrics$WAIC
  metrics$delta_DIC_vs_baseline <- metrics$DIC - baseline$metrics$DIC
  metrics$delta_WAIC_vs_gamma <- metrics$WAIC - summaries$gamma_comparison$metrics$WAIC
  metrics$max_heat_IRR_change_pct_vs_baseline <- vapply(names(scenarios), function(scenario) {
    h <- summaries[[scenario]]$heat
    idx <- match(h$lincomb_id, baseline$heat$lincomb_id)
    stopifnot(!anyNA(idx), nrow(h) == 42L)
    100 * max(abs(h$IRR_median / baseline$heat$IRR_median[idx] - 1))
  }, numeric(1))
  metrics$max_heat_probability_change_pp_vs_baseline <- vapply(names(scenarios), function(scenario) {
    h <- summaries[[scenario]]$heat
    idx <- match(h$lincomb_id, baseline$heat$lincomb_id)
    100 * max(abs(h$P_IRR_gt_1 - baseline$heat$P_IRR_gt_1[idx]))
  }, numeric(1))
  for (scenario in names(scenarios)) {
    h <- summaries[[scenario]]$heat
    writeLines(c(paste("CORRECTED MODEL PRIOR SENSITIVITY:", scenario),
      readLines(file.path(output_dir, scenario, "priors.txt")), "",
      "Conditional cluster intercept IRRs:", capture.output(print(summaries[[scenario]]$clusters, row.names = FALSE)),
      "", "Full one-SD heat contrasts by SES cluster and climate zone:",
      capture.output(print(h, row.names = FALSE))),
      file.path(output_dir, paste0(scenario, "_prior_results.txt")))
  }
  write.csv(metrics, file.path(output_dir, "comparison.csv"), row.names = FALSE)
  write.csv(heat, file.path(output_dir, "heat_contrasts.csv"), row.names = FALSE)
  write.csv(clusters, file.path(output_dir, "conditional_cluster_irrs.csv"), row.names = FALSE)
  report <- c("CORRECTED NEGATIVE-BINOMIAL PRIOR SENSITIVITY", "",
    "All four scenarios exclude the common heat slope and all 13 SES heat interactions.",
    "Identical 54,910 observations, 309 ZIPs, offsets, 50-node grid, fixed priors, NB size prior and PC correlation prior.",
    "Only the four random-effect precision priors vary, jointly:",
    "pc_baseline: P(SD > 1) = 0.01; exact current fit reused.",
    "pc_stronger: P(SD > 0.5) = 0.01.", "pc_weaker: P(SD > 2) = 0.01.",
    "gamma_comparison: Gamma(shape=1, rate=0.00005) on precision; comparator only.",
    "Every scenario retains the same PC zone-correlation prior and PC NB-size prior.",
    "The main corrected model and current maps are unchanged.", "",
    "Smaller WAIC/DIC is preferred by those criteria. Changes compare posterior point summaries.",
    "Probability changes are percentage points, using INLA full-contrast marginals.",
    "This comparison is fit-only: no excess-rate posterior draws are generated.",
    "The main baseline's excess-rate draws are preserved; no draws are created or used in these comparisons.",
    "RW2 precision controls curvature, not unpenalized linear trends; interpret complete heat contrasts.",
    "These results assess the listed joint precision-prior changes, not every possible prior choice.", "",
    capture.output(print(metrics, row.names = FALSE, digits = 6)), "",
    "Detailed outputs: heat_contrasts.csv (42 contrasts per scenario) and conditional_cluster_irrs.csv.",
    "Each scenario contains the saved fit, active fitted priors, summaries and validation checks.")
  writeLines(report, file.path(output_dir, "comparison.txt"))
  writeLines(c(report, unlist(lapply(names(scenarios), function(scenario)
    c("", readLines(file.path(output_dir, paste0(scenario, "_prior_results.txt"))))))),
    file.path(output_dir, "prior_results.txt"))
  jsonlite::write_json(list(status = "completed", scenarios = scenarios,
    invariants_identical = TRUE, observations = 54910L, ZIPs = 309L, fixed_coefficients = 26L,
    SES_heat_slopes = 0L, full_heat_contrasts_per_scenario = 42L,
    posterior_sampling = FALSE, newly_generated_excess_draws = 0L, baseline_reused = TRUE,
    current_model_md5 = unname(tools::md5sum(file.path(baseline_dir, "nbinomial_model.rds")))),
    file.path(output_dir, "sensitivity_validation.json"), pretty = TRUE, auto_unbox = TRUE, digits = 16)
  message("All four corrected-model prior sensitivity comparisons completed.")
  invisible(metrics)
}



likelihood_zero_diagnostics <- function(fit, data, family) {
  stopifnot(family %in% c("poisson", "nbinomial"))
  count_mean <- as.numeric(fit$summary.fitted.values$mean) * data$Population
  stopifnot(length(count_mean) == nrow(data), all(is.finite(count_mean)), all(count_mean > 0),
    all(is.finite(data$D2)), all(data$D2 >= 0), all(data$D2 == floor(data$D2)))
  size <- if (family == "nbinomial") {
    index <- grep("^size for", rownames(fit$summary.hyperpar))
    stopifnot(length(index) == 1L)
    fit$summary.hyperpar[index, "mean"]
  } else NA_real_
  probability <- function(q) if (family == "poisson")
    stats::ppois(q, lambda = count_mean) else stats::pnbinom(q, mu = count_mean, size = size)
  p_zero <- probability(0)
  rows <- data.frame(family = family, zip = as.character(data$zip),
    cluster = as.character(data$cluster),
    climate_zone = as.character(data$climate_zone_group), month = as.character(data$month),
    observed_count = data$D2, fitted_mean_count = count_mean, predicted_zero_probability_plugin = p_zero)
  # The group index for missing climate zones is a computational placeholder;
  # the missing deviation index identifies them and must not be counted as zone 1.
  if (!is.null(fit$nonlinear_heat$zone)) rows$climate_zone <- as.character(fit$nonlinear_heat$zone)
  rows$climate_zone[is.na(rows$climate_zone)] <- "missing"
  summarize <- function(indices) data.frame(n = length(indices),
    observed_zeros = sum(data$D2[indices] == 0), expected_zeros_plugin = sum(p_zero[indices]),
    observed_zero_fraction = mean(data$D2[indices] == 0), predicted_zero_fraction_plugin = mean(p_zero[indices]),
    zero_gap_percentage_points = 100 * (mean(data$D2[indices] == 0) - mean(p_zero[indices])))
  grouped <- do.call(rbind, lapply(c("cluster", "zip", "climate_zone", "month"), function(group) {
    groups <- split(seq_len(nrow(data)), rows[[group]])
    do.call(rbind, lapply(names(groups), function(value)
      data.frame(family = family, group = group, level = value, summarize(groups[[value]]), row.names = NULL)))
  }))
  cumulative <- lapply(c(0, 1, 2, 4), probability)
  probabilities <- cbind(cumulative[[1]], cumulative[[2]] - cumulative[[1]],
    cumulative[[3]] - cumulative[[2]], cumulative[[4]] - cumulative[[3]], 1 - cumulative[[4]])
  observed_bin <- cut(data$D2, breaks = c(-Inf, 0, 1, 2, 4, Inf),
    labels = c("0", "1", "2", "3-4", "5+"))
  count_bins <- data.frame(family = family, count_bin = levels(observed_bin),
    observed_count = as.numeric(table(observed_bin)), expected_count_plugin = colSums(probabilities))
  count_bins$observed_fraction <- count_bins$observed_count / nrow(data)
  count_bins$predicted_fraction_plugin <- count_bins$expected_count_plugin / nrow(data)
  variance <- if (family == "poisson") count_mean else count_mean + count_mean^2 / size
  metrics <- data.frame(family = family, DIC = fit$dic$dic, WAIC = fit$waic$waic,
    summarize(seq_len(nrow(data))), NB_size_posterior_mean = size,
    mean_squared_Pearson_residual_plugin = mean((data$D2 - count_mean)^2 / variance))
  stopifnot(all(is.finite(p_zero)), all(p_zero >= 0 & p_zero <= 1),
    abs(sum(count_bins$expected_count_plugin) - nrow(data)) < 1e-6)
  list(metrics = metrics, grouped = grouped, count_bins = count_bins, rows = rows)
}

likelihood_common_priors <- function(fit) {
  hyper <- function(x) list(prior = as.character(x$prior), param = as.numeric(x$param),
    fixed = as.logical(x$fixed))
  list(fixed = fit$all.hyper$fixed, random = lapply(fit$all.hyper$random, function(term)
    list(term = term$hyperid, precision = hyper(term$hyper[[1L]]),
      correlation = if (term$hyperid == "acute_heat_rw2_dev") hyper(term$group.hyper[[1L]]) else NULL)))
}

baseline_comparison_inputs <- function() {
  baseline_file <- "Los_Angeles_discrete_none/climate_zone_added_nb_pc/nbinomial_model.rds"
  fit <- readRDS(baseline_file)
  saved <- readRDS("Los_Angeles_discrete_none/climate_zone_added_nb_pc/summaries.rds")
  stopifnot(identical(fit$.args$family, "nbinomial"), fit$mode$mode.status == 0L,
    length(fit$names.fixed) == 26L, nrow(fit$.args$data) == 54910L,
    length(unique(fit$.args$data$zip)) == 309L,
    identical(as.numeric(fit$.args$E), as.numeric(fit$.args$data$Population)),
    identical(saved$metrics$DIC, fit$dic$dic), identical(saved$metrics$WAIC, fit$waic$waic))
  precision <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
  correlation <- exchangeable_pc_prior(5)$hyper
  for (term in fit$all.hyper$fixed) stopifnot(term$prior.mean == 0,
    term$prior.prec == if (term$label == "(Intercept)") 0 else 0.001)
  for (term in fit$all.hyper$random) {
    stopifnot(identical(as.character(term$hyper[[1L]]$prior), precision$prec$prior),
      identical(as.numeric(term$hyper[[1L]]$param), precision$prec$param))
    if (term$hyperid == "acute_heat_rw2_dev") stopifnot(
      identical(as.character(term$group.hyper[[1L]]$prior), correlation$rho$prior))
  }
  size <- fit$all.hyper$family[[1L]]$hyper[[1L]]
  size <- list(prior = as.character(size$prior), param = as.numeric(size$param), fixed = as.logical(size$fixed))
  stopifnot(identical(size, list(prior = "pc.mgamma", param = 7, fixed = FALSE)))
  inputs <- list(data = fit$.args$data, exposure = fit$.args$E,
    formula_text = deparse(fit$.args$formula), fixed_prior = fit$.args$control.fixed,
    precision_hyper = precision, zone_correlation_hyper = correlation,
    common_priors = likelihood_common_priors(fit),
    nonlinear_heat = fit$nonlinear_heat, names_fixed = fit$names.fixed,
    baseline_md5 = unname(tools::md5sum(baseline_file)), INLA_version = as.character(packageVersion("INLA")))
  list(inputs = inputs, reference = list(zero = likelihood_zero_diagnostics(fit, inputs$data, "nbinomial")))
}

run_likelihood_comparison <- function(output_dir = "Los_Angeles_discrete_none/climate_zone_added_poisson_pc") {
  suppressPackageStartupMessages(library(INLA))
  source("shared_functions.R")
  baseline_file <- "Los_Angeles_discrete_none/climate_zone_added_nb_pc/nbinomial_model.rds"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  before <- protected_output_hashes(output_dir)
  protected <- names(before)
  input_file <- file.path(output_dir, "matched_inputs.rds")
  reference_file <- file.path(output_dir, "nb_baseline_diagnostics.rds")
  if (!file.exists(input_file) || !file.exists(reference_file)) {
    baseline <- baseline_comparison_inputs()
    saveRDS(baseline$inputs, input_file)
    saveRDS(baseline$reference$zero, reference_file)
    rm(baseline); gc(verbose = FALSE)
  }
  inputs <- readRDS(input_file)
  stopifnot(identical(inputs$baseline_md5, unname(tools::md5sum(baseline_file))),
    identical(inputs$INLA_version, as.character(packageVersion("INLA"))))
  model_file <- file.path(output_dir, "poisson_model.rds")
  if (!file.exists(model_file)) {
    message("Fitting matched Poisson model: same data, offsets, fixed effects and random-effect priors; no draws.")
    fit <- fit_inla_model(inputs, output_dir, family = "poisson")
    fit$nonlinear_heat <- inputs$nonlinear_heat
    fit$matched_likelihood_comparison <- list(baseline_md5 = inputs$baseline_md5,
      elapsed_seconds = fit$analysis_elapsed_seconds, excess_rate_draws = 0L)
    saveRDS(fit, model_file)
  } else fit <- readRDS(model_file)
  stopifnot(identical(fit$names.fixed, inputs$names_fixed),
    identical(fit$.args$data, inputs$data), identical(as.numeric(fit$.args$E), as.numeric(inputs$exposure)),
    identical(deparse(fit$.args$formula), inputs$formula_text),
    identical(likelihood_common_priors(fit), inputs$common_priors),
    identical(fit$matched_likelihood_comparison$baseline_md5, inputs$baseline_md5),
    fit$mode$mode.status == 0L)
  poisson <- likelihood_zero_diagnostics(fit, inputs$data, "poisson")
  saveRDS(poisson, file.path(output_dir, "poisson_diagnostics.rds"))
  write_fit_details(fit, output_dir, summary_file = "poisson_fit_summary.txt",
    prior_file = "poisson_fitted_priors.txt")
  rm(fit, inputs); gc(verbose = FALSE)
  nb <- readRDS(reference_file)
  comparison <- rbind(nb$metrics, poisson$metrics)
  comparison$delta_DIC_vs_NB <- comparison$DIC - nb$metrics$DIC
  comparison$delta_WAIC_vs_NB <- comparison$WAIC - nb$metrics$WAIC
  write.csv(comparison, file.path(output_dir, "comparison.csv"), row.names = FALSE)
  write.csv(rbind(nb$grouped, poisson$grouped), file.path(output_dir, "zero_diagnostics_by_group.csv"), row.names = FALSE)
  write.csv(rbind(nb$count_bins, poisson$count_bins), file.path(output_dir, "count_distribution_diagnostics.csv"), row.names = FALSE)
  write.csv(rbind(nb$rows, poisson$rows), file.path(output_dir, "zero_diagnostics_by_observation.csv"), row.names = FALSE)
  after <- tools::md5sum(protected)
  stopifnot(identical(before, after))
  jsonlite::write_json(list(status = "completed", observations = nrow(nb$rows),
    ZIPs = length(unique(nb$rows$zip)), fixed_coefficients = 26L,
    identical_predictor_data_exposure_and_common_priors = TRUE, baseline_reused = TRUE,
    protected_files_unchanged = TRUE, protected_file_count = length(before),
    baseline_model_md5 = unname(tools::md5sum(baseline_file)),
    poisson_model_md5 = unname(tools::md5sum(model_file)), newly_generated_excess_rate_draws = 0L,
    zero_diagnostic = "In-sample plug-in: posterior mean rate times population; NB posterior mean size. No joint-posterior integration or predictive p-values.",
    INLA_version = as.character(packageVersion("INLA"))),
    file.path(output_dir, "validation.json"), pretty = TRUE, auto_unbox = TRUE, digits = 16)
  writeLines(c("MATCHED POISSON VERSUS CURRENT BASELINE NEGATIVE BINOMIAL", "",
    "The exact current NB fit is reused. Only the count likelihood changes in the Poisson fit.",
    "Same 54,910 observations, offsets, 26 fixed coefficients, ZIP iid, shared/zone RW2 and exchangeable group prior.",
    "Fixed Normal priors and four PC precision priors are unchanged; Poisson has no NB size parameter.",
    "Both fits use INLA's default approximation/integration settings and two outer threads.",
    "The Poisson fit skips saved joint-sampling configurations and derived heat contrasts; these are not needed for the requested diagnostics.",
    "No excess-rate draws or map exports are generated. Baseline results, maps and methods files are hash-verified unchanged.", "",
    "Smaller DIC and WAIC indicate better fit by those criteria; differences are Poisson minus NB.",
    "Zero checks are in-sample plug-in diagnostics, not fully integrated posterior-predictive tests.",
    "Predicted counts use population times posterior mean fitted rate, and posterior mean size for NB.",
    "The zero gap is observed minus predicted in percentage points; positive means underprediction of zeros.",
    "The mean squared Pearson residual is descriptive and is not adjusted for fitted degrees of freedom.", "",
    capture.output(print(comparison, row.names = FALSE, digits = 8)), "",
    "Grouped zero diagnostics: zero_diagnostics_by_group.csv (SES cluster, ZIP, climate zone and month).",
    "Count-bin checks: count_distribution_diagnostics.csv (0, 1, 2, 3-4, 5+)."),
    file.path(output_dir, "comparison.txt"))
  print(comparison, row.names = FALSE)
  invisible(comparison)
}
