# Shared calculations for hri_ed_analysis.R and acute_heat.R.
# Sourcing this file only defines functions; it does not load data or fit models.

# SVI-informed composite with the original ranks, domains and missing-value rule.
compute_vulnerability_score <- function(df) {
  df |>
    dplyr::mutate(
      socioeconomic = rowMeans(cbind(
        dplyr::percent_rank(poverty_rate),
        1 - dplyr::percent_rank(median_income),
        dplyr::percent_rank(unemployment_rate),
        dplyr::percent_rank(education_low),
        dplyr::percent_rank(uninsured_rate)), na.rm = TRUE),
      household = rowMeans(cbind(
        dplyr::percent_rank(pct_age65),
        dplyr::percent_rank(pct_age17),
        dplyr::percent_rank(disability_rate),
        dplyr::percent_rank(language_isolation),
        dplyr::percent_rank(no_phone_rate),
        dplyr::percent_rank(elderly_alone_rate)), na.rm = TRUE),
      minority_status = dplyr::percent_rank(1 - pct_white),
      housing = rowMeans(cbind(
        dplyr::percent_rank(housing_overcrowding),
        dplyr::percent_rank(renter_burden_rate),
        dplyr::percent_rank(alt_housing_rate),
        dplyr::percent_rank(old_housing_rate)), na.rm = TRUE),
      occupational_exposure = dplyr::percent_rank(outdoor_worker_rate),
      vulnerability_score = rowMeans(cbind(
        socioeconomic, household, minority_status, housing,
        occupational_exposure), na.rm = TRUE))
}

prepare_heat_data <- function(clustered_zctas, daily_ed) {
  if (inherits(clustered_zctas, "sf")) {
    clustered_zctas <- sf::st_drop_geometry(clustered_zctas)
  }
  clustered_zctas <- clustered_zctas |>
    dplyr::mutate(GEOID = as.character(GEOID), cluster = as.character(cluster))
  if (anyNA(clustered_zctas$GEOID) || anyNA(clustered_zctas$cluster) ||
      anyDuplicated(clustered_zctas$GEOID)) {
    stop("Cluster input requires unique ZCTA identifiers and nonmissing clusters.")
  }
  cluster_vuln <- clustered_zctas |>
    compute_vulnerability_score() |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(vulnerability_score = mean(vulnerability_score, na.rm = TRUE),
                     n = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(vulnerability_score)
  if (any(!is.finite(cluster_vuln$vulnerability_score))) {
    stop("Undefined cluster vulnerability.")
  }
  reference_cluster <- cluster_vuln$cluster[1]
  cluster_labels <- dplyr::distinct(clustered_zctas, cluster, cluster_label)
  if (anyDuplicated(cluster_labels$cluster)) stop("Multiple labels for one cluster.")

  daily_ed <- daily_ed |>
    dplyr::mutate(cluster = clustered_zctas$cluster[
      match(as.character(zip), clustered_zctas$GEOID)])
  heat_mean <- mean(daily_ed$Max_HI_Value, na.rm = TRUE)
  heat_sd <- sd(daily_ed$Max_HI_Value, na.rm = TRUE)
  if (!is.finite(heat_sd) || heat_sd <= 0) stop("Heat must have positive finite SD.")
  daily_ed <- daily_ed |>
    dplyr::mutate(heat_z = (Max_HI_Value - heat_mean) / heat_sd) |>
    dplyr::group_by(zip) |>
    dplyr::mutate(chronic_heat = mean(heat_z, na.rm = TRUE),
                   acute_heat = heat_z - chronic_heat) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      day_of_week = factor(lubridate::wday(date, label = TRUE), ordered = FALSE),
      doy = as.integer(doy), month = factor(lubridate::month(date)),
      cluster = stats::relevel(factor(cluster), ref = reference_cluster)) |>
    dplyr::select(zip, date, D2, Population, cluster, heat_z, chronic_heat,
                   acute_heat, month, day_of_week, doy)
  numeric_inputs <- c("D2", "Population", "heat_z", "chronic_heat", "acute_heat", "doy")
  if (anyNA(daily_ed) || any(!is.finite(as.matrix(daily_ed[numeric_inputs]))) ||
      any(daily_ed$Population <= 0) || any(daily_ed$D2 < 0) ||
      any(daily_ed$D2 != floor(daily_ed$D2))) {
    stop("Invalid model rows; resolve missing/nonfinite values, population or counts first.")
  }
  list(daily_ed = daily_ed, clustered_zctas = clustered_zctas,
       cluster_vuln = cluster_vuln, cluster_labels = cluster_labels,
       reference_cluster = reference_cluster, heat_mean = heat_mean, heat_sd = heat_sd)
}

load_heat_data <- function(input_dir, study_name = "Los Angeles") {
  prepare_heat_data(
    readRDS(file.path(input_dir, paste0(study_name, "_ACS_zcta_clustered_for_ED.rds"))),
    readRDS(file.path(input_dir, paste0(study_name, "_daily_ed_no_heat.rds"))))
}

# Keep one definition of the predictor for both likelihoods and both scripts.
heat_model_formula <- function() {
  D2 ~ cluster + chronic_heat + acute_heat + acute_heat:cluster + month +
    day_of_week + f(doy, model = "rw2", scale.model = TRUE, constr = TRUE)
}

make_heat_lincombs <- function(data) {
  cluster_ids <- levels(data$cluster)
  reference_cluster <- cluster_ids[1]
  # Remove only the random-effect term to get the fixed-effect coefficient names.
  fixed_formula <- stats::update(heat_model_formula(),
    ~ . - f(doy, model = "rw2", scale.model = TRUE, constr = TRUE))
  design_names <- colnames(stats::model.matrix(fixed_formula, data = data))
  setNames(lapply(cluster_ids, function(k) {
    weights <- list(acute_heat = 1)
    if (k != reference_cluster) {
      candidates <- c(paste0("cluster", k, ":acute_heat"), paste0("acute_heat:cluster", k))
      interaction_name <- intersect(candidates, design_names)
      if (length(interaction_name) != 1L) stop("Cannot identify interaction for cluster ", k)
      weights[[interaction_name]] <- 1
    }
    do.call(INLA::inla.make.lincomb, weights)[[1]]
  }), paste0("cluster", cluster_ids))
}

# Hash the actual ordered model inputs rather than input filenames. The same
# prepared data and model specification therefore reuse the same fit in either script.
heat_object_md5 <- function(object) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(object, path, compress = FALSE)
  unname(tools::md5sum(path))
}

fit_heat_model <- function(prepared, cache_dir, family = c("poisson", "nbinomial"),
                           force_refit = FALSE) {
  family <- match.arg(family)
  data <- prepared$daily_ed
  model_formula <- heat_model_formula()
  heat_lincombs <- make_heat_lincombs(data)
  compute_settings <- list(dic = TRUE, waic = TRUE, config = TRUE)
  predictor_settings <- list(compute = TRUE)
  signature <- heat_object_md5(list(
    data = data, formula = deparse(model_formula), family = family,
    lincombs = heat_lincombs, contrasts = getOption("contrasts"),
    control.compute = compute_settings, control.predictor = predictor_settings,
    INLA_version = as.character(utils::packageVersion("INLA")), cache_version = 1L))
  cache_file <- file.path(cache_dir, paste0("heat_", family, "_fit.rds"))
  cached <- if (!force_refit && file.exists(cache_file)) readRDS(cache_file) else NULL
  cluster_terms <- names(heat_lincombs)
  if (!is.null(cached) && identical(cached$signature, signature) &&
      all(cluster_terms %in% names(cached$fit$marginals.lincomb.derived))) {
    message("Using shared cached ", family, " heat model.")
    return(cached$fit)
  }
  message("Fitting ", family, " heat model with cluster-specific linear combinations...")
  fit <- INLA::inla(
    model_formula, family = family, data = data, E = data$Population,
    lincomb = heat_lincombs, control.compute = compute_settings,
    control.predictor = predictor_settings)
  if (!all(cluster_terms %in% names(fit$marginals.lincomb.derived))) {
    stop("INLA did not return the requested cluster-slope marginals.")
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(signature = signature, fit = fit), cache_file)
  fit
}

summarise_heat_slopes <- function(model, cluster_ids) {
  cluster_terms <- paste0("cluster", cluster_ids)
  if (!all(cluster_terms %in% names(model$marginals.lincomb.derived))) {
    stop("The model lacks cluster-slope marginals; use fit_heat_model() to obtain them.")
  }
  result <- dplyr::bind_rows(lapply(seq_along(cluster_ids), function(j) {
    term <- cluster_terms[j]
    marginal <- model$marginals.lincomb.derived[[term]]
    beta_q <- INLA::inla.qmarginal(c(0.025, 0.5, 0.975), marginal)
    data.frame(cluster = cluster_ids[j],
      beta_mean = INLA::inla.emarginal(identity, marginal),
      # Retain INLA's stored mean as well, for the original HRI output convention.
      beta_inla_mean = model$summary.lincomb.derived[term, "mean"],
      beta_median = beta_q[2], beta_low = beta_q[1], beta_high = beta_q[3],
      P_beta_gt_0 = 1 - INLA::inla.pmarginal(0, marginal),
      RR_median = exp(beta_q[2]), RR_low = exp(beta_q[1]), RR_high = exp(beta_q[3]))
  }))
  stopifnot(all(is.finite(result$P_beta_gt_0)),
            all(result$P_beta_gt_0 >= 0 & result$P_beta_gt_0 <= 1),
            all(is.finite(result$RR_high)), all(result$RR_low > 0))
  result
}
