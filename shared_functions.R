# Single export standard for every geographic map in this project.
map_width_inches <- 10
map_height_inches <- 8
map_dpi <- 300

save_project_map <- function(filename, plot) {
  extension <- tools::file_ext(filename)
  temporary <- tempfile(pattern = "map_export_", tmpdir = dirname(filename),
                         fileext = paste0(".", extension))
  on.exit(unlink(temporary))
  ggplot2::ggsave(temporary, plot = plot, width = map_width_inches,
                  height = map_height_inches, units = "in", dpi = map_dpi,
                  bg = "white", device = if (tolower(extension) == "pdf")
                    grDevices::cairo_pdf else NULL)
  stopifnot(file.exists(temporary), file.info(temporary)$size > 0)
  # Replace the completed export atomically; overwriting an open destination
  # through file.copy() can fail on Windows even when replacement is allowed.
  stopifnot(file.rename(temporary, filename))
}

# Bound posterior-sampling memory by restarting an R worker periodically.
# The R RNG state is carried between workers; INLA's seed is set per batch.
heat_posterior_sampler <- function(model, model_file = NULL, draws_per_worker = 500L) {
  worker <- NULL
  worker_draws <- 0L
  setup_worker <- function(paths, path) {
    .libPaths(paths)
    assign(".heat_sampling_model", readRDS(path), envir = .GlobalEnv)
    loadNamespace("INLA")
    NULL
  }
  draw_worker <- function(n, selection, seed, rng) {
    assign(".Random.seed", rng, envir = .GlobalEnv)
    samples <- INLA::inla.posterior.sample(n, get(".heat_sampling_model", .GlobalEnv),
      selection = selection, seed = seed, num.threads = "1:1", add.names = FALSE)
    list(samples = samples, rng = get(".Random.seed", .GlobalEnv))
  }
  # Avoid serializing the enclosing model and socket connections with functions.
  environment(setup_worker) <- environment(draw_worker) <- baseenv()
  close_worker <- function() {
    if (!is.null(worker)) try(parallel::stopCluster(worker), silent = TRUE)
    worker <<- NULL
  }
  sample_batch <- function(n, selection, seed) {
    if (is.null(model_file)) {
      return(INLA::inla.posterior.sample(n, model, selection = selection,
        seed = seed, num.threads = "1:1", add.names = FALSE))
    }
    batch_rng <- get(".Random.seed", .GlobalEnv)
    if (is.null(worker) || worker_draws + n > draws_per_worker) {
      close_worker()
      message("Starting fresh posterior-sampling worker (up to ", draws_per_worker, " draws).")
      worker <<- parallel::makePSOCKcluster(1L, outfile = "")
      parallel::clusterCall(worker, setup_worker, .libPaths(),
        normalizePath(model_file, winslash = "/", mustWork = TRUE))
      worker_draws <<- 0L
    }
    result <- parallel::clusterCall(worker, draw_worker,
      n, selection, seed, batch_rng)[[1L]]
    assign(".Random.seed", result$rng, envir = .GlobalEnv)
    worker_draws <<- worker_draws + n
    result$samples
  }
  list(sample = sample_batch, close = close_worker)
}

# Shared calculations for hri_ed_analysis.R and acute_heat.R.
# Sourcing this file only defines functions; it does not load data or fit models.

# Fixed probability scale shared by all probability maps.
probability_colour_scale <- function(probabilities, name) {
  stopifnot(all(is.finite(probabilities)), all(probabilities >= 0 & probabilities <= 1))
  ggplot2::scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0.5, limits = c(0, 1),
    breaks = seq(0, 1, 0.2), labels = function(x) sprintf("%.1f", x), oob = scales::squish,
    na.value = "grey92", name = name,
    guide = ggplot2::guide_colorbar(nbin = 1024, barheight = grid::unit(80, "mm")))
}

missing_result_hatching <- function(polygons) {
  if (!nrow(polygons)) return(sf::st_sf(geometry = sf::st_sfc(crs = sf::st_crs(polygons))))
  projected <- sf::st_transform(sf::st_make_valid(polygons), 3310)
  bounds <- sf::st_bbox(projected)
  height <- unname(bounds["ymax"] - bounds["ymin"])
  starts <- seq(unname(bounds["xmin"]) - height, unname(bounds["xmax"]), by = 1500)
  lines <- sf::st_sfc(lapply(starts, function(x) sf::st_linestring(rbind(
    c(x, bounds["ymin"]), c(x + height, bounds["ymax"])))), crs = 3310)
  clipped <- suppressWarnings(sf::st_intersection(lines, sf::st_union(sf::st_geometry(projected))))
  clipped <- clipped[!sf::st_is_empty(clipped) & sf::st_dimension(clipped) == 1L]
  sf::st_sf(geometry = sf::st_transform(clipped, sf::st_crs(polygons)))
}

missing_result_hatch_layer <- function(geometry) {
  if (is.null(geometry$hatching) || !nrow(geometry$hatching)) return(NULL)
  ggplot2::geom_sf(data = geometry$hatching, colour = "grey72", linewidth = 0.18,
    inherit.aes = FALSE, show.legend = FALSE)
}

mainland_probability_geometry <- function(geometry, assignments) {
  parts <- sf::st_cast(sf::st_geometry(geometry$boundary), "POLYGON")
  mainland <- parts[which.max(sf::st_area(parts))]
  polygons <- suppressWarnings(sf::st_intersection(geometry$zctas, mainland))
  polygons <- polygons[sf::st_dimension(polygons) == 2L, ]
  polygons <- dplyr::mutate(polygons, GEOID = as.character(ZCTA5CE20)) |>
    dplyr::left_join(dplyr::select(assignments, GEOID, cluster), by = "GEOID")
  polygons <- dplyr::filter(polygons, !is.na(cluster))
  # Hatch the gray mainland area not covered by any modeled ZCTA.
  missing <- suppressWarnings(sf::st_difference(
    sf::st_transform(sf::st_sf(geometry = mainland), 3310),
    sf::st_union(sf::st_transform(sf::st_make_valid(polygons), 3310))))
  missing <- sf::st_transform(missing, sf::st_crs(polygons))
  list(boundary = sf::st_sf(geometry = mainland), polygons = polygons,
       hatching = missing_result_hatching(missing),
       clusters = dplyr::summarise(dplyr::group_by(polygons, cluster), .groups = "drop"))
}

draw_probability_map <- function(geometry, estimates, legend, subtitle, by = "cluster") {
  polygons <- dplyr::left_join(geometry$polygons, estimates, by = by)
  stopifnot(!anyNA(polygons$probability))
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = geometry$boundary, fill = "grey92", colour = "grey60", linewidth = 0.2) +
    missing_result_hatch_layer(geometry) +
    ggplot2::geom_sf(data = polygons, ggplot2::aes(fill = probability), colour = "white", linewidth = 0.10) +
    ggplot2::geom_sf(data = geometry$clusters, fill = NA, colour = "grey25", linewidth = 0.25) +
    probability_colour_scale(estimates$probability, legend) +
    ggplot2::coord_sf(datum = NA, expand = FALSE) +
    ggplot2::labs(title = "Los Angeles posterior exceedance probabilities", subtitle = subtitle, caption = NULL) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(legend.position = "right", plot.margin = ggplot2::margin(12, 12, 12, 12))
}

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

  # Zero is the source's unavailable-heat sentinel. Exclude missing exposure
  # days before standardization and fitting; retain true zero ED counts.
  daily_ed <- daily_ed |>
    dplyr::mutate(Max_HI_Value = dplyr::na_if(Max_HI_Value, 0))
  missing_heat <- sum(is.na(daily_ed$Max_HI_Value))
  message("Excluding ", missing_heat, " daily ED rows with unavailable heat.")
  daily_ed <- daily_ed |>
    dplyr::filter(!is.na(Max_HI_Value)) |>
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

heat_model_formula <- function() {
  D2 ~ cluster + chronic_heat + acute_heat + acute_heat:cluster + month +
    day_of_week + f(doy, model = "rw2", scale.model = TRUE, constr = TRUE)
}

heat_fixed_formula <- function(model_formula) {
  labels <- attr(stats::terms(model_formula), "term.labels")
  stats::reformulate(labels[!grepl("^f\\(", labels)], response = "D2",
    env = environment(model_formula))
}

make_heat_lincombs <- function(data, model_formula = heat_model_formula()) {
  cluster_ids <- levels(data$cluster)
  reference_cluster <- cluster_ids[1]
  # Remove only the random-effect term to get the fixed-effect coefficient names.
  fixed_formula <- heat_fixed_formula(model_formula)
  design_names <- colnames(stats::model.matrix(fixed_formula, data = data))
  if (!"acute_heat" %in% design_names) return(list())
  setNames(lapply(cluster_ids, function(k) {
    weights <- list(acute_heat = 1)
    if (k != reference_cluster) {
      candidates <- c(paste0("cluster", k, ":acute_heat"), paste0("acute_heat:cluster", k))
      interaction_name <- intersect(candidates, design_names)
      if (length(interaction_name) > 1L) stop("Ambiguous interaction for cluster ", k)
      if (length(interaction_name)) weights[[interaction_name]] <- 1
    }
    do.call(INLA::inla.make.lincomb, weights)[[1]]
  }), paste0("cluster", cluster_ids))
}

# Coefficients of the full log-IRR, including any RW2 linear trends.
full_heat_contrast_weights <- function(cluster, reference_cluster, climate_zone,
                                       nonlinear_heat, design_names,
                                       heat_from = 0, heat_to = 1) {
  grid <- nonlinear_heat$grid
  stopifnot(heat_from >= min(grid), heat_to <= max(grid), heat_to > heat_from)
  nodes <- findInterval(c(heat_from, heat_to),
    (head(grid, -1L) + tail(grid, -1L)) / 2) + 1L
  weights <- if ("acute_heat" %in% design_names)
    list(acute_heat = heat_to - heat_from) else list()
  if (cluster != reference_cluster) {
    term <- intersect(c(paste0("cluster", cluster, ":acute_heat"),
      paste0("acute_heat:cluster", cluster)), design_names)
    if (length(term) > 1L) stop("Ambiguous interaction for cluster ", cluster)
    if (length(term)) weights[[term]] <- heat_to - heat_from
  }
  curve_weights <- numeric(length(grid))
  curve_weights[nodes[1]] <- curve_weights[nodes[1]] - 1
  curve_weights[nodes[2]] <- curve_weights[nodes[2]] + 1
  weights[[nonlinear_heat$shared_field]] <- curve_weights
  if (!is.na(climate_zone) && nzchar(climate_zone)) {
    zone <- match(climate_zone, nonlinear_heat$zone_levels)
    stopifnot(!is.na(zone))
    deviation_weights <- numeric(length(grid) * length(nonlinear_heat$zone_levels))
    block <- (zone - 1L) * length(grid) + seq_along(grid)
    deviation_weights[block] <- curve_weights
    weights[[nonlinear_heat$deviation_field]] <- deviation_weights
  }
  weights
}

make_full_heat_lincombs <- function(data, nonlinear_heat, model_formula = heat_model_formula()) {
  stopifnot(isTRUE(nonlinear_heat$pooled), length(nonlinear_heat$zone) == nrow(data))
  combinations <- unique(data.frame(cluster = as.character(data$cluster),
    climate_zone = nonlinear_heat$zone, stringsAsFactors = FALSE))
  combinations$lincomb_id <- paste0("heat_1SD_cluster", combinations$cluster,
    "_cz", ifelse(is.na(combinations$climate_zone), "missing", combinations$climate_zone))
  stopifnot(!anyDuplicated(combinations$lincomb_id))
  fixed_formula <- heat_fixed_formula(model_formula)
  design_names <- colnames(stats::model.matrix(fixed_formula, data = data))
  lincombs <- lapply(seq_len(nrow(combinations)), function(j) {
    weights <- full_heat_contrast_weights(combinations$cluster[j],
      levels(data$cluster)[1], combinations$climate_zone[j], nonlinear_heat, design_names)
    do.call(INLA::inla.make.lincomb, weights)[[1]]
  })
  names(lincombs) <- combinations$lincomb_id
  list(lincombs = lincombs, combinations = combinations,
    heat_from = 0, heat_to = 1, method = "INLA linear-combination marginals")
}

heat_object_md5 <- function(object) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(object, path, compress = FALSE)
  unname(tools::md5sum(path))
}

fit_heat_model <- function(prepared, cache_dir, family = c("poisson", "nbinomial"),
                           force_refit = FALSE, model_formula = heat_model_formula(), cache_tag = "",
                           nonlinear_heat = NULL) {
  family <- match.arg(family)
  data <- prepared$daily_ed
  heat_lincombs <- make_heat_lincombs(data, model_formula)
  full_heat <- if (!is.null(nonlinear_heat)) make_full_heat_lincombs(data, nonlinear_heat, model_formula) else NULL
  if (!is.null(full_heat)) heat_lincombs <- c(heat_lincombs, full_heat$lincombs)
  compute_settings <- list(dic = TRUE, waic = TRUE, config = TRUE)
  predictor_settings <- list(compute = TRUE)
  signature <- heat_object_md5(list(
    data = data, formula = deparse(model_formula), family = family,
    lincombs = heat_lincombs, contrasts = getOption("contrasts"),
    control.compute = compute_settings, control.predictor = predictor_settings,
    INLA_version = as.character(utils::packageVersion("INLA")), cache_version = 1L))
  cache_file <- file.path(cache_dir, paste0("heat_", family, cache_tag, "_fit.rds"))
  cached <- if (!force_refit && file.exists(cache_file)) readRDS(cache_file) else NULL
  cluster_terms <- names(heat_lincombs)
  if (!is.null(cached) && identical(cached$signature, signature) &&
      all(cluster_terms %in% names(cached$fit$marginals.lincomb.derived))) {
    message("Using shared cached ", family, " heat model.")
    return(cached$fit)
  }
  message("Fitting ", family, " heat model with heat-contrast linear combinations...")
  fit <- INLA::inla(
    # E supplies the population offset: expected count = Population * exp(eta).
    model_formula, family = family, data = data, E = data$Population,
    lincomb = heat_lincombs, control.compute = compute_settings,
    control.predictor = predictor_settings)
  if (!all(cluster_terms %in% names(fit$marginals.lincomb.derived))) {
    stop("INLA did not return all requested heat-contrast marginals.")
  }
  if (!is.null(full_heat)) {
    fit$full_heat_irr <- full_heat[setdiff(names(full_heat), "lincombs")]
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
