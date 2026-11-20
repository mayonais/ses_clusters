# NB climate-zone model with ZIP iid and PC hyperpriors. Acute heat enters only
# through the shared and correlated zone RW2 curves, without linear SES slopes.
.libPaths(c("C:/Users/yysma/AppData/Local/R/win-library/4.4", .libPaths()))

# Two-sided PC prior on the five-zone exchangeable correlation, base rho = 0.
# Each side has prior probability 1/2. Distance is sqrt(-log(det(R(rho))));
# R(rho) has diagonal 1 and common off-diagonal rho. This preserves the full
# support (-1/(K-1), 1), unlike substituting an AR(1) correlation prior.
exchangeable_pc_prior <- function(K, u = 0.5, alpha = 0.1) {
  stopifnot(K >= 2, u > 0, u < 1, alpha > 0, alpha < 0.5)
  distance_u <- sqrt(-(K - 1) * log1p(-u) - log1p((K - 1) * u))
  lambda <- -log(2 * alpha) / distance_u
  log_density <- function(theta) {
    a <- K - 1
    log_sum <- pmax(theta, log(a)) + log1p(exp(-abs(theta - log(a))))
    D <- K * log_sum - theta - K * log(K)
    deriv <- K * plogis(theta - log(a)) - 1
    near <- abs(theta) < 1e-3
    t <- theta[near]
    D[near] <- a / (2 * K) * t^2 + a * (K - 2) / (6 * K^2) * t^3 +
      a * (K^2 - 6 * K + 6) / (24 * K^3) * t^4
    deriv[near] <- a / K * t + a * (K - 2) / (2 * K^2) * t^2 +
      a * (K^2 - 6 * K + 6) / (6 * K^3) * t^3
    d <- sqrt(pmax(0, D))
    jac <- abs(deriv) / (2 * d)
    jac[theta == 0] <- sqrt(a / (2 * K))
    log(lambda / 2) - lambda * d + log(jac)
  }
  # INLA interpolates a log-density table on its internal theta scale.
  # theta = log((1+(K-1)*rho)/(1-rho)); the Jacobian is already included.
  bounds <- c(-300, 120)
  outer_step <- 0.05
  theta <- sort(unique(c(seq(bounds[1], -10, by = outer_step), seq(-10, 10, by = 0.002),
    seq(10, bounds[2], by = outer_step), 0)))
  logp <- log_density(theta)
  stopifnot(all(is.finite(logp)))
  mass <- integrate(function(x) exp(log_density(x)), bounds[1], 0, rel.tol = 1e-9)$value +
    integrate(function(x) exp(log_density(x)), 0, bounds[2], rel.tol = 1e-9)$value
  theta_u <- log((1 + (K - 1) * u) / (1 - u))
  tail <- integrate(function(x) exp(log_density(x)), theta_u, bounds[2], rel.tol = 1e-9)$value
  interpolated <- splinefun(theta, logp)
  interpolation_mass <- integrate(function(x) exp(interpolated(x)), bounds[1], 0,
    subdivisions = 1000L, rel.tol = 1e-7)$value +
    integrate(function(x) exp(interpolated(x)), 0, bounds[2],
      subdivisions = 1000L, rel.tol = 1e-7)$value
  stopifnot(abs(mass - 1) < 1e-7, abs(tail - alpha) < 1e-7,
    abs(interpolation_mass - 1) < 1e-6)
  list(hyper = list(rho = list(prior = paste(c("table:", sprintf("%.17g", theta),
    sprintf("%.17g", logp)), collapse = " "), initial = 0)),
    table = data.frame(theta = theta, log_density = logp),
    specification = list(type = "PC exchangeable correlation", K = K, base = 0,
      u = u, alpha = alpha, calibration = sprintf("P(rho > %g) = %g", u, alpha),
      branch_probability = 0.5, lambda = lambda, mass = mass,
      interpolated_mass = interpolation_mass, calibrated_tail = tail,
      support = c(-1 / (K - 1), 1), internal_table_bounds = range(theta)))
}

load_shared_heat_inputs <- function(input_dir = "Los_Angeles_discrete_none") {
  source("shared_functions.R", local = TRUE)
  files <- c(layout = file.path(input_dir, "heat_model_layout.rds"),
    grid = file.path(input_dir, "heat_grid.rds"))
  if (!all(file.exists(files))) stop("Missing shared heat input: ",
    paste(files[!file.exists(files)], collapse = ", "))
  layout <- readRDS(files[["layout"]])
  grid <- readRDS(files[["grid"]])
  stopifnot(identical(layout$schema_version, 1L), identical(grid$schema_version, 1L),
    !"D2" %in% names(layout$data), "D2" %in% layout$data_columns,
    identical(setdiff(layout$data_columns, "D2"), names(layout$data)),
    identical(layout$row_layout_md5, grid$row_layout_md5),
    identical(heat_object_md5(layout$data[c("zip", "date")]), layout$row_layout_md5),
    nrow(layout$data) == length(layout$exposure), length(layout$fixed_names) == 26L,
    !any(grepl("acute_heat", layout$fixed_names)),
    isTRUE(grid$nonlinear_heat$pooled), length(grid$nonlinear_heat$grid) == 50L,
    length(grid$nonlinear_heat$zone) == nrow(layout$data),
    length(grid$nonlinear_heat$observed_index) == nrow(layout$data))
  list(layout = layout, nonlinear_heat = grid$nonlinear_heat, files = files,
    md5 = setNames(unname(tools::md5sum(files)), names(files)))
}

preserved_nb_inputs <- function(precision_prior = list(prior = "pc.prec", param = c(1, 0.01))) {
  stopifnot(precision_prior$prior %in% c("pc.prec", "loggamma"),
    length(precision_prior$param) == 2L, all(is.finite(precision_prior$param)),
    all(precision_prior$param > 0))
  library(INLA)
  library(dplyr)
  source("shared_functions.R")
  source("cluster_ar.R")
  shared <- load_shared_heat_inputs()
  layout <- shared$layout
  data <- layout$data
  meta <- shared$nonlinear_heat
  fixed_formula <- layout$fixed_formula
  # Keep rows, contrasts, offsets and the heat grid, but use reconstructed counts
  # from suppressed_data.R, including the one-visit all-zero-month override.
  prepared <- load_heat_data("Los_Angeles_discrete_none")
  for (field in setdiff(names(prepared$daily_ed), "D2"))
    stopifnot(isTRUE(all.equal(data[[field]], prepared$daily_ed[[field]], check.attributes = FALSE)))
  data$D2 <- prepared$daily_ed$D2
  data <- data[layout$data_columns]
  data$zip_id <- match(as.character(data$zip), sort(unique(as.character(data$zip))))
  data$climate_zone_group <- match(meta$zone, meta$zone_levels)
  data$climate_zone_group[is.na(meta$zone)] <- 1L
  reference <- cluster_heat_reference(data)
  stopifnot(all(reference >= min(meta$grid)), all(reference <= max(meta$grid)))
  meta$counterfactual_index <- findInterval(reference,
    (head(meta$grid, -1L) + tail(meta$grid, -1L)) / 2) + 1L
  data$acute_heat_rw2_cf <- meta$counterfactual_index
  meta$ar_reference <- "SES-cluster-month mean heat"
  meta$ses_heat_slopes <- meta$common_heat_slope <- FALSE
  meta$group_model <- "exchangeable"
  meta$model_variant <- "climate_zone_added_nb_pc_no_linear_heat_slopes"
  correlation <- exchangeable_pc_prior(length(meta$zone_levels))
  formula_env <- new.env(parent = .GlobalEnv)
  formula_env$precision_hyper <- list(prec = precision_prior)
  formula_env$zone_correlation_hyper <- correlation$hyper
  rhs <- fixed_formula[[3L]]
  terms <- list(
    quote(f(doy, model = "rw2", scale.model = TRUE, constr = TRUE, hyper = precision_hyper)),
    quote(f(acute_heat_rw2_shared, model = "rw2", n = 50L, values = seq_len(50L),
      scale.model = TRUE, constr = TRUE, hyper = precision_hyper)),
    quote(f(acute_heat_rw2_dev, model = "rw2", n = 50L, values = seq_len(50L),
      scale.model = TRUE, constr = TRUE, hyper = precision_hyper,
      group = climate_zone_group, ngroup = 5L,
      control.group = list(model = "exchangeable", hyper = zone_correlation_hyper))),
    quote(f(zip_id, model = "iid", hyper = precision_hyper)))
  for (term in terms) rhs <- call("+", rhs, term)
  formula <- as.formula(call("~", as.name("D2"), rhs), env = formula_env)
  new_design <- model.matrix(heat_fixed_formula(formula), data)
  stopifnot(identical(colnames(new_design), layout$fixed_names),
    identical(heat_object_md5(as.numeric(new_design)), layout$fixed_design_numeric_md5),
    identical(attr(new_design, "contrasts"), layout$factor_contrasts),
    ncol(new_design) == 26L,
    identical(layout$exposure, as.numeric(data$Population)),
    length(unique(data$zip_id)) == 309L, nrow(data) == 54910L)
  result <- list(data = data, nonlinear_heat = meta, formula = formula,
    fixed_prior = list(mean = 0, prec = 0.001, mean.intercept = 0, prec.intercept = 0),
    family_prior = list(hyper = list(theta = list(prior = "pc.mgamma", param = 7))),
    correlation = correlation, precision_prior = precision_prior, source_files = shared$files,
    source_md5 = shared$md5, reference_cluster = prepared$reference_cluster,
    heat_sd_f = prepared$heat_sd, fixed_design_md5 = heat_object_md5(new_design))
  message("Verified 26 fixed coefficients: all 14 linear acute-heat slopes removed; other fixed effects retained.")
  result
}

restore_model_formula <- function(inputs) {
  formula <- inputs[["formula"]]
  if (!is.null(formula)) return(formula)
  formula_env <- new.env(parent = .GlobalEnv)
  formula_env$precision_hyper <- inputs$precision_hyper
  formula_env$zone_correlation_hyper <- inputs$zone_correlation_hyper
  as.formula(paste(inputs$formula_text, collapse = "\n"), env = formula_env)
}

fit_inla_model <- function(inputs, output_dir, family = "nbinomial",
  fixed_prior = inputs$fixed_prior, lincomb = NULL, save_config = FALSE, start = NULL) {
  formula <- restore_model_formula(inputs)
  exposure <- if (is.null(inputs$exposure)) inputs$data$Population else inputs$exposure
  family_prior <- if (family == "nbinomial") inputs$family_prior else list()
  if (is.null(start)) start <- list()
  gc(verbose = FALSE)
  elapsed <- system.time(fit <- INLA::inla(formula, family = family,
    data = inputs$data, E = exposure, control.fixed = fixed_prior,
    control.family = family_prior, lincomb = lincomb, control.mode = start,
    control.compute = list(dic = TRUE, waic = TRUE, config = save_config),
    control.predictor = list(compute = TRUE), num.threads = "2:1", verbose = TRUE,
    working.directory = file.path(normalizePath(output_dir, winslash = "/", mustWork = TRUE), "inla_work")))[["elapsed"]]
  stopifnot(fit$mode$mode.status == 0L, all(is.finite(fit$summary.fixed$mean)),
    all(is.finite(fit$summary.hyperpar$mean)))
  fit$analysis_elapsed_seconds <- elapsed
  fit
}

summarise_heat_fit <- function(fit, full) {
  heat <- full$combinations
  values <- lapply(heat$lincomb_id, function(id) {
    marginal <- fit$marginals.lincomb.derived[[id]]
    stopifnot(!is.null(marginal))
    q <- INLA::inla.qmarginal(c(0.025, 0.5, 0.975), marginal)
    c(IRR_median = exp(q[2]), IRR_lower = exp(q[1]), IRR_upper = exp(q[3]),
      P_IRR_gt_1 = 1 - INLA::inla.pmarginal(0, marginal))
  })
  heat <- cbind(heat, as.data.frame(do.call(rbind, values)))
  stopifnot(all(is.finite(as.matrix(heat[c("IRR_median", "IRR_lower", "IRR_upper", "P_IRR_gt_1")]))),
    all(heat$IRR_lower > 0), all(heat$P_IRR_gt_1 >= 0 & heat$P_IRR_gt_1 <= 1))
  list(heat = heat,
    fixed = data.frame(term = rownames(fit$summary.fixed), fit$summary.fixed, row.names = NULL),
    hyper = data.frame(term = rownames(fit$summary.hyperpar), fit$summary.hyperpar, row.names = NULL),
    mode_status = fit$mode$mode.status)
}

write_fit_details <- function(fit, output_dir, results = NULL,
  summary_file = "diagnostics.txt", prior_file = "fitted_priors.txt") {
  capture.output(summary(fit), file = file.path(output_dir, summary_file))
  capture.output(str(fit$all.hyper, max.level = 6), file = file.path(output_dir, prior_file))
  if (!is.null(results)) {
    saveRDS(results, file.path(output_dir, "summaries.rds"))
    for (name in c("heat", "fixed", "hyper"))
      write.csv(results[[name]], file.path(output_dir, paste0(name, "_details.csv")), row.names = FALSE)
  }
}

protected_output_hashes <- function(output_dir, extra = character()) {
  paths <- c(list.files("Los_Angeles_discrete_none/climate_zone_added_nb_pc", full.names = TRUE),
    list.files("maps/nb_pc_prior", full.names = TRUE, recursive = TRUE),
    list.files("methods", full.names = TRUE, recursive = TRUE), extra)
  paths <- sort(unique(normalizePath(paths[!dir.exists(paths)], winslash = "/", mustWork = TRUE)))
  hashes <- tools::md5sum(paths)
  file <- file.path(output_dir, "protected_files.rds")
  if (file.exists(file)) stopifnot(identical(hashes, readRDS(file))) else saveRDS(hashes, file)
  hashes
}

fit_preserved_nb <- function(output_dir = "Los_Angeles_discrete_none/climate_zone_added_nb_pc",
  fit_only = FALSE, preflight = FALSE, save_config = TRUE, restart_fit = FALSE,
  precision_prior = list(prior = "pc.prec", param = c(1, 0.01))) {
  inputs <- preserved_nb_inputs(precision_prior = precision_prior)
  full <- make_full_heat_lincombs(inputs$data, inputs$nonlinear_heat, inputs$formula)
  stopifnot(nrow(full$combinations) == 42L)
  signature_inputs <- list(data = inputs$data, formula = deparse(inputs$formula),
    precision = if (identical(precision_prior$prior, "pc.prec")) precision_prior$param else precision_prior,
    correlation = inputs$correlation$specification,
    fixed = inputs$fixed_prior, family = inputs$family_prior, source = inputs$source_md5,
    INLA = as.character(packageVersion("INLA")), version = 3L)
  signature <- heat_object_md5(signature_inputs)
  if (preflight) {
    print(inputs$formula)
    print(inputs$correlation$specification)
    message("Preflight passed: 54910 rows, 309 ZIPs, no linear acute-heat slopes, 42 curve-only contrasts.")
    return(invisible(inputs))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  signature_file <- file.path(output_dir, "signature.txt")
  model_file <- file.path(output_dir, "nbinomial_model.rds")
  if (file.exists(signature_file) && file.exists(model_file))
    stopifnot(identical(readLines(signature_file), signature))
  writeLines(signature, signature_file)
  initial_mode <- list()
  if (restart_fit) {
    previous <- readRDS(model_file)
    stopifnot(previous$mode$mode.status == 0L, identical(previous$.args$data, inputs$data))
    initial_mode <- list(theta = previous$mode$theta, x = previous$mode$x,
      restart = TRUE, fixed = FALSE)
    rm(previous); gc(verbose = FALSE)
    message("Restarting from the saved hyperparameter and latent modes; priors remain unchanged.")
  }
  if (!file.exists(model_file) || restart_fit) {
    lincombs <- c(make_heat_lincombs(inputs$data, inputs$formula), full$lincombs)
    message("Fitting preserved climate_zone_added NB model; precision prior: ",
      precision_prior$prior, " (", paste(precision_prior$param, collapse = ", "), ").")
    fit <- fit_inla_model(inputs, output_dir, lincomb = lincombs,
      save_config = save_config, start = initial_mode)
    fit$nonlinear_heat <- inputs$nonlinear_heat
    fit$full_heat_irr <- full
    fit$preserved_nb_specification <- list(source_files = inputs$source_files,
      source_md5 = inputs$source_md5, fixed_design_md5 = inputs$fixed_design_md5,
      correlation = inputs$correlation$specification, precision_prior = precision_prior,
      elapsed = fit$analysis_elapsed_seconds,
      restarted_from_saved_mode = restart_fit)
    stopifnot(fit$mode$mode.status == 0L, all(is.finite(fit$summary.fixed$mean)),
      all(is.finite(fit$summary.hyperpar$mean)),
      identical(fit$names.fixed, colnames(model.matrix(heat_fixed_formula(inputs$formula), inputs$data))))
    saveRDS(fit, model_file)
  } else fit <- readRDS(model_file)
  diagnostics <- summarise_zip_fit(fit, inputs$data, "nbinomial", inputs$reference_cluster)
  results <- c(list(metrics = diagnostics$metrics, clusters = diagnostics$clusters),
    summarise_heat_fit(fit, full), list(specification = fit$preserved_nb_specification,
    sampling_info = list(n_observations = nrow(inputs$data), names.fixed = fit$names.fixed)))
  write_fit_details(fit, output_dir, results)
  write.csv(inputs$correlation$table, file.path(output_dir, "correlation_pc_prior.csv"), row.names = FALSE)
  precision_description <- if (identical(precision_prior$prior, "pc.prec"))
    sprintf("Four precision priors: pc.prec, P(1/sqrt(tau) > %g) = %g.",
      precision_prior$param[1], precision_prior$param[2]) else
    sprintf("Four precision priors: Gamma(shape=%g, rate=%g), INLA loggamma on log precision.",
      precision_prior$param[1], precision_prior$param[2])
  writeLines(c("PRESERVED CLIMATE_ZONE_ADDED: NEGATIVE BINOMIAL PRIOR SPECIFICATION",
    "26 fixed coefficients: common acute heat and all 13 SES heat interactions removed.",
    "Original rows, population offset, 50-node heat grid, scaling and RW2 sum-to-zero constraints retained; reconstructed counts updated.",
    "Added iid ZIP random intercept; zone replication replaced by exchangeable group correlation.",
    precision_description,
    "Zone correlation: two-sided exchangeable PC prior, base rho=0; P(rho>0.5)=0.1; each sign has mass 0.5.",
    "Correlation support: (-0.25, 1); implemented as an internal-scale log-density table including its Jacobian.",
    "NB size: pc.mgamma, parameter 7; non-intercept fixed coefficients: Normal(0, variance=1000); intercept: flat.",
    if (fit_only) "Fit-only summaries: no excess-rate posterior draws generated." else
      "Excess-rate reference: SES-cluster-month mean heat; equal ZIP-day weighting; 2000 posterior draws.",
    paste("Shared input:", inputs$source_files),
    paste("Shared input MD5:", names(inputs$source_md5), inputs$source_md5),
    paste("Reduced fixed-design MD5:", inputs$fixed_design_md5),
    paste("INLA version:", packageVersion("INLA")), paste("R version:", R.version.string)),
    file.path(output_dir, "priors.txt"))
  rm(fit, inputs); gc(verbose = FALSE)
  if (fit_only) return(invisible(results))
  export_preserved_nb_outputs(results_dir = output_dir)
  invisible(results)
}

validate_nb_pc_model <- function(
  results_dir = "Los_Angeles_discrete_none/climate_zone_added_nb_pc",
  with_maps = TRUE, map_dir = "maps/nb_pc_prior",
  precision_prior = list(prior = "pc.prec", param = c(1, 0.01)),
  with_excess_rates = with_maps, check_sampling = TRUE) {
  suppressPackageStartupMessages(library(INLA))
  suppressPackageStartupMessages(library(dplyr))
  source("shared_functions.R", local = TRUE)
  fit <- readRDS(file.path(results_dir, "nbinomial_model.rds"))
  results <- readRDS(file.path(results_dir, "summaries.rds"))
  checks <- list()
  check <- function(name, value, detail = NULL) {
    stopifnot(isTRUE(value))
    checks[[length(checks) + 1L]] <<- list(name = name, passed = TRUE, detail = detail)
    message("Passed: ", name)
  }

  # Independent two-zone benchmark: the generic exchangeable PC prior reduces to
  # INLA's native PC correlation prior after the correct theta Jacobian.
  pc2 <- exchangeable_pc_prior(2)
  theta2 <- c(-6, -2, -0.5, 0.5, 2, 6)
  rho2 <- tanh(theta2 / 2)
  native <- INLA::inla.pc.dcor0(rho2, u = 0.5, alpha = 0.2, log = TRUE) +
    log(0.5 * (1 - rho2^2))
  tabulated <- splinefun(pc2$table$theta, pc2$table$log_density)(theta2)
  check("Correlation PC prior matches native INLA two-zone benchmark",
    max(abs(native - tabulated)) < 1e-6, max(abs(native - tabulated)))
  pc5 <- exchangeable_pc_prior(5)
  theta5 <- c(-3, -0.5, 0.4, 1.6, 5)
  rho5 <- (exp(theta5) - 1) / (exp(theta5) + 4)
  distance5 <- vapply(rho5, function(rho) {
    R <- matrix(rho, 5, 5); diag(R) <- 1
    stopifnot(min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) > 0)
    sqrt(-as.numeric(determinant(R, logarithm = TRUE)$modulus))
  }, numeric(1))
  step <- 1e-5
  direct_distance <- function(t) {
    rho <- (exp(t) - 1) / (exp(t) + 4)
    R <- matrix(rho, 5, 5); diag(R) <- 1
    sqrt(-as.numeric(determinant(R, logarithm = TRUE)$modulus))
  }
  jacobian5 <- vapply(theta5, function(t)
    abs((direct_distance(t + step) - direct_distance(t - step)) / (2 * step)), numeric(1))
  direct_logp <- log(pc5$specification$lambda / 2) - pc5$specification$lambda * distance5 + log(jacobian5)
  table_logp <- splinefun(pc5$table$theta, pc5$table$log_density)(theta5)
  check("Five-zone correlation density agrees with matrix determinant and numerical Jacobian",
    max(abs(direct_logp - table_logp)) < 1e-5, max(abs(direct_logp - table_logp)))

  shared <- load_shared_heat_inputs()
  layout <- shared$layout
  check("Shared input provenance matches the fitted model",
    identical(fit$preserved_nb_specification$source_files, shared$files) &&
      identical(fit$preserved_nb_specification$source_md5, shared$md5))
  new_design <- model.matrix(heat_fixed_formula(fit$.args$formula), fit$.args$data)
  has_linear_heat <- "acute_heat" %in% fit$names.fixed
  check("All 14 linear heat slopes removed; shared fixed design preserved",
    !has_linear_heat && identical(colnames(new_design), layout$fixed_names) &&
      identical(heat_object_md5(as.numeric(new_design)), layout$fixed_design_numeric_md5) &&
      identical(attr(new_design, "contrasts"), layout$factor_contrasts) &&
      identical(fit$names.fixed, layout$fixed_names),
    list(coefficients = ncol(new_design), observations = nrow(new_design)))
  check("Shared ZIP/date row order preserved",
    identical(fit$.args$data[c("zip", "date")], layout$data[c("zip", "date")]))
  if (!has_linear_heat) {
    prepared <- load_heat_data("Los_Angeles_discrete_none")
    check("Fitted counts match the current suppressed-count reconstruction",
      identical(as.numeric(fit$.args$data$D2), as.numeric(prepared$daily_ed$D2)))
    check("Linear heat metadata disabled",
      identical(fit$nonlinear_heat$ses_heat_slopes, FALSE) &&
        identical(fit$nonlinear_heat$common_heat_slope, FALSE))
    rm(prepared)
  }
  check("Shared population offsets and heat grid preserved",
    identical(as.numeric(fit$.args$E), layout$exposure) &&
      identical(fit$nonlinear_heat$grid, shared$nonlinear_heat$grid))
  check("Negative-binomial likelihood and converged mode",
    identical(fit$.args$family, "nbinomial") && fit$mode$mode.status == 0L)
  rm(shared, layout, new_design); gc(verbose = FALSE)

  # Check the Gaussian backend used for excess-rate sampling against the saved
  # INLA conditional variances at three distinct posterior configurations.
  if (check_sampling) {
  configurations <- fit$misc$configs
  covariance_errors <- vapply(unique(c(1L, ceiling(configurations$nconfig / 2),
    configurations$nconfig)), function(i) {
    config <- configurations$config[[i]]
    gaussian <- heat_matrix_qsample(n = 1L, Q = config$Q, mu = config$mean,
      constr = configurations$constr, seed = 20260910L, logdens = TRUE,
      compute.mean = TRUE, .debug = TRUE)
    saved_variance <- as.numeric(Matrix::diag(config$Qinv))
    max(abs(gaussian$variance - saved_variance) / pmax(abs(saved_variance), 1e-8))
  }, numeric(1))
  check("Gaussian sampling backend matches saved constrained INLA variances",
    max(covariance_errors) < 1e-5,
    list(configurations = 3L, max_relative_variance_error = max(covariance_errors)))
  }

  random <- fit$all.hyper$random
  check("Exactly four random-effect precision terms",
    setequal(vapply(random, function(x) x$hyperid, character(1)),
      c("doy", "acute_heat_rw2_shared", "acute_heat_rw2_dev", "zip_id")))
  for (term in random) {
    precision <- term$hyper[[1L]]
    check(paste("Expected precision prior", term$hyperid),
      identical(as.character(precision$prior), precision_prior$prior) &&
        isTRUE(all.equal(as.numeric(precision$param), precision_prior$param)))
    if (term$hyperid == "acute_heat_rw2_dev") {
      correlation <- term$group.hyper[[1L]]
      check("Exchangeable zone correlation uses the calibrated PC table",
        identical(as.character(correlation$prior), pc5$hyper$rho$prior) &&
          identical(fit$nonlinear_heat$group_model, "exchangeable") &&
          isTRUE(all.equal(fit$preserved_nb_specification$correlation$calibrated_tail, 0.1)) &&
          identical(fit$preserved_nb_specification$correlation$u, 0.5))
    }
  }
  dispersion <- fit$all.hyper$family[[1L]]$hyper[[1L]]
  check("PC prior on negative-binomial size",
    identical(as.character(dispersion$prior), "pc.mgamma") &&
      identical(as.numeric(dispersion$param), 7))
  check("42 finite cluster-by-zone full heat contrasts",
    nrow(results$heat) == 42L && all(is.finite(as.matrix(results$heat[
      c("IRR_median", "IRR_lower", "IRR_upper", "P_IRR_gt_1")]))) &&
      all(results$heat$IRR_lower > 0) && all(results$heat$P_IRR_gt_1 >= 0 & results$heat$P_IRR_gt_1 <= 1))

  # Independently check every saved linear-combination coefficient. Historical
  # fits include fixed slopes; the current fit uses only shared/zone curve nodes.
  grid <- fit$nonlinear_heat$grid
  nodes <- findInterval(c(0, 1), (head(grid, -1L) + tail(grid, -1L)) / 2) + 1L
  weights_match <- vapply(seq_len(nrow(results$heat)), function(i) {
    row <- results$heat[i, ]
    expected <- c(if (has_linear_heat) c("acute_heat:1" = 1),
      setNames(c(-1, 1), paste0("acute_heat_rw2_shared:", nodes)))
    if (paste0("cluster", row$cluster, ":acute_heat") %in% fit$names.fixed)
      expected <- c(expected, setNames(1, paste0("cluster", row$cluster, ":acute_heat:1")))
    if (!is.na(row$climate_zone)) {
      zone <- match(row$climate_zone, fit$nonlinear_heat$zone_levels)
      expected <- c(expected, setNames(c(-1, 1),
        paste0("acute_heat_rw2_dev:", (zone - 1L) * length(grid) + nodes)))
    }
    actual <- unlist(lapply(fit$.args$lincomb[[row$lincomb_id]], function(component) {
      field <- names(component)[1L]
      part <- component[[1L]]
      index <- if (is.null(part$idx)) 1L else part$idx
      setNames(as.numeric(part$weight), paste0(field, ":", index))
    }), use.names = TRUE)
    identical(expected[order(names(expected))], actual[order(names(actual))])
  }, logical(1))
  check("All 42 saved contrasts match the fitted fixed terms and shared/zone curves", all(weights_match))

  # Separately approximated INLA marginal means need not commute exactly with
  # linear combinations. Record the discrepancy as a numerical diagnostic rather
  # than asserting equality between these two different marginal approximations.
  fixed_mean <- setNames(fit$summary.fixed$mean, rownames(fit$summary.fixed))
  shared_mean <- fit$summary.random$acute_heat_rw2_shared$mean
  dev_mean <- fit$summary.random$acute_heat_rw2_dev$mean
  reconstructed <- vapply(seq_len(nrow(results$heat)), function(i) {
    row <- results$heat[i, ]
    slope <- if (has_linear_heat) fixed_mean[["acute_heat"]] else 0
    interaction <- paste0("cluster", row$cluster, ":acute_heat")
    if (interaction %in% names(fixed_mean)) slope <- slope + fixed_mean[[interaction]]
    value <- slope + diff(shared_mean[nodes])
    if (!is.na(row$climate_zone)) {
      zone <- match(row$climate_zone, fit$nonlinear_heat$zone_levels)
      value <- value + diff(dev_mean[(zone - 1L) * length(grid) + nodes])
    }
    value
  }, numeric(1))
  derived_mean <- fit$summary.lincomb.derived[results$heat$lincomb_id, "mean"]
  contrast_error <- max(abs(reconstructed - derived_mean))
  message("Diagnostic: largest difference between separately approximated means = ",
    signif(contrast_error, 6), " log-IRR units.")

  if (with_excess_rates) {
    ar <- readRDS(file.path(results_dir, "cluster_excess_rates.rds"))
    check("2000 finite signed joint excess-rate draws for each SES cluster",
      identical(dim(ar$draws), c(14L, 2000L)) && all(is.finite(ar$draws)) &&
        all(ar$summary$n_samples == 2000L) &&
        identical(as.numeric(rowMeans(ar$draws > 0)), as.numeric(ar$summary$P_AR_gt_0)))
    check("Excess-rate provenance matches the corrected model",
      identical(ar$provenance$model_md5, unname(tools::md5sum(file.path(results_dir, "nbinomial_model.rds")))) &&
        identical(ar$settings$counterfactual, "SES-cluster-month mean heat"))
  }
  if (with_maps) {
    map_names <- c("acute_heat_irr_1sd_map", "acute_heat_probability_map", "excess_rate_probability_map")
    files <- as.vector(outer(file.path(map_dir, map_names), c("png", "pdf"), paste, sep = "."))
    check("Three requested maps exported as PNG and PDF", all(file.exists(files)) && all(file.info(files)$size > 0))
    map_files <- data.frame(path = files, md5 = unname(tools::md5sum(files)), size = file.info(files)$size)
  } else map_files <- NULL
  validation_path <- file.path(results_dir, "validation.json")
  plot_settings <- if (file.exists(validation_path))
    jsonlite::read_json(validation_path)$plot_settings else NULL
  report <- list(status = "completed", checks = checks, map_files = map_files,
    sampling_backend_checked = check_sampling,
    plot_settings = plot_settings,
    separate_marginal_mean_discrepancy_log_irr = contrast_error,
    model_md5 = unname(tools::md5sum(file.path(results_dir, "nbinomial_model.rds"))),
    correlation_pc_prior = pc5$specification, precision_prior = precision_prior,
    interpretation = "Finite implementation checks passed; approximate Bayesian inference and the specified prior calibration remain modeling assumptions.")
  jsonlite::write_json(report, file.path(results_dir, "validation.json"), auto_unbox = TRUE, pretty = TRUE, digits = 16)
  writeLines(c("PRESERVED CLIMATE_ZONE_ADDED NB VALIDATION", vapply(checks, function(x) paste("PASSED:",x$name), character(1))),
    file.path(results_dir, "validation.txt"))
  message("Validation completed.")

  invisible(report)
}

export_preserved_nb_outputs <- function(
  results_dir = "Los_Angeles_discrete_none/climate_zone_added_nb_pc",
  map_dir = "maps/nb_pc_prior") {
  source("shared_functions.R", local = TRUE)
  source("cluster_ar.R", local = TRUE)
  model_file <- file.path(results_dir, "nbinomial_model.rds")
  stopifnot(file.exists(model_file), file.exists(file.path(results_dir, "summaries.rds")))
  run_saved_heat_excess_rates(source_model = model_file, output_dir = results_dir)
  map_environment <- new.env(parent = environment())
  map_environment$.map_functions_only <- TRUE
  source("acute_heat_irr_maps.R", local = map_environment)
  map_environment$export_nb_pc_figures(output_dir = map_dir, results_dir = results_dir)
  validate_nb_pc_model(results_dir = results_dir, with_maps = TRUE, map_dir = map_dir)
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


# Fit the same current predictor with Poisson counts, keeping the NB fit read-only.
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
    family_prior = list(hyper = list(theta = list(prior = "pc.mgamma", param = 7))),
    nonlinear_heat = fit$nonlinear_heat, names_fixed = fit$names.fixed,
    baseline_md5 = unname(tools::md5sum(baseline_file)), INLA_version = as.character(packageVersion("INLA")),
    size_prior = size)
  list(inputs = inputs, reference = list(zero = likelihood_zero_diagnostics(fit, inputs$data, "nbinomial"),
    heat = saved$heat, fixed = saved$fixed, hyper = saved$hyper, mode_status = fit$mode$mode.status))
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

run_coefficient_sensitivity <- function(
  output_dir = "prior_sensitivity_nb/no_linear_heat_slopes/coefficient_variance", restart_100 = FALSE) {
  suppressPackageStartupMessages(library(INLA))
  source("shared_functions.R")
  baseline_file <- "Los_Angeles_discrete_none/climate_zone_added_nb_pc/nbinomial_model.rds"
  comparison_dir <- output_dir
  input_file <- file.path(comparison_dir, "matched_inputs.rds")
  reference_file <- file.path(comparison_dir, "baseline_variance_1000.rds")
  extra <- character()
  if (restart_100) {
    original_dir <- file.path(comparison_dir, "variance_100")
    original_file <- file.path(original_dir, "nbinomial_model.rds")
    if (!all(file.exists(c(input_file, reference_file, original_file))))
      stop("Run --coefficient-sensitivity before restarting its variance-100 fit.")
    output_dir <- file.path(comparison_dir, "variance_100_rerun")
    extra <- c(list.files(comparison_dir, full.names = TRUE), list.files(original_dir, full.names = TRUE),
      list.files(file.path(comparison_dir, "variance_10"), full.names = TRUE))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  before <- protected_output_hashes(output_dir, extra)
  protected <- names(before)
  if (!file.exists(input_file) || !file.exists(reference_file)) {
    baseline <- baseline_comparison_inputs()
    saveRDS(baseline$inputs, input_file)
    saveRDS(baseline$reference, reference_file)
    rm(baseline); gc(verbose = FALSE)
  }
  inputs <- readRDS(input_file)
  stopifnot(identical(inputs$baseline_md5, unname(tools::md5sum(baseline_file))),
    identical(inputs$INLA_version, as.character(packageVersion("INLA"))))
  full <- make_full_heat_lincombs(inputs$data, inputs$nonlinear_heat, restore_model_formula(inputs))
  stopifnot(nrow(full$combinations) == 42L)
  summaries <- list(variance_1000 = readRDS(reference_file))
  variances <- if (restart_100) 100 else c(100, 10)
  for (variance in variances) {
    scenario <- paste0("variance_", variance)
    scenario_dir <- if (restart_100) output_dir else file.path(output_dir, scenario)
    dir.create(scenario_dir, showWarnings = FALSE)
    model_file <- file.path(scenario_dir, "nbinomial_model.rds")
    fixed_prior <- inputs$fixed_prior
    fixed_prior$prec <- 1 / variance
    start <- NULL
    if (restart_100) {
      original <- readRDS(original_file)
      stopifnot(original$mode$mode.status == 0L, original$coefficient_sensitivity$variance == 100,
        identical(original$.args$data, inputs$data), identical(deparse(original$.args$formula), inputs$formula_text))
      start <- list(theta = original$mode$theta, x = original$mode$x, restart = TRUE, fixed = FALSE)
      stopifnot(length(start$theta) == 6L, all(is.finite(start$theta)), all(is.finite(start$x)))
      saveRDS(start, file.path(output_dir, "starting_solution.rds"))
      rm(original); gc(verbose = FALSE)
    }
    if (!file.exists(model_file)) {
      message("COEFFICIENT SENSITIVITY: variance ", variance, "; all PC priors unchanged; no draws.")
      fit <- fit_inla_model(inputs, scenario_dir, fixed_prior = fixed_prior, lincomb = full$lincombs, start = start)
      fit$nonlinear_heat <- inputs$nonlinear_heat
      fit$coefficient_sensitivity <- list(variance = variance, baseline_md5 = inputs$baseline_md5,
        elapsed_seconds = fit$analysis_elapsed_seconds, newly_generated_excess_rate_draws = 0L)
      if (restart_100) fit$warm_start_source <- list(file = original_file, md5 = unname(tools::md5sum(original_file)),
        theta_and_latent_modes = TRUE, restart_optimization = TRUE, hyperparameters_fixed = FALSE)
      saveRDS(fit, model_file)
    } else fit <- readRDS(model_file)
    expected_priors <- inputs$common_priors
    for (j in seq_along(expected_priors$fixed))
      if (expected_priors$fixed[[j]]$label != "(Intercept)") expected_priors$fixed[[j]]$prior.prec <- 1 / variance
    size_prior <- fit$all.hyper$family[[1L]]$hyper[[1L]]
    stopifnot(fit$mode$mode.status == 0L, identical(fit$names.fixed, inputs$names_fixed),
      identical(fit$.args$data, inputs$data), identical(as.numeric(fit$.args$E), as.numeric(inputs$exposure)),
      identical(deparse(fit$.args$formula), inputs$formula_text),
      identical(likelihood_common_priors(fit), expected_priors),
      identical(list(prior = as.character(size_prior$prior), param = as.numeric(size_prior$param),
        fixed = as.logical(size_prior$fixed)), inputs$size_prior),
      identical(fit$coefficient_sensitivity$baseline_md5, inputs$baseline_md5),
      fit$coefficient_sensitivity$variance == variance,
      all(is.finite(fit$summary.fixed$mean)), all(is.finite(fit$summary.hyperpar$mean)))
    if (restart_100) stopifnot(
      identical(fit$warm_start_source$md5, unname(tools::md5sum(original_file))),
      identical(as.numeric(fit$.args$control.mode$theta), as.numeric(start$theta)),
      identical(as.numeric(fit$.args$control.mode$x), as.numeric(start$x)),
      isTRUE(fit$.args$control.mode$restart), !isTRUE(fit$.args$control.mode$fixed),
      length(fit$dic$local.dic) == nrow(inputs$data), length(fit$waic$local.waic) == nrow(inputs$data),
      all(is.finite(fit$dic$local.dic)), all(is.finite(fit$waic$local.waic)),
      all(is.finite(as.matrix(fit$summary.fixed))), all(is.finite(as.matrix(fit$summary.hyperpar))),
      abs(sum(fit$dic$local.dic) - fit$dic$dic) < 1e-6,
      abs(sum(fit$waic$local.waic) - fit$waic$waic) < 1e-6)
    result <- c(list(zero = likelihood_zero_diagnostics(fit, inputs$data, "nbinomial")), summarise_heat_fit(fit, full))
    write_fit_details(fit, scenario_dir, result)
    if (restart_100) {
      write.csv(result$zero$grouped, file.path(output_dir, "zero_diagnostics_by_group.csv"), row.names = FALSE)
    } else {
      jsonlite::write_json(list(status = "completed", coefficient_variance = variance, mode_status = fit$mode$mode.status,
        fixed_coefficients = 26L, changed_non_intercept_priors = 25L, PC_priors_unchanged = TRUE,
        data_predictor_and_exposure_unchanged = TRUE, full_heat_contrasts = 42L, newly_generated_excess_rate_draws = 0L),
        file.path(scenario_dir, "validation.json"), pretty = TRUE, auto_unbox = TRUE)
    }
    summaries[[scenario]] <- result
    rm(fit, expected_priors, size_prior); gc(verbose = FALSE)
    stopifnot(identical(before, tools::md5sum(protected)))
  }
  if (restart_100) {
    original <- readRDS(file.path(original_dir, "summaries.rds"))
    baseline <- readRDS(file.path(comparison_dir, "baseline_variance_1000.rds"))
    results <- list(baseline_variance_1000 = baseline, original_variance_100 = original, rerun_variance_100 = result)
    comparison <- do.call(rbind, lapply(names(results), function(scenario)
      data.frame(scenario = scenario, results[[scenario]]$zero$metrics, row.names = NULL)))
    for (j in seq_along(results)) {
      index <- match(results[[j]]$heat$lincomb_id, original$heat$lincomb_id)
      stopifnot(!anyNA(index), length(index) == 42L)
      comparison$max_heat_IRR_change_pct_vs_original[j] <- 100 * max(abs(results[[j]]$heat$IRR_median / original$heat$IRR_median[index] - 1))
      comparison$max_heat_probability_change_pp_vs_original[j] <- 100 * max(abs(results[[j]]$heat$P_IRR_gt_1 - original$heat$P_IRR_gt_1[index]))
    }
    write.csv(comparison, file.path(output_dir, "comparison.csv"), row.names = FALSE)
    writeLines(c("VARIANCE 100: RESTART FROM ITS FITTED SOLUTION", "",
      "Same data, predictor and priors as the original variance-100 fit.",
      "Both saved theta and latent modes initialize the rerun; restart=TRUE and fixed=FALSE allow further optimization.",
      "No joint-sampling configurations, excess-rate posterior draws or map exports are generated.",
      "Original sensitivity results, NB baseline, current maps and methods files are preserved.",
      "The comparison quantifies reproducibility; finite output alone does not establish numerical stability.", "",
      capture.output(print(comparison, row.names = FALSE, digits = 9))), file.path(output_dir, "comparison.txt"))
    stopifnot(identical(before, tools::md5sum(protected)))
    jsonlite::write_json(list(status = "completed", mode_status = result$mode_status,
      same_model_and_priors = TRUE, warm_start_used = TRUE, optimization_restarted = TRUE,
      hyperparameters_fixed = FALSE, finite_pointwise_DIC_WAIC_checked = TRUE,
      protected_files_unchanged = TRUE, protected_file_count = length(before), newly_generated_excess_rate_draws = 0L,
      original_model_md5 = unname(tools::md5sum(original_file)), baseline_model_md5 = inputs$baseline_md5),
      file.path(output_dir, "validation.json"), pretty = TRUE, auto_unbox = TRUE, digits = 16)
    print(comparison, row.names = FALSE)
    return(invisible(comparison))
  }
  combine <- function(fun) do.call(rbind, lapply(names(summaries), function(scenario)
    data.frame(scenario = scenario, fun(summaries[[scenario]]), row.names = NULL)))
  comparison <- combine(function(x) x$zero$metrics)
  comparison$coefficient_variance <- c(1000, 100, 10)
  baseline <- summaries$variance_1000
  comparison$delta_DIC_vs_baseline <- comparison$DIC - baseline$zero$metrics$DIC
  comparison$delta_WAIC_vs_baseline <- comparison$WAIC - baseline$zero$metrics$WAIC
  for (j in seq_along(summaries)) {
    heat <- summaries[[j]]$heat
    index <- match(heat$lincomb_id, baseline$heat$lincomb_id)
    stopifnot(length(index) == 42L, !anyNA(index))
    comparison$max_heat_IRR_change_pct[j] <- 100 * max(abs(heat$IRR_median / baseline$heat$IRR_median[index] - 1))
    comparison$max_heat_probability_change_pp[j] <- 100 * max(abs(heat$P_IRR_gt_1 - baseline$heat$P_IRR_gt_1[index]))
    comparison$heat_CI_exclusion_changes[j] <- sum((heat$IRR_lower > 1 | heat$IRR_upper < 1) !=
      (baseline$heat$IRR_lower[index] > 1 | baseline$heat$IRR_upper[index] < 1))
  }
  stopifnot(all(is.finite(comparison$DIC)), all(is.finite(comparison$WAIC)))
  write.csv(comparison, file.path(output_dir, "comparison.csv"), row.names = FALSE)
  for (name in c("heat", "fixed", "hyper")) write.csv(combine(function(x) x[[name]]),
    file.path(output_dir, paste0(name, "_comparison.csv")), row.names = FALSE)
  write.csv(combine(function(x) x$zero$grouped), file.path(output_dir, "zero_diagnostics_by_group.csv"), row.names = FALSE)
  write.csv(combine(function(x) x$zero$count_bins), file.path(output_dir, "count_distribution_diagnostics.csv"), row.names = FALSE)
  stopifnot(identical(before, tools::md5sum(protected)))
  jsonlite::write_json(list(status = "completed", variances = c(1000, 100, 10), baseline_reused = TRUE,
    baseline_model_md5 = inputs$baseline_md5, protected_files_unchanged = TRUE, protected_file_count = length(before),
    newly_generated_excess_rate_draws = 0L, observations = 54910L, ZIPs = 309L,
    INLA_version = as.character(packageVersion("INLA"))),
    file.path(output_dir, "validation.json"), pretty = TRUE, auto_unbox = TRUE, digits = 16)
  writeLines(c("NB COEFFICIENT-PRIOR VARIANCE SENSITIVITY", "",
    "Variance 1000 is the exact saved baseline; only variances 100 and 10 are newly fitted.",
    "Only the Normal variances of the 25 non-intercept fixed coefficients change; prior means remain zero.",
    "The flat intercept, four PC precision priors, exchangeable PC correlation and PC NB-size priors are unchanged.",
    "Same observations, population exposures, heat grid, shared/zone heat curves and ZIP iid intercept.",
    "No posterior draws or maps are generated. Baseline results, current maps and methods are hash-verified unchanged.",
    "Smaller DIC/WAIC is preferred. Heat changes use the complete 42 cluster-zone contrasts.",
    "Zero checks use posterior mean fitted rates and posterior mean NB size (in-sample plug-in, not fully integrated predictive tests).", "",
    capture.output(print(comparison, row.names = FALSE, digits = 8))), file.path(output_dir, "comparison.txt"))
  print(comparison, row.names = FALSE)
  invisible(comparison)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  supported_options <- c("--preflight", "--fit-only", "--outputs-only", "--validate-only", "--sensitivity", "--likelihood-comparison", "--coefficient-sensitivity", "--rerun-variance-100", "--help")
  if (length(args) > 1L || any(!args %in% supported_options))
    stop("Choose one supported option; run --help for the current model and sensitivity commands.")
  if ("--help" %in% args) {
    cat("Run from the project root:\n",
      "  Rscript hri_ed_analysis.R                  Fit, sample, export, and validate\n",
      "  Rscript hri_ed_analysis.R --preflight      Check the model specification\n",
      "  Rscript hri_ed_analysis.R --fit-only       Fit and save model summaries\n",
      "  Rscript hri_ed_analysis.R --outputs-only   Reuse the saved fit for outputs\n",
      "  Rscript hri_ed_analysis.R --validate-only  Check the saved model and outputs\n", sep = "")
    cat("  Rscript hri_ed_analysis.R --sensitivity    Compare priors using the corrected model\n")
    cat("  Rscript hri_ed_analysis.R --likelihood-comparison  Matched Poisson/NB and zero checks, no draws\n")
    cat("  Rscript hri_ed_analysis.R --coefficient-sensitivity  Normal coefficient variances 100 and 10, no draws\n")
    cat("  Rscript hri_ed_analysis.R --rerun-variance-100  Restart variance 100 from its fitted solution, no draws\n")
  } else if ("--rerun-variance-100" %in% args) {
    run_coefficient_sensitivity(restart_100 = TRUE)
  } else if ("--coefficient-sensitivity" %in% args) {
    run_coefficient_sensitivity()
  } else if ("--likelihood-comparison" %in% args) {
    run_likelihood_comparison()
  } else if ("--sensitivity" %in% args) {
    run_nb_prior_sensitivity()
  } else if ("--validate-only" %in% args) {
    validate_nb_pc_model()
  } else if ("--outputs-only" %in% args) {
    export_preserved_nb_outputs()
  } else {
    fit_preserved_nb(preflight = "--preflight" %in% args, fit_only = "--fit-only" %in% args)
  }
}
