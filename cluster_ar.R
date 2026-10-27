# Posterior inference for the signed cluster-average excess rate per 10,000.
# Data retain the fitted row order; ZIP-days receive equal weight.
# Express each SES cluster-month's mean heat on the overall ZIP-centered scale.
cluster_heat_reference <- function(data) {
  stopifnot(all(c("cluster", "month", "heat_z", "chronic_heat") %in% names(data)),
    !anyNA(data[c("cluster", "month", "heat_z", "chronic_heat")]))
  as.numeric(ave(data$heat_z, data$cluster, data$month, FUN = mean)) - data$chronic_heat
}

cluster_ar_exceedance <- function(model, data, reference_cluster,
                                  n_samples = 2000L, batch_size = 50L,
                                  seed = 20260910L, nonlinear_heat = NULL,
                                  sampling_model_file = NULL) {
  stopifnot(n_samples > 0L, n_samples == as.integer(n_samples),
            batch_size > 0L, batch_size == as.integer(batch_size),
            seed > 0L, seed + n_samples < .Machine$integer.max,
            nrow(data) == if (!is.null(model$n_observations))
              model$n_observations else nrow(model$summary.linear.predictor),
            all(is.finite(data$heat_z)), all(is.finite(data$chronic_heat)),
            all(is.finite(data$acute_heat)),
            max(abs(data$acute_heat - (data$heat_z - data$chronic_heat))) < 1e-10,
            !anyNA(data$cluster))
  clusters <- levels(droplevels(factor(data$cluster)))
  reference_cluster <- as.character(reference_cluster)
  stopifnot(reference_cluster %in% clusters)
  cluster_index <- match(as.character(data$cluster), clusters)
  cluster_n <- tabulate(cluster_index, nbins = length(clusters))
  reference_acute <- if (identical(nonlinear_heat$ar_reference, "SES-cluster-month mean heat"))
    cluster_heat_reference(data) else rep(0, nrow(data))
  heat_delta <- data$acute_heat - reference_acute
  mean_reference_heat_z <- as.numeric(rowsum(data$chronic_heat + reference_acute,
    cluster_index, reorder = TRUE)) / cluster_n
  
  interaction_clusters <- clusters[clusters != reference_cluster]
  has_interaction <- paste0("cluster", interaction_clusters, ":acute_heat") %in% model$names.fixed
  interaction_clusters <- interaction_clusters[has_interaction]
  interaction_terms <- paste0("cluster", interaction_clusters, ":acute_heat")
  if (!length(interaction_clusters)) interaction_terms <- character()
  has_fixed_heat <- "acute_heat" %in% model$names.fixed
  slope_terms <- c(if (has_fixed_heat) "acute_heat", interaction_terms)
  stopifnot(all(slope_terms %in% model$names.fixed))
  selection <- c(list(Predictor = seq_len(nrow(data))),
                 setNames(rep(list(1L), length(slope_terms)), slope_terms))
  if (!is.null(nonlinear_heat)) {
    n_grid <- length(nonlinear_heat$grid)
    shared_field <- nonlinear_heat$shared_field
    has_cz <- !is.na(nonlinear_heat$zone) & nzchar(nonlinear_heat$zone)
    pooled <- !is.null(nonlinear_heat$deviation_field)
    zone_index <- match(nonlinear_heat$zone, if (pooled)
      nonlinear_heat$zone_levels else names(nonlinear_heat$fields))
    stopifnot(length(zone_index) == nrow(data), all(!has_cz | !is.na(zone_index)), length(nonlinear_heat$observed_index) == nrow(data), length(nonlinear_heat$counterfactual_index) == nrow(data))
    stopifnot(all(nonlinear_heat$observed_index %in% seq_len(n_grid)), all(nonlinear_heat$counterfactual_index %in% seq_len(n_grid)))
    stopifnot(all(reference_acute >= min(nonlinear_heat$grid)),
      all(reference_acute <= max(nonlinear_heat$grid)))
    reference_index <- findInterval(reference_acute, (head(nonlinear_heat$grid, -1L) +
      tail(nonlinear_heat$grid, -1L)) / 2) + 1L
    stopifnot(all(nonlinear_heat$counterfactual_index == reference_index))
    if (pooled) {
      selection <- c(selection, setNames(list(seq_len(n_grid)), shared_field),
        setNames(list(seq_len(n_grid * length(nonlinear_heat$zone_levels))),
          nonlinear_heat$deviation_field))
    } else {
      nonlinear_fields <- c(shared_field, unname(nonlinear_heat$fields))
      selection <- c(selection,
        setNames(rep(list(seq_len(n_grid)), length(nonlinear_fields)), nonlinear_fields))
    }
  }
  cluster_AR <- matrix(NA_real_, nrow = length(clusters), ncol = n_samples,
                       dimnames = list(paste0("cluster", clusters), NULL))
  set.seed(seed)
  sampler <- heat_posterior_sampler(model, model_file = sampling_model_file)
  on.exit(sampler$close(), add = TRUE)
  for (start in seq.int(1L, n_samples, by = batch_size)) {
    end <- min(start + batch_size - 1L, n_samples)
    message("AR posterior draws ", start, "-", end, " of ", n_samples)
    samples <- sampler$sample(end - start + 1L, selection = selection,
      seed = seed + start - 1L)
    latent_names <- rownames(samples[[1]]$latent)
    predictor_index <- match(paste0("Predictor:", seq_len(nrow(data))), latent_names)
    shared_index <- if (has_fixed_heat) match("acute_heat:1", latent_names) else integer()
    deviation_index <- if (length(interaction_terms))
      match(paste0(interaction_terms, ":1"), latent_names) else integer()
    stopifnot(!anyNA(predictor_index), !anyNA(shared_index), !anyNA(deviation_index))
    if (!is.null(nonlinear_heat)) {
      if (!is.null(shared_field)) {
        shared_rw2 <- match(paste0(shared_field, ":", seq_len(n_grid)), latent_names)
        stopifnot(!anyNA(shared_rw2))
      }
      rw2_index <- if (pooled) {
        matrix(match(paste0(nonlinear_heat$deviation_field, ":",
          seq_len(n_grid * length(nonlinear_heat$zone_levels))), latent_names),
          nrow = n_grid)
      } else {
        vapply(nonlinear_heat$fields, function(field) match(paste0(field, ":", seq_len(n_grid)), latent_names), integer(n_grid))
      }
      stopifnot(!anyNA(rw2_index))
      observed_rw2 <- counterfactual_rw2 <- rep(NA_integer_, nrow(data))
      observed_rw2[has_cz] <- rw2_index[cbind(nonlinear_heat$observed_index[has_cz], zone_index[has_cz])]
      counterfactual_rw2[has_cz] <- rw2_index[cbind(nonlinear_heat$counterfactual_index[has_cz], zone_index[has_cz])]
    }
    for (j in seq_along(samples)) {
      latent <- samples[[j]]$latent
      beta <- if (has_fixed_heat) rep(latent[shared_index], length(clusters)) else numeric(length(clusters))
      if (length(interaction_clusters)) {
        beta[clusters %in% interaction_clusters] <- beta[clusters %in% interaction_clusters] +
          latent[deviation_index]
      }
      eta <- latent[predictor_index]
      # Signed contrast: aggregate positive AND negative ZIP-day contributions.
      # expm1 retains precision for temperatures close to the reference.
      log_rate_contrast <- beta[cluster_index] * heat_delta
      if (!is.null(nonlinear_heat)) {
        nonlinear_contrast <- numeric(nrow(data))
        # Shared curve affects every ZIP, including those with unknown zones.
        if (!is.null(shared_field)) {
          nonlinear_contrast <-
            latent[shared_rw2[nonlinear_heat$observed_index]] -
            latent[shared_rw2[nonlinear_heat$counterfactual_index]]
        }
        # Known zones also receive their climate-zone-specific deviation.
        nonlinear_contrast[has_cz] <- nonlinear_contrast[has_cz] +
          latent[observed_rw2[has_cz]] - latent[counterfactual_rw2[has_cz]]
        log_rate_contrast <- log_rate_contrast + nonlinear_contrast
      }
      AR <- 10000 * exp(eta) * -expm1(-log_rate_contrast)
      stopifnot(all(is.finite(AR)))
      cluster_AR[, start + j - 1L] <-
        as.numeric(rowsum(AR, cluster_index, reorder = TRUE)) / cluster_n
    }
    rm(samples)
    gc(verbose = FALSE)
  }
  n_positive <- rowSums(cluster_AR > 0)
  intervals <- t(apply(cluster_AR, 1L, quantile, probs = c(0.025, 0.975)))
  summary <- data.frame(
    cluster = rownames(cluster_AR), n_zip_days = cluster_n,
    mean_reference_heat_z = mean_reference_heat_z, AR_mean = rowMeans(cluster_AR),
    AR_q025 = intervals[, 1L], AR_q975 = intervals[, 2L],
    P_AR_gt_0 = n_positive / n_samples,
    n_positive = n_positive, n_samples = n_samples,
    row.names = NULL)
  list(summary = summary, draws = cluster_AR,
       settings = list(seed = seed, n_samples = n_samples, batch_size = batch_size,
                       weighting = "equal ZIP-day", counterfactual = if (is.null(nonlinear_heat$ar_reference)) "ZIP-mean heat" else nonlinear_heat$ar_reference,
                       reference_acute_heat = reference_acute,
                       units = "excess rate per 10000 population per day",
                       INLA_version = as.character(utils::packageVersion("INLA")), nonlinear_heat = nonlinear_heat))
}
