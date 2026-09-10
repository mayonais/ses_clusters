# Posterior inference for the signed cluster-average excess rate per 10,000.
# data must retain the fitted model's row order. The reference heat is the
# cluster mean of chronic_heat, with equal weight for each observed ZIP-day.
cluster_ar_exceedance <- function(model, data, reference_cluster,
                                  n_samples = 4000L, batch_size = 200L,
                                  seed = 20260910L) {
  stopifnot(n_samples > 0L, n_samples == as.integer(n_samples),
            batch_size > 0L, batch_size == as.integer(batch_size),
            seed > 0L, seed + n_samples < .Machine$integer.max,
            nrow(data) == nrow(model$summary.linear.predictor),
            all(is.finite(data$heat_z)), all(is.finite(data$chronic_heat)),
            !anyNA(data$cluster))
  clusters <- levels(droplevels(factor(data$cluster)))
  reference_cluster <- as.character(reference_cluster)
  stopifnot(reference_cluster %in% clusters)
  cluster_index <- match(as.character(data$cluster), clusters)
  cluster_n <- tabulate(cluster_index, nbins = length(clusters))
  H_bar <- as.numeric(rowsum(data$chronic_heat, cluster_index, reorder = TRUE)) /
    cluster_n
  heat_delta <- data$heat_z - H_bar[cluster_index]

  interaction_terms <- paste0("cluster", clusters[clusters != reference_cluster],
                              ":acute_heat")
  slope_terms <- c("acute_heat", interaction_terms)
  stopifnot(all(slope_terms %in% model$names.fixed))
  selection <- c(list(Predictor = seq_len(nrow(data))),
                 setNames(rep(list(1L), length(slope_terms)), slope_terms))
  cluster_AR <- matrix(NA_real_, nrow = length(clusters), ncol = n_samples,
                       dimnames = list(paste0("cluster", clusters), NULL))
  set.seed(seed)
  for (start in seq.int(1L, n_samples, by = batch_size)) {
    end <- min(start + batch_size - 1L, n_samples)
    samples <- INLA::inla.posterior.sample(
      end - start + 1L, model, selection = selection,
      seed = seed + start - 1L, num.threads = "1:1", add.names = FALSE)
    latent_names <- rownames(samples[[1]]$latent)
    predictor_index <- match(paste0("Predictor:", seq_len(nrow(data))), latent_names)
    shared_index <- match("acute_heat:1", latent_names)
    deviation_index <- match(paste0(interaction_terms, ":1"), latent_names)
    stopifnot(!anyNA(predictor_index), !is.na(shared_index), !anyNA(deviation_index))
    for (j in seq_along(samples)) {
      latent <- samples[[j]]$latent
      beta <- rep(latent[shared_index], length(clusters))
      beta[clusters != reference_cluster] <- beta[clusters != reference_cluster] +
        latent[deviation_index]
      eta <- latent[predictor_index]
      # Signed contrast: aggregate positive AND negative ZIP-day contributions.
      # expm1 retains precision for temperatures close to the reference.
      AR <- 10000 * exp(eta) * -expm1(-beta[cluster_index] * heat_delta)
      stopifnot(all(is.finite(AR)))
      cluster_AR[, start + j - 1L] <-
        as.numeric(rowsum(AR, cluster_index, reorder = TRUE)) / cluster_n
    }
    rm(samples)
  }
  n_positive <- rowSums(cluster_AR > 0)
  intervals <- t(apply(cluster_AR, 1L, quantile, probs = c(0.025, 0.975)))
  summary <- data.frame(
    cluster = rownames(cluster_AR), n_zip_days = cluster_n,
    H_bar = H_bar, AR_mean = rowMeans(cluster_AR),
    AR_q025 = intervals[, 1L], AR_q975 = intervals[, 2L],
    P_AR_gt_0 = n_positive / n_samples,
    n_positive = n_positive, n_samples = n_samples,
    row.names = NULL)
  list(summary = summary, draws = cluster_AR,
       settings = list(seed = seed, n_samples = n_samples, batch_size = batch_size,
                       weighting = "equal ZIP-day", counterfactual = "cluster-mean heat",
                       units = "excess rate per 10000 population per day",
                       INLA_version = as.character(utils::packageVersion("INLA"))))
}
