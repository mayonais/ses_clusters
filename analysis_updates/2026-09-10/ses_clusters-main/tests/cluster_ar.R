# Run from ses_clusters-main: Rscript tests/cluster_ar.R
library(INLA)
source("cluster_ar.R")
set.seed(42)
data <- data.frame(
  cluster = factor(rep(c("2", "1"), each = 80), levels = c("2", "1")),
  acute_heat = rep(seq(-2, 2, length.out = 80), 2))
data$heat_z <- data$acute_heat
data$chronic_heat <- 0
data$cold <- as.numeric(data$acute_heat < 0)
data$y <- rpois(nrow(data), exp(1 + 1.4 * data$cold + 0.15 * data$acute_heat))
model <- inla(y ~ cluster + acute_heat + acute_heat:cluster + cold,
  family = "poisson", data = data,
  control.compute = list(config = TRUE),
  control.predictor = list(compute = TRUE), num.threads = "1:1")

n <- 400L
batch <- 100L
seed <- 123L
result <- cluster_ar_exceedance(model, data, "2", n, batch, seed)
# Independently compute signed contrasts from the same joint posterior draws.
# This checks the returned distribution, not only a probability close to 0 or 1.
selection <- list(Predictor = seq_len(nrow(data)), acute_heat = 1L,
                  "cluster1:acute_heat" = 1L)
manual <- matrix(NA_real_, 2, n)
clipped <- matrix(NA_real_, 2, n)
daily_fraction <- matrix(NA_real_, 2, n)
set.seed(seed)
for (start in seq.int(1L, n, by = batch)) {
  samples <- inla.posterior.sample(batch, model, selection = selection,
    seed = seed + start - 1L, num.threads = "1:1", add.names = TRUE)
  for (j in seq_along(samples)) {
    latent <- samples[[j]]$latent
    eta <- latent[paste0("Predictor:", seq_len(nrow(data))), 1L]
    beta <- latent["acute_heat:1", 1L] +
      ifelse(data$cluster == "1", latent["cluster1:acute_heat:1", 1L], 0)
    signed <- 10000 * (exp(eta) - exp(eta - beta * data$heat_z))
    for (c in seq_along(levels(data$cluster))) {
      rows <- data$cluster == levels(data$cluster)[c]
      manual[c, start + j - 1L] <- mean(signed[rows])
      clipped[c, start + j - 1L] <- mean(pmax(signed[rows], 0))
      daily_fraction[c, start + j - 1L] <- mean(signed[rows] > 0)
    }
  }
}
stopifnot(max(abs(result$draws - manual)) < 1e-8,
          identical(result$summary$P_AR_gt_0, unname(rowMeans(manual > 0))),
          any(abs(rowMeans(manual > 0) - rowMeans(daily_fraction)) > 0.1),
          any(abs(rowMeans(manual > 0) - rowMeans(clipped > 0)) > 0.1),
          any(manual < 0), any(manual > 0))

# A reference-temperature contrast is exactly zero: strict exceedance is false.
zero_data <- data
zero_data$heat_z <- 0
zero <- cluster_ar_exceedance(model, zero_data, "2", 20L, 10L, seed)
stopifnot(all(zero$draws == 0), all(zero$summary$P_AR_gt_0 == 0))

# The short final batch must be included, and a fixed seed must be reproducible.
a <- cluster_ar_exceedance(model, data, "2", 23L, 10L, seed)
b <- cluster_ar_exceedance(model, data, "2", 23L, 10L, seed)
stopifnot(identical(a$draws, b$draws), ncol(a$draws) == 23L,
          all(a$summary$n_samples == 23L))
print(result$summary[, c("cluster", "AR_mean", "P_AR_gt_0")])
cat("PASS: joint signed cluster contrasts, non-6 reference, zero contrast, batching and reproducibility.\n")
