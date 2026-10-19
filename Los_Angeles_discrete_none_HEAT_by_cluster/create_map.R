# Run from this input folder to redraw the three probability maps.
local({
  original_dir <- getwd()
  on.exit(setwd(original_dir))
  setwd("..")
# Redraw all three mainland probability maps with identical panel dimensions.
library(ggplot2)
source("shared_functions.R")
root <- "Los_Angeles_discrete_none_HEAT_by_cluster"
assignments <- readRDS(file.path(root, "Los Angeles_ACS_zcta_clustered_for_ED.rds"))
if (inherits(assignments, "sf")) assignments <- sf::st_drop_geometry(assignments)
assignments$cluster <- as.character(assignments$cluster)
geometry <- mainland_probability_geometry(
  readRDS(file.path(root, "acute_heat_poisson/los_angeles_zcta_2020.rds")), assignments)
plots <- list()
paths <- character()
for (folder in c("acute_heat_poisson", "acute_heat_nb")) {
  output <- file.path(root, folder)
  summary <- read.csv(file.path(output, "acute_heat_cluster_summary.csv"))
  estimates <- data.frame(cluster = as.character(summary$cluster), probability = summary$P_beta_gt_0)
  model <- if (folder == "acute_heat_poisson") "Poisson" else "Negative-binomial"
  plot <- draw_probability_map(geometry, estimates,
     expression(Pr(beta[acute*","*k] > 0 ~ "|" ~ Y)),
     paste0(model, ": positive acute-heat coefficient\nCluster-specific acute-heat slope\nMainland Los Angeles, May-October 2018"))
  saved <- readRDS(file.path(output, "acute_heat_ggplots.rds"))
  saved$probability_map <- plot
  saveRDS(saved, file.path(output, "acute_heat_ggplots.rds"))
  plots[[folder]] <- plot
  map_output <- file.path("maps", folder)
  dir.create(map_output, recursive = TRUE, showWarnings = FALSE)
  paths <- c(paths, file.path(map_output, "acute_heat_probability_map_zcta"))
}
ar <- read.csv(file.path(root, "Los Angeles_poisson_cluster_AR_exceedance.csv"))
estimates <- data.frame(cluster = sub("^cluster", "", ar$cluster), probability = ar$P_AR_gt_0)
plots$AR <- draw_probability_map(geometry, estimates,
  expression(pi[c] == Pr(AR[c](theta) > 0 ~ "|" ~ Y)),
  "Poisson: positive cluster-average excess rate\nSigned mean rate contrast against the cluster-mean heat reference\nMainland Los Angeles, May-October 2018")
dir.create("maps", showWarnings = FALSE)
paths <- c(paths, file.path("maps", "AR_gt_0_gradient"))
saveRDS(plots$AR, file.path(root, "AR_gt_0_gradient.rds"))
grobs <- lapply(plots, ggplotGrob)
widths <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$widths))
heights <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$heights))
for (i in seq_along(grobs)) {
  grobs[[i]]$widths <- widths
  grobs[[i]]$heights <- heights
  for (ext in c("png", "pdf")) {
    save_project_map(paste0(paths[i], ".", ext), grobs[[i]])
  }
}
cat("All three maps: mainland geometry, identical 10 x 8 inch layout, no footer, fixed 0-1 blue-white-red scale centered at 0.5.\n")

})
