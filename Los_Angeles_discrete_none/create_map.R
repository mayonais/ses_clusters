# Run from this input folder to redraw the three probability maps.
local({
  original_dir <- getwd()
  on.exit(setwd(original_dir))
  script_paths <- as.character(unlist(lapply(sys.frames(), function(frame) frame$ofile)))
  script_paths <- script_paths[basename(script_paths) == "create_map.R"]
  command_path <- grep("^--file=", commandArgs(), value = TRUE)
  script_path <- if (length(script_paths)) tail(script_paths, 1) else
    if (length(command_path)) sub("^--file=", "", command_path[1]) else character()
  project_dir <- if (length(script_path)) {
    normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
  } else if (file.exists("shared_functions.R")) getwd() else normalizePath("..")
  setwd(project_dir)
  if (!interactive() && grDevices::dev.cur() == 1L) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
# Redraw all three mainland probability maps with identical panel dimensions.
library(ggplot2)
source("shared_functions.R")
model_variant <- "climate_zone"  # "climate_zone", "climate_zone_added", or "linear_cluster"
model_family <- "poisson"  # "poisson", "nbinomial", or "both"
root <- file.path("Los_Angeles_discrete_none")
results_root <- if (model_variant %in% c("climate_zone_added", "climate_zone")) file.path(root, model_variant) else file.path(root, "linear_cluster")
map_root <- file.path("maps", model_variant)
assignments <- readRDS(file.path(root, "Los Angeles_ACS_zcta_clustered_for_ED.rds"))
if (inherits(assignments, "sf")) assignments <- sf::st_drop_geometry(assignments)
assignments$cluster <- as.character(assignments$cluster)
geometry <- mainland_probability_geometry(
  readRDS(file.path(root, "linear_cluster", "acute_heat_poisson/los_angeles_zcta_2020.rds")), assignments)
plots <- list()
paths <- character()
probability_folders <- switch(model_family,
  poisson = "acute_heat_poisson", nbinomial = "acute_heat_nb",
  both = c("acute_heat_poisson", "acute_heat_nb"), stop("Unknown model_family"))
if (model_variant %in% c("climate_zone_added", "climate_zone")) {
  # Full-response acute-heat maps share posterior draws and aligned panels.
  # Delegate both exports to their existing script so this entry point cannot
  # replace the nonlinear probability map with a linear-slope probability.
  full_heat_environment <- new.env(parent = environment())
  full_heat_environment$irr_model_families <- switch(model_family,
    poisson = "poisson", nbinomial = "nbinomial", both = c("poisson", "nbinomial"))
  full_heat_environment$model_variant <- model_variant
  source("acute_heat_irr_maps.R", local = full_heat_environment)
  probability_folders <- character()
}
for (folder in probability_folders) {
  output <- file.path(results_root, folder)
  summary <- read.csv(file.path(output, "acute_heat_cluster_summary.csv"))
  estimates <- data.frame(cluster = as.character(summary$cluster), probability = summary$P_beta_gt_0)
  model <- if (folder == "acute_heat_poisson") "Poisson" else "Negative-binomial"
  plot <- draw_probability_map(geometry, estimates,
     expression(Pr(beta[acute*","*k] > 0 ~ "|" ~ Y)),
     if (model_variant %in% c("climate_zone_added", "climate_zone")) {
       paste0(model, " + shared and pooled climate-zone RW2 (nonlinear addition)\n",
         "Positive linear acute-heat component, adjusted for RW2 curves\n",
         "Mainland Los Angeles, May-October 2018")
     } else paste0(model, ": positive acute-heat coefficient\n",
       "Cluster-specific acute-heat slope\nMainland Los Angeles, May-October 2018"))
  saved <- readRDS(file.path(output, "acute_heat_ggplots.rds"))
  saved$probability_map <- plot
  saveRDS(saved, file.path(output, "acute_heat_ggplots.rds"))
  plots[[folder]] <- plot
  map_output <- file.path(map_root, folder)
  dir.create(map_output, recursive = TRUE, showWarnings = FALSE)
  paths <- c(paths, file.path(map_output, "acute_heat_probability_map_zcta"))
}
ar <- read.csv(file.path(results_root, "Los Angeles_poisson_cluster_AR_exceedance.csv"))
estimates <- data.frame(cluster = sub("^cluster", "", ar$cluster), probability = ar$P_AR_gt_0)
plots$AR <- draw_probability_map(geometry, estimates,
  expression(Pr(AR[c](theta) > 0 ~ "|" ~ Y)),
  if (model_variant %in% c("climate_zone_added", "climate_zone")) {
    paste0("Poisson model |\n",
      "Cluster-average excess ED rate: observed heat versus\n",
      "each SES cluster\u2019s monthly 2018 mean heat")
  } else "Poisson: positive cluster-average excess rate\nSigned mean rate contrast against each ZIP's 2018 mean heat\nMainland Los Angeles, May-October 2018") +
  labs(title = "Probability of Positive Cluster-Average Excess ED Rate")
ar_probability_scale <- plots$AR$scales$get_scales("fill")
ar_probability_scale$labels <- scales::label_percent(accuracy = 1)
dir.create(map_root, recursive = TRUE, showWarnings = FALSE)
paths <- c(paths, file.path(map_root, "AR_gt_0_gradient"))
saveRDS(plots$AR, file.path(results_root, "AR_gt_0_gradient.rds"))
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
cat("Probability maps: mainland geometry, aligned 10 x 8 inch layout, fixed 0-1 blue-white-red scale centered at 0.5.\n")

})
