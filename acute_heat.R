# Choose "poisson", "nbinomial", or "both"; then source this whole script.
model_family <- "poisson"
model_variant <- "climate_zone"  # "climate_zone", "climate_zone_added", or "linear_cluster"
# Linear maps use exp(beta_acute,k). Nonlinear maps use the complete heat
# contrast (ZIP mean + one global heat SD versus ZIP mean), including both RW2s.
library(dplyr)
library(ggplot2)
library(INLA)

script_files <- as.character(unlist(lapply(sys.frames(), function(frame) frame$ofile)))
script_files <- script_files[basename(script_files) == "acute_heat.R"]
command_file <- grep("^--file=", commandArgs(), value = TRUE)
project_dir <- if (length(script_files)) {
  dirname(normalizePath(tail(script_files, 1), mustWork = TRUE))
} else if (length(command_file) && sub("^--file=", "", command_file[1]) != "-") {
  dirname(normalizePath(sub("^--file=", "", command_file[1]), mustWork = TRUE))
} else {
  getwd() 
}

# ---- Settings ---------------------------------------------------------------
study_name <- "Los Angeles"
input_dir <- file.path(project_dir, "Los_Angeles_discrete_none")
model_prefix <- "Los_Angeles_discrete_none_HEAT_by_cluster"
model_families <- switch(model_family, both = c("poisson", "nbinomial"),
                         poisson = "poisson", nbinomial = "nbinomial",
                         stop('model_family must be "poisson", "nbinomial", or "both".'))
model_cache_dir <- file.path(input_dir, "model_cache")
results_root <- if (model_variant %in% c("climate_zone_added", "climate_zone")) {
  file.path(input_dir, model_variant)
} else file.path(input_dir, "linear_cluster")
# FALSE reuses a compatible fit and refits automatically if its signature changes.
# Set TRUE only when you want to force a fresh fit, independently of plot changes.
force_refit <- FALSE
geometry_year <- 2020L
source(file.path(project_dir, "shared_functions.R"), local = TRUE)

# ---- Shared analysis --------------------------------------------------------
prepared <- load_heat_data(input_dir, study_name)
clustered_zctas <- prepared$clustered_zctas
cluster_ids <- levels(prepared$daily_ed$cluster)
for (family in model_families) {
model_label <- if (family == "poisson") "Poisson" else "Negative-binomial"
output_dir <- file.path(results_root, if (family == "poisson") "acute_heat_poisson" else "acute_heat_nb")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
map_output_dir <- file.path(project_dir, "maps", model_variant, basename(output_dir))
dir.create(map_output_dir, recursive = TRUE, showWarnings = FALSE)
if (model_variant %in% c("climate_zone_added", "climate_zone")) {
  model_suffix <- if (family == "poisson") "poisson" else "nb"
  model_file <- file.path(results_root, paste0(model_prefix,
    "_suppressed_INLA_", model_suffix, "_model.rds"))
  if (!file.exists(model_file)) stop("Run hri_ed_analysis.R first: missing ", model_file)
  acute_fit <- readRDS(model_file)
  stopifnot(!is.null(acute_fit$nonlinear_heat$shared_field),
    isTRUE(acute_fit$nonlinear_heat$pooled),
    nrow(acute_fit$summary.linear.predictor) == nrow(prepared$daily_ed))
} else {
  acute_fit <- fit_heat_model(prepared, model_cache_dir, family = family,
                              force_refit = force_refit)
}
if (!"acute_heat" %in% acute_fit$names.fixed) {
  stopifnot(model_variant %in% c("climate_zone_added", "climate_zone"))
  message("No fixed heat slope: exporting RW2 IRRs and probabilities below.")
  next
}
acute_summary <- summarise_heat_slopes(acute_fit, cluster_ids) %>%
  left_join(prepared$cluster_labels, by = "cluster") %>%
  arrange(as.numeric(cluster)) %>%
  mutate(cluster_display = factor(paste("Cluster", cluster),
                                  levels = rev(paste("Cluster", cluster))))
write.csv(acute_summary %>% select(-cluster_display),
          file.path(output_dir, "acute_heat_cluster_summary.csv"), row.names = FALSE)

# ---- Plots ------------------------------------------------------------------
# RR point estimate: posterior median; intervals: equal-tail 95% credible intervals.
rr_plot <- ggplot(acute_summary, aes(y = cluster_display)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_segment(aes(x = RR_low, xend = RR_high, yend = cluster_display), linewidth = 0.6) +
  geom_point(aes(x = RR_median), size = 2.8, colour = "#2166AC") +
  scale_x_log10() +
  labs(title = "Cluster-specific acute-heat relative risks",
       subtitle = paste0(model_label, " model: posterior medians and 95% credible intervals"),
       x = expression("Incidence rate ratio " * exp(beta[acute*","*k]) * " (log scale)"), y = NULL,
       caption = paste0("Per one SD increase in acute heat, using the SD of heat exposure across all LA ZCTA-day observations,\n",
                        "holding chronic heat fixed. Dashed line: IRR = 1.")) +
  theme_minimal(base_size = 11)

save_plot <- function(plot, name, width = 9, height = 7) {
  ggsave(file.path(output_dir, paste0(name, ".png")), plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, paste0(name, ".pdf")), plot, width = width, height = height, bg = "white")
}
if (model_variant == "linear_cluster") save_plot(rr_plot, "acute_heat_relative_risk")

# ---- Probability choropleth: same ZCTA vintage as the existing map ------------
geometry_file <- file.path(input_dir, "linear_cluster", basename(output_dir),
                          paste0("los_angeles_zcta_", geometry_year, ".rds"))
if (file.exists(geometry_file)) {
  map_geometry <- readRDS(geometry_file)
} else {
  message("Loading Census county and ZCTA geometry (first run may download files)...")
  options(tigris_use_cache = TRUE)
  boundary <- tigris::counties(state = "CA", cb = TRUE, year = geometry_year) %>%
    filter(NAME == study_name) %>% sf::st_make_valid() %>% sf::st_transform(3310)
  zctas <- tigris::zctas(cb = TRUE, year = geometry_year) %>%
    sf::st_make_valid() %>% sf::st_transform(3310)
  zctas <- sf::st_filter(zctas, boundary)
  zctas <- suppressWarnings(sf::st_intersection(zctas, sf::st_geometry(boundary)))
  map_geometry <- list(boundary = boundary, zctas = zctas)
  saveRDS(map_geometry, geometry_file)
}
map_data <- map_geometry$zctas %>%
  mutate(GEOID = as.character(ZCTA5CE20)) %>%
  left_join(clustered_zctas %>% select(GEOID, cluster), by = "GEOID") %>%
  left_join(acute_summary %>% select(cluster, P_beta_gt_0), by = "cluster")
if (length(setdiff(clustered_zctas$GEOID, map_data$GEOID))) {
  stop("Some classified ZCTAs are absent from the county-intersection geometry.")
}
assigned <- !is.na(map_data$cluster)
stopifnot(!anyNA(map_data$P_beta_gt_0[assigned]),
          all(map_data$P_beta_gt_0[assigned] == acute_summary$P_beta_gt_0[
            match(map_data$cluster[assigned], acute_summary$cluster)]))
geometry_coverage <- map_data %>%
  mutate(map_coverage = ifelse(sf::st_dimension(map_data) == 2L,
                               "polygon", "boundary_only")) %>%
  sf::st_drop_geometry() %>% select(GEOID, map_coverage)
zcta_probability_table <- clustered_zctas %>%
  select(GEOID, cluster, cluster_label) %>%
  left_join(acute_summary %>% select(cluster, P_beta_gt_0), by = "cluster") %>%
  left_join(geometry_coverage, by = "GEOID")
write.csv(zcta_probability_table,
          file.path(output_dir, "acute_heat_zcta_probabilities.csv"), row.names = FALSE)
map_polygons <- map_data[sf::st_dimension(map_data) == 2L, ]
unmapped <- setdiff(cluster_ids, unique(map_polygons$cluster))
if (length(unmapped)) stop("No map geometry for cluster(s): ", paste(unmapped, collapse = ", "))
cluster_boundaries <- map_polygons %>% filter(!is.na(P_beta_gt_0)) %>%
  group_by(cluster) %>% summarise(.groups = "drop")
if (model_variant == "linear_cluster") {
probability_map <- draw_probability_map(
  mainland_probability_geometry(map_geometry, clustered_zctas),
  data.frame(cluster = acute_summary$cluster, probability = acute_summary$P_beta_gt_0),
  expression(Pr(beta[acute*","*k] > 0 ~ "|" ~ Y)),
  if (model_variant %in% c("climate_zone_added", "climate_zone")) {
    paste0(model_label, " + shared and pooled climate-zone RW2 (nonlinear addition)\n",
      "Positive linear acute-heat component, adjusted for RW2 curves\n",
      "Mainland Los Angeles, May-October 2018")
  } else paste0(model_label, ": positive acute-heat coefficient\n",
    "Cluster-specific acute-heat slope\nMainland Los Angeles, May-October 2018"))
for (ext in c("png", "pdf")) {
  save_project_map(file.path(map_output_dir, paste0("acute_heat_probability_map_zcta.", ext)), probability_map)
}

saveRDS(list(probability_map = probability_map, rr_plot = rr_plot),
         file.path(output_dir, "acute_heat_ggplots.rds"))
if (interactive()) {
  print(probability_map)
  print(rr_plot)
}
} else {
  # The IRR script below also draws the matching full-response probability map
  # from the same INLA contrast marginals and aligns the two map panels.
  saveRDS(list(rr_plot = rr_plot), file.path(output_dir, "acute_heat_ggplots.rds"))
}
message("Finished. Figures and summary saved in: ", output_dir)
}
# Generate IRR maps
rm(acute_fit)
gc(verbose = FALSE)
irr_environment <- new.env(parent = environment())
irr_environment$irr_model_families <- model_families
irr_environment$model_variant <- model_variant
source(file.path(project_dir, "acute_heat_irr_maps.R"), local = irr_environment)
