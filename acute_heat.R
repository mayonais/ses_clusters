# Choose "poisson", "nbinomial", or "both"; then source this whole script.
model_family <- "both"
# RR_k = exp(beta_acute,k) per one global heat SD, holding chronic heat fixed.
# P(beta_acute,k > 0 | Y) = P(RR_k > 1 | Y).
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
input_dir <- file.path(project_dir, "Los_Angeles_discrete_none_HEAT_by_cluster")
model_families <- switch(model_family, both = c("poisson", "nbinomial"),
                         poisson = "poisson", nbinomial = "nbinomial",
                         stop('model_family must be "poisson", "nbinomial", or "both".'))
model_cache_dir <- file.path(input_dir, "model_cache")
# FALSE reuses a compatible fit and refits automatically if its signature changes.
# Set TRUE only when you want to force a fresh fit, independently of plot changes.
force_refit <- FALSE
geometry_year <- 2020L
source(file.path(project_dir, "analysis_functions.R"), local = TRUE)

# ---- Shared analysis --------------------------------------------------------
prepared <- load_heat_data(input_dir, study_name)
clustered_zctas <- prepared$clustered_zctas
cluster_ids <- levels(prepared$daily_ed$cluster)
for (family in model_families) {
model_label <- if (family == "poisson") "Poisson" else "Negative-binomial"
output_dir <- file.path(input_dir, if (family == "poisson") "acute_heat_poisson" else "acute_heat_nb")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
acute_fit <- fit_heat_model(prepared, model_cache_dir, family = family,
                            force_refit = force_refit)
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
save_plot(rr_plot, "acute_heat_relative_risk")

# ---- Probability choropleth: same ZCTA vintage as the existing map ------------
geometry_file <- file.path(output_dir, paste0("los_angeles_zcta_", geometry_year, ".rds"))
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
probability_map <- ggplot() +
  geom_sf(data = map_geometry$boundary, fill = "grey92", colour = "grey60", linewidth = 0.2) +
  geom_sf(data = filter(map_polygons, !is.na(P_beta_gt_0)), aes(fill = P_beta_gt_0),
          colour = "white", linewidth = 0.10) +
  geom_sf(data = cluster_boundaries, fill = NA, colour = "grey25", linewidth = 0.25) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0.5, limits = c(0, 1), breaks = seq(0, 1, 0.2),
                       labels = function(x) sprintf("%.1f", x), na.value = "grey92",
                       name = expression(P(beta[acute*","*k] > 0 ~ "|" ~ Y))) +
  coord_sf(datum = NA, expand = FALSE) +
  labs(title = paste(study_name, "acute-heat exceedance probabilities"),
       subtitle = paste0(model_label, " posterior probability that the acute-heat coefficient exceeds 0")) +
  theme_void(base_size = 11) + theme(legend.position = "right")
save_plot(probability_map, "acute_heat_probability_map_zcta", width = 10, height = 8)

saveRDS(list(probability_map = probability_map, rr_plot = rr_plot),
         file.path(output_dir, "acute_heat_ggplots.rds"))
print(probability_map)
print(rr_plot)
message("Finished. Figures and summary saved in: ", output_dir)
}
# Generate IRR maps
irr_environment <- new.env(parent = environment())
irr_environment$irr_model_families <- model_families
source(file.path(project_dir, "acute_heat_irr_maps.R"), local = irr_environment)
