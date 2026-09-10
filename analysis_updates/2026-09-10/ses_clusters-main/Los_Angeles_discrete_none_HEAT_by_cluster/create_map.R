library(tidyverse)
library(sf)
library(tigris)
library(scales)

state <- "CA"
county <- "Los Angeles"
file_name <- "Los Angeles"
folder_name <- "Los_Angeles_discrete_none_HEAT_by_cluster"
year <- 2020

options(tigris_use_cache = TRUE)


clustered_zctas <- readRDS("Los Angeles_ACS_zcta_clustered_for_ED.rds")

# Run from this output directory after hri_ed.R has generated the probabilities.
cluster_exceedance <- read_csv(
  paste0(file_name, "_poisson_cluster_AR_exceedance.csv"),
  show_col_types = FALSE) %>%
  transmute(cluster = as.integer(sub("^cluster", "", cluster)), P_AR_gt_0)

# Use the medians and vulnerability scores already exported by hri_ed.R.
profile_variables <- tribble(
  ~variable, ~display_label, ~theme,
  "pct_white", "White (NH)", "Race/ethnicity",
  "pct_black", "Black (NH)", "Race/ethnicity",
  "pct_asian", "Asian (NH)", "Race/ethnicity",
  "pct_hispanic", "Hispanic", "Race/ethnicity",
  "median_income", "Median income", "Socioeconomic",
  "poverty_rate", "Poverty", "Socioeconomic",
  "unemployment_rate", "Unemployment", "Socioeconomic",
  "education_low", "Low education", "Socioeconomic",
  "uninsured_rate", "Uninsured", "Socioeconomic",
  "pct_age17", "Age <18", "Household",
  "pct_age65", "Age 65+", "Household",
  "disability_rate", "Disability", "Household",
  "language_isolation", "Language isolation", "Household",
  "no_phone_rate", "No phone", "Household",
  "elderly_alone_rate", "Older adults alone", "Household",
  "renter_burden_rate", "Rent burden", "Housing",
  "housing_overcrowding", "Overcrowding", "Housing",
  "alt_housing_rate", "Alternative housing", "Housing",
  "old_housing_rate", "Pre-1980 housing", "Housing",
  "outdoor_worker_rate", "Outdoor-work proxy", "Work")

profile_summary <- read_csv(
  paste0(file_name, "_poisson_cluster_results.csv"),
  show_col_types = FALSE) %>%
  select(cluster, cluster_label, vulnerability_score, n_zcta = n.y,
         all_of(profile_variables$variable)) %>%
  arrange(vulnerability_score, cluster)
profile_reference <- profile_summary$cluster[1]
profile_row_labels <- with(profile_summary, setNames(
  paste0(cluster, ifelse(cluster == profile_reference, " [ref]", ""),
         "  (", n_zcta, ")"), cluster))

# Center and scale each variable across cluster medians, weighting clusters equally.
profile_cells <- profile_summary %>%
  pivot_longer(all_of(profile_variables$variable),
               names_to = "variable", values_to = "cluster_median") %>%
  left_join(profile_variables, by = "variable") %>%
  group_by(variable) %>%
  mutate(
    n_clusters_observed = sum(!is.na(cluster_median)),
    median_center = mean(cluster_median, na.rm = TRUE),
    median_sd = sd(cluster_median, na.rm = TRUE),
    standardized_median = if_else(
      median_sd > 0, (cluster_median - median_center) / median_sd, NA_real_)) %>%
  ungroup()

# Ward ordering changes only the display, not the supplied cluster assignments.
profile_matrix <- profile_cells %>%
  select(cluster, variable, standardized_median) %>%
  pivot_wider(names_from = variable, values_from = standardized_median) %>%
  column_to_rownames("cluster") %>%
  as.matrix()
profile_row_tree <- hclust(dist(profile_matrix), method = "ward.D2")
profile_row_order <- as.integer(rownames(profile_matrix)[profile_row_tree$order])
profile_column_order <- unlist(lapply(unique(profile_variables$theme), function(group) {
  variables <- profile_variables$variable[profile_variables$theme == group]
  if (length(variables) == 1L) return(variables)
  tree <- hclust(dist(t(profile_matrix[, variables, drop = FALSE])), method = "ward.D2")
  variables[tree$order]
}), use.names = FALSE)
profile_cells <- profile_cells %>%
  mutate(
    row_order = match(cluster, profile_row_order),
    column_order = match(variable, profile_column_order)) %>%
  arrange(row_order, column_order)
write_csv(profile_cells, "cluster_profiles_heatmap.csv")

profile_plot_data <- profile_cells %>%
  mutate(
    cluster = factor(cluster, levels = rev(profile_row_order)),
    variable = factor(variable, levels = profile_column_order),
    theme = factor(theme, levels = unique(profile_variables$theme)))
profile_limit <- max(abs(profile_cells$standardized_median), na.rm = TRUE)
profile_heatmap <- ggplot(profile_plot_data, aes(variable, cluster)) +
  geom_tile(aes(fill = standardized_median), color = "white", linewidth = 0.2) +
  facet_grid(. ~ theme, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = setNames(
    profile_variables$display_label, profile_variables$variable), expand = c(0, 0)) +
  scale_y_discrete(labels = profile_row_labels, expand = c(0, 0)) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-profile_limit, profile_limit), na.value = "grey70",
    name = "Standardized cluster median") +
  labs(x = NULL, y = "Cluster (number of ZCTAs)") +
  guides(fill = guide_colorbar(
    title.position = "top", title.hjust = 0.5,
    barwidth = unit(6, "cm"), barheight = unit(0.35, "cm"))) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.x = unit(0.12, "cm"),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10),
    axis.ticks = element_blank(),
    legend.position = "bottom")
ggsave("cluster_profiles_heatmap.pdf", profile_heatmap,
       width = 9, height = 5.5, bg = "white")
ggsave("cluster_profiles_heatmap.png", profile_heatmap,
       width = 9, height = 5.5, dpi = 300, bg = "white")

# =====================================================================
# GET ZCTA GEOMETRIES
# =====================================================================

boundary <- counties(state = state, cb = TRUE) %>%
  filter(NAME == county) %>%
  st_make_valid()

# Display only the largest connected county polygon (the mainland).
# This clips the map, not the data used to estimate cluster probabilities.
boundary <- st_cast(boundary, "POLYGON")
boundary <- boundary[which.max(st_area(boundary)), ]

zcta_shapes <- zctas(cb = TRUE, year = year) %>%
  st_transform(st_crs(boundary)) %>%
  st_make_valid()

zcta_cropped <- st_intersection(zcta_shapes, boundary) %>%
  st_make_valid() %>%
  left_join(clustered_zctas %>%
      st_drop_geometry() %>%
      select(GEOID, cluster, cluster_label),
    by = c("ZCTA5CE20" = "GEOID")) %>%
  left_join(cluster_exceedance, by = "cluster")

# dissolve zctas into clusters for boundaries
cluster_shapes <- zcta_cropped %>%
  filter(!is.na(cluster)) %>%
  group_by(cluster, cluster_label) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

# =====================================================================
# MAP
# =====================================================================
map <- ggplot() +
  geom_sf(data = boundary, fill = "grey70", color = "white", linewidth = 0.1) +
  geom_sf(data = zcta_cropped, aes(fill = P_AR_gt_0), color = NA) +
  geom_sf(data = cluster_shapes, fill = NA, color = "black", linewidth = 0.25) +
  coord_sf(expand = FALSE) + 
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0.5,
    na.value = "grey70",
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = percent_format(accuracy = 1), name = expression(Pr(AR[c] > 0 ~ "|" ~ D))) + labs(
      title = paste("Cluster AR Exceedance Across Mainland", file_name),
      subtitle = "Posterior probability of positive cluster-average excess rate") +
  theme_minimal() + theme(
    legend.position = "right",
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(0.35, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    plot.title = element_text(size = 12, margin = margin(b = 8)),
    plot.subtitle = element_text(size = 9, margin = margin(b = 15)),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7))

ggsave(filename = "AR_gt_0_gradient.png", plot = map,
  width = 8, height = 6, dpi = 300, bg = "white")
