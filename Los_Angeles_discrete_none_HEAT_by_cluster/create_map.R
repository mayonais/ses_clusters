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

cluster_exceedance <- tibble(cluster = 1:14,
  P_AR_gt_0 = c(0.505, 0.514, 0.514, 0.578, 0.726, 0.653, 0.550,
    0.784, 0.525, 0.494, 0.519, 0.495, 0.512, 0.642))

# =====================================================================
# GET ZCTA GEOMETRIES
# =====================================================================

boundary <- counties(state = state, cb = TRUE) %>%
  filter(NAME == county) %>%
  st_make_valid()

zcta_shapes <- zctas(cb = TRUE, year = year) %>%
  st_transform(st_crs(boundary)) %>%
  st_make_valid()

zcta_cropped <- st_intersection(zcta_shapes, boundary) %>%
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
  geom_sf(data = boundary, fill = "grey90", color = "white", linewidth = 0.1) +
  geom_sf(data = zcta_cropped, aes(fill = P_AR_gt_0), color = NA) +
  geom_sf(data = cluster_shapes, fill = NA, color = "black", linewidth = 0.25) +
  coord_sf(expand = FALSE) + scale_fill_gradient(
    low = "#E8F5E9", high = "#006D2C", limits = c(0.49, 0.79),
    breaks = c(0.50, 0.60, 0.70, 0.78),
    labels = percent_format(accuracy = 1), name = "P(AR > 0)") + labs(
    title = paste("Heat Exceedance Probabilities Across", file_name, "SES Clusters"),) +
  theme_minimal() + theme(aspect.ratio = 1.6,
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
  width = 8, height = 6, dpi = 300)