# Map the saved acute-heat posterior median IRRs
library(dplyr)
library(ggplot2)
library(sf)

script_files <- as.character(unlist(lapply(sys.frames(), function(frame) frame$ofile)))
script_files <- script_files[basename(script_files) == "acute_heat_irr_maps.R"]
command_file <- grep("^--file=", commandArgs(), value = TRUE)
project_dir <- if (length(script_files)) {
  dirname(normalizePath(tail(script_files, 1), mustWork = TRUE))
} else if (length(command_file) && sub("^--file=", "", command_file[1]) != "-") {
  dirname(normalizePath(sub("^--file=", "", command_file[1]), mustWork = TRUE))
} else {
  getwd()
}
source(file.path(project_dir, "shared_functions.R"))
input_dir <- file.path(project_dir, "Los_Angeles_discrete_none_HEAT_by_cluster")
model_dirs <- c(Poisson = "acute_heat_poisson", "Negative binomial" = "acute_heat_nb")
if (exists("irr_model_families", inherits = FALSE)) {
  model_dirs <- model_dirs[c("poisson", "nbinomial") %in% irr_model_families]
}

read_irr_map_inputs <- function(model_dir) {
  estimates <- read.csv(file.path(model_dir, "acute_heat_cluster_summary.csv"),
                        colClasses = c(cluster = "character"))
  assignments <- read.csv(file.path(model_dir, "acute_heat_zcta_probabilities.csv"),
                          colClasses = c(GEOID = "character", cluster = "character"))
  geometry <- readRDS(file.path(model_dir, "los_angeles_zcta_2020.rds"))
  stopifnot(all(c("cluster", "cluster_label", "RR_median", "RR_low", "RR_high") %in%
                  names(estimates)),
            !anyDuplicated(estimates$cluster), !anyDuplicated(assignments$GEOID),
            !anyNA(estimates[c("cluster", "cluster_label")]),
            !anyNA(assignments[c("GEOID", "cluster")]),
            all(is.finite(as.matrix(estimates[c("RR_median", "RR_low", "RR_high")]))),
            all(estimates$RR_low > 0),
            all(estimates$RR_low <= estimates$RR_median),
            all(estimates$RR_median <= estimates$RR_high),
            setequal(estimates$cluster, assignments$cluster))
  polygons <- geometry$zctas %>%
    mutate(GEOID = as.character(ZCTA5CE20)) %>%
    left_join(select(assignments, GEOID, cluster), by = "GEOID") %>%
    left_join(select(estimates, cluster, RR_median), by = "cluster")
  stopifnot(!length(setdiff(assignments$GEOID, polygons$GEOID)),
            !anyNA(polygons$RR_median[!is.na(polygons$cluster)]))
  polygons <- polygons[sf::st_dimension(polygons) == 2L & !is.na(polygons$cluster), ]
  stopifnot(setequal(estimates$cluster, polygons$cluster))
  estimates <- estimates %>% arrange(RR_median, as.numeric(cluster))
  list(estimates = estimates, polygons = polygons, boundary = geometry$boundary)
}

map_inputs <- lapply(model_dirs, function(folder) read_irr_map_inputs(file.path(input_dir, folder)))

green_palette <- c("#F7FCF5", "#E2F2DC", "#B9DFAE", "#67A96B", "#1E5A32")
irr_scale <- function(limits, midpoint) {
  scale_fill_gradientn(
    colours = green_palette,
    values = scales::rescale(
      c(limits[1], midpoint, limits[2]),
      from = limits
    ),
    limits = limits,
    guide = "none"
  )
}

irr_maps <- lapply(names(map_inputs), function(model_name) {
  x <- map_inputs[[model_name]]
  estimates <- x$estimates
  cluster_order <- estimates$cluster
  
  irr_limits <- range(estimates$RR_median)
  if (diff(irr_limits) == 0) irr_limits <- irr_limits + c(-0.1, 0.1)
  
  irr_mid <- median(estimates$RR_median)
  
  legend_scale <- irr_scale(irr_limits, irr_mid)
  legend_scale$train(estimates$RR_median)
  swatches <- setNames(legend_scale$map(estimates$RR_median), cluster_order)
  
  legend_labels <- setNames(sprintf("%s\nIRR %.4f  (95%% CrI %.4f to %.4f)",
                                    estimates$cluster_label, estimates$RR_median,
                                    estimates$RR_low, estimates$RR_high), cluster_order)
  
  boundaries <- x$polygons %>% group_by(cluster) %>% summarise(.groups = "drop")
  bbox <- sf::st_bbox(x$boundary)
  
  ggplot() +
    geom_sf(data = x$boundary, fill = "grey92", colour = "grey60", linewidth = 0.2) +
    geom_sf(data = x$polygons, aes(fill = RR_median),
            colour = "white", linewidth = 0.10) +
    geom_sf(data = boundaries, fill = NA, colour = "black", linewidth = 0.30) +
    
    geom_point(data = estimates, aes(colour = cluster),
               x = mean(bbox[c("xmin", "xmax")]), y = mean(bbox[c("ymin", "ymax")]),
               size = 0, alpha = 0, inherit.aes = FALSE, show.legend = TRUE) +
    
    irr_scale(irr_limits, irr_mid) +
    
    scale_colour_manual(values = swatches, breaks = cluster_order,
                        labels = legend_labels, name = "Cluster profile and 95% credible interval",
                        guide = guide_legend(order = 1, ncol = 1,
                                             override.aes = list(alpha = 1, shape = 15, size = 4))) +
    
    coord_sf(datum = NA, expand = FALSE) +
    
    labs(title = "Los Angeles acute-heat incidence rate ratios",
         subtitle = paste0(model_name, " model | Posterior medians and 95% credible intervals")) +
    
    theme_void(base_size = 11) +
    
    theme(legend.position = "right", legend.box = "vertical",
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 9, lineheight = 1.05),
          legend.key.height = grid::unit(6, "mm"),
          legend.key.width = grid::unit(5, "mm"),
          legend.key.spacing.y = grid::unit(1.5, "mm"),
          legend.spacing.y = grid::unit(1.5, "mm"),
          plot.title.position = "plot",
          plot.margin = margin(10, 12, 10, 12))
})
names(irr_maps) <- names(map_inputs)

for (model_name in names(irr_maps)) {
  map_dir <- file.path(project_dir, "maps", model_dirs[[model_name]])
  dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (extension in c("png", "pdf")) {
    save_project_map(file.path(map_dir, paste0("acute_heat_irr_map.", extension)),
                     irr_maps[[model_name]])
  }
  
  saveRDS(irr_maps[[model_name]], file.path(input_dir, model_dirs[[model_name]], "acute_heat_irr_map.rds"))
  print(irr_maps[[model_name]])
}

message("Finished. IRR maps saved in: ", paste(unname(model_dirs), collapse = ", "))
