# Map acute-heat IRRs: linear slopes or full nonlinear heat contrasts.
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
input_dir <- file.path(project_dir, "Los_Angeles_discrete_none")
model_prefix <- "Los_Angeles_discrete_none_HEAT_by_cluster"
if (!exists("model_variant", inherits = FALSE)) model_variant <- "climate_zone"
results_root <- if (model_variant %in% c("climate_zone_added", "climate_zone")) {
  file.path(input_dir, model_variant)
} else file.path(input_dir, "linear_cluster")
model_dirs <- c(Poisson = "acute_heat_poisson", "Negative binomial" = "acute_heat_nb")
if (exists("irr_model_families", inherits = FALSE)) {
  model_dirs <- model_dirs[c("poisson", "nbinomial") %in% irr_model_families]
}

read_irr_map_inputs <- function(model_dir) {
  if (model_variant %in% c("climate_zone_added", "climate_zone")) {
    assignments <- readRDS(file.path(input_dir, "Los Angeles_ACS_zcta_clustered_for_ED.rds"))
    assignments <- sf::st_drop_geometry(assignments) %>%
      transmute(GEOID = as.character(GEOID), cluster = as.character(cluster))
    geometry <- readRDS(file.path(input_dir, "linear_cluster", basename(model_dir), "los_angeles_zcta_2020.rds"))
    return(read_full_heat_irr_inputs(model_dir, assignments, geometry))
  }
  estimates <- read.csv(file.path(model_dir, "acute_heat_cluster_summary.csv"),
                        colClasses = c(cluster = "character"))
  assignments <- read.csv(file.path(model_dir, "acute_heat_zcta_probabilities.csv"),
                          colClasses = c(GEOID = "character", cluster = "character"))
  geometry <- readRDS(file.path(input_dir, "linear_cluster", basename(model_dir), "los_angeles_zcta_2020.rds"))
  if (model_variant %in% c("climate_zone_added", "climate_zone")) {
    return(read_full_heat_irr_inputs(model_dir, assignments, geometry))
  }
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

# Full nonlinear contrast: one global heat SD above each ZIP's mean versus
# its mean (acute_heat 1 versus 0), with all other covariates held fixed.
# Population exposure cancels in this within-ZIP rate ratio.
irr_heat_from <- 0
irr_heat_to <- 1
read_full_heat_irr_inputs <- function(model_dir, assignments, geometry) {
  suffix <- if (basename(model_dir) == "acute_heat_poisson") "poisson" else "nb"
  model_file <- file.path(results_root, paste0(model_prefix,
    "_suppressed_INLA_", suffix, "_model.rds"))
  model <- readRDS(model_file)
  meta <- model$nonlinear_heat
  has_fixed_heat <- "acute_heat" %in% model$names.fixed
  info <- model$full_heat_irr
  if (is.null(info) || !all(info$combinations$lincomb_id %in%
      names(model$marginals.lincomb.derived))) {
    stop("The saved fit lacks full-response IRR linear-combination marginals. ",
      "Request them with fit_heat_model(..., nonlinear_heat = nonlinear_heat) first; ",
      "no IRR sampling is performed.")
  }
  stopifnot(isTRUE(meta$pooled), info$heat_from == irr_heat_from,
    info$heat_to == irr_heat_to)
  prepared <- load_heat_data(input_dir)
  data <- prepared$daily_ed
  stopifnot(nrow(data) == length(meta$zone),
    nrow(data) == nrow(model$summary.linear.predictor))
  zip_zones <- distinct(data.frame(GEOID = as.character(data$zip),
    cluster = as.character(data$cluster), climate_zone = meta$zone))
  stopifnot(!anyDuplicated(zip_zones$GEOID))
  combinations <- info$combinations
  key <- function(cluster, zone) paste(cluster, ifelse(is.na(zone), "missing", zone), sep = ":")
  combo_index <- match(key(zip_zones$cluster, zip_zones$climate_zone),
    key(combinations$cluster, combinations$climate_zone))
  stopifnot(!anyNA(combo_index))
  log_quantiles <- t(vapply(combinations$lincomb_id, function(id)
    INLA::inla.qmarginal(c(0.025, 0.5, 0.975), model$marginals.lincomb.derived[[id]]),
    numeric(3)))
  probabilities <- vapply(combinations$lincomb_id, function(id)
    1 - INLA::inla.pmarginal(0, model$marginals.lincomb.derived[[id]]), numeric(1))
  estimates <- zip_zones %>% mutate(
    RR_low = exp(log_quantiles[combo_index, 1L]),
    RR_median = exp(log_quantiles[combo_index, 2L]),
    RR_high = exp(log_quantiles[combo_index, 3L]),
    P_IRR_gt_1 = probabilities[combo_index],
    acute_heat_from = irr_heat_from, acute_heat_to = irr_heat_to,
    heat_sd_original_units = prepared$heat_sd)
  stopifnot(all(is.finite(estimates$RR_high)), all(estimates$RR_low > 0),
    all(estimates$P_IRR_gt_1 >= 0 & estimates$P_IRR_gt_1 <= 1))
  rm(model)
  gc(verbose = FALSE)
  stopifnot(!anyDuplicated(estimates$GEOID), setequal(estimates$GEOID, assignments$GEOID),
    all(estimates$cluster == assignments$cluster[match(estimates$GEOID, assignments$GEOID)]))
  estimates <- left_join(estimates, prepared$cluster_labels, by = "cluster")
  write.csv(estimates, file.path(model_dir, "acute_heat_full_response_irr_zcta.csv"), row.names = FALSE)
  probability_table <- assignments %>%
    select(-any_of(c("P_IRR_gt_1", "acute_heat_from", "acute_heat_to"))) %>%
    left_join(select(estimates, GEOID, P_IRR_gt_1, acute_heat_from, acute_heat_to), by = "GEOID")
  write.csv(probability_table, file.path(model_dir, "acute_heat_zcta_probabilities.csv"), row.names = FALSE)
  stopifnot(!length(setdiff(estimates$GEOID, as.character(geometry$zctas$ZCTA5CE20))))
  # Use the probability map's mainland extent for both full-response maps.
  mainland <- mainland_probability_geometry(geometry, assignments)
  polygons <- mainland$polygons %>%
    left_join(select(estimates, -cluster), by = "GEOID")
  stopifnot(!anyNA(polygons$RR_median))
  list(estimates = estimates, polygons = polygons, boundary = mainland$boundary,
    hatching = mainland$hatching, full_response = TRUE, has_fixed_heat = has_fixed_heat)
}

# Export additional comparison plots from the saved one-SD INLA contrasts.
export_grouped_heat_irrs <- function(model_dir, variant, plot_dir, input_dir,
                                     model_name = "Poisson") {
  stopifnot(variant %in% c("climate_zone", "climate_zone_added"))
  estimates <- read.csv(file.path(model_dir, "acute_heat_full_response_irr_zcta.csv"),
    colClasses = c(GEOID = "character", cluster = "character", climate_zone = "character"))
  stopifnot(all(estimates$acute_heat_from == 0), all(estimates$acute_heat_to == 1),
    !anyDuplicated(estimates$GEOID), all(is.finite(estimates$RR_high)),
    all(estimates$RR_low > 0), all(estimates$RR_low <= estimates$RR_median),
    all(estimates$RR_median <= estimates$RR_high))
  keys <- if (variant == "climate_zone") "climate_zone" else c("cluster", "climate_zone")
  # Identical model contrasts are repeated across ZIPs; retain one, never average intervals.
  consistency <- estimates %>% group_by(across(all_of(keys))) %>%
    summarise(across(c(RR_low, RR_median, RR_high), ~ diff(range(.x))), .groups = "drop")
  stopifnot(all(as.matrix(select(consistency, RR_low, RR_median, RR_high)) < 1e-8))
  table <- estimates %>% group_by(across(all_of(keys))) %>%
    summarise(n_zctas = n_distinct(GEOID),
      across(c(RR_low, RR_median, RR_high, P_IRR_gt_1), first), .groups = "drop")
  zones <- sort(unique(as.numeric(na.omit(estimates$climate_zone))))
  zone_levels <- c(paste("Climate zone", zones),
    if (anyNA(estimates$climate_zone)) "Climate zone unavailable")
  table$zone_label <- factor(ifelse(is.na(table$climate_zone), "Climate zone unavailable",
    paste("Climate zone", table$climate_zone)), levels = zone_levels)
  caption <- paste0("One county-wide heat SD above each ZIP's overall May-October mean versus that mean.\n",
    "Points are posterior medians; bars are 95% credible intervals. Dashed line: IRR = 1.")
  if (variant == "climate_zone_added") {
    assignments <- sf::st_drop_geometry(readRDS(file.path(input_dir,
      "Los Angeles_ACS_zcta_clustered_for_ED.rds")))
    vulnerability <- compute_vulnerability_score(assignments) %>%
      mutate(cluster = as.character(cluster)) %>% group_by(cluster) %>%
      summarise(vulnerability_score = mean(vulnerability_score),
        cluster_label = first(cluster_label), .groups = "drop") %>% arrange(vulnerability_score)
    table <- left_join(table, vulnerability, by = "cluster")
    stopifnot(!anyNA(table$vulnerability_score))
    table$plot_row <- factor(table$cluster, levels = vulnerability$cluster)
    row_labels <- function(x) paste("Cluster", x)
    caption <- paste0(caption, "\nClusters run from least vulnerable (bottom) to most vulnerable (top); blank rows have no modeled ZIP.")
  } else {
    table$plot_row <- factor(table$zone_label, levels = rev(zone_levels))
    row_labels <- identity
  }
  plot <- ggplot(table, aes(y = plot_row)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
    geom_segment(aes(x = RR_low, xend = RR_high, yend = plot_row),
      linewidth = 0.55, colour = "#2166AC") +
    geom_point(aes(x = RR_median), size = 2.2, colour = "#2166AC") +
    scale_x_log10(breaks = scales::breaks_log(n = 5), labels = scales::label_number()) +
    scale_y_discrete(labels = row_labels, drop = FALSE) +
    labs(title = if (variant == "climate_zone") "Acute-heat IRRs by climate zone" else
      "Acute-heat IRRs by vulnerability cluster and climate zone",
      subtitle = paste0(model_name, " model | ", variant),
      x = "Incidence rate ratio for a 1-SD heat increase (log scale)",
      y = if (variant == "climate_zone_added") "Vulnerability cluster" else NULL,
      caption = caption) + theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
      plot.caption = element_text(hjust = 0, size = 9), plot.title.position = "plot",
      plot.margin = margin(12, 18, 12, 12))
  if (variant == "climate_zone_added") plot <- plot + facet_wrap(~zone_label, ncol = 3)
  basename <- if (variant == "climate_zone") "acute_heat_irr_by_climate_zone" else
    "acute_heat_irr_by_cluster_climate_zone"
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(select(table, -plot_row), file.path(model_dir, paste0(basename, ".csv")), row.names = FALSE)
  saveRDS(plot, file.path(model_dir, paste0(basename, ".rds")))
  for (extension in c("png", "pdf")) {
    temporary <- tempfile(pattern = "irr_plot_", tmpdir = plot_dir, fileext = paste0(".", extension))
    on.exit(unlink(temporary), add = TRUE)
    ggsave(temporary, plot,
      width = if (variant == "climate_zone") 10 else 12,
      height = if (variant == "climate_zone") 6 else 9.5, units = "in", dpi = 300,
      bg = "white", device = if (extension == "pdf") grDevices::cairo_pdf else NULL)
    stopifnot(file.exists(temporary), file.info(temporary)$size > 0,
      file.rename(temporary, file.path(plot_dir, paste0(basename, ".", extension))))
    gc(verbose = FALSE)
  }
  invisible(table)
}

map_inputs <- lapply(model_dirs, function(folder) read_irr_map_inputs(file.path(results_root, folder)))

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
  if (isTRUE(x$full_response)) {
    boundaries <- x$polygons %>% group_by(cluster) %>% summarise(.groups = "drop")
    full_irr_limits <- c(min(0.9999, estimates$RR_median), max(1.0001, estimates$RR_median))
    full_irr_breaks <- seq(ceiling(full_irr_limits[1] / 0.2),
      floor(full_irr_limits[2] / 0.2)) * 0.2
    return(ggplot() +
      geom_sf(data = x$boundary, fill = "grey92", colour = "grey60", linewidth = 0.2) +
      missing_result_hatch_layer(x) +
      geom_sf(data = x$polygons, aes(fill = RR_median), colour = "white", linewidth = 0.10) +
      geom_sf(data = boundaries, fill = NA, colour = "black", linewidth = 0.30) +
      # The observed minimum is blue, IRR 1 is white, the maximum is red.
      scale_fill_gradientn(colours = c("#2166AC", "white", "#B2182B"),
        values = scales::rescale(c(full_irr_limits[1], 1, full_irr_limits[2]), from = full_irr_limits),
        limits = full_irr_limits,
        breaks = full_irr_breaks, labels = function(x) sprintf("%.1f", x),
        name = "Acute-heat IRR\n(1 = no change)",
        guide = guide_colourbar(barheight = grid::unit(65, "mm"))) +
      coord_sf(datum = NA, expand = FALSE) +
      labs(title = "Los Angeles acute-heat incidence rate ratios",
        subtitle = paste0(model_name, " model | Posterior median IRRs for a 1-SD heat increase\n",
          "from each ZIP's average"),
        caption = "SD is calculated county-wide. ZIPs with missing climate zones use the shared RW2 curve only.") +
      theme_void(base_size = 11) +
      theme(legend.position = "right", plot.title.position = "plot",
        plot.caption = element_text(size = 8, hjust = 0),
        plot.margin = margin(10, 12, 10, 12)))
  }
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
         subtitle = if (model_variant %in% c("climate_zone_added", "climate_zone")) {
           paste0(model_name, " + shared and pooled climate-zone RW2 (nonlinear addition)\n",
             "Linear acute-heat component: posterior medians and 95% credible intervals")
         } else paste0(model_name, " model | Posterior medians and 95% credible intervals")) +
    
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
  map_dir <- file.path(project_dir, "maps", model_variant, model_dirs[[model_name]])
  dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
  
  plots <- list(acute_heat_irr_map = irr_maps[[model_name]])
  if (isTRUE(map_inputs[[model_name]]$full_response)) {
    x <- map_inputs[[model_name]]
    geometry <- list(boundary = x$boundary, hatching = x$hatching, polygons = x$polygons %>%
      select(-P_IRR_gt_1), clusters = x$polygons %>% group_by(cluster) %>% summarise(.groups = "drop"))
    probability_map <- draw_probability_map(geometry,
      transmute(x$estimates, GEOID, probability = P_IRR_gt_1),
      expression(Pr(IRR["1SD"] > 1 ~ "|" ~ Y)),
      paste0(model_name, " model | Posterior probability of an increased rate\n",
        "for a 1-SD heat increase from each ZIP's average"), by = "GEOID") +
      labs(caption = "SD is calculated county-wide. ZIPs with missing climate zones use the shared RW2 curve only.")
    # Express the posterior probability as a percentage.
    probability_scale <- probability_map$scales$get_scales("fill")
    probability_scale$labels <- scales::label_percent(accuracy = 1)
    plots$acute_heat_probability_map_zcta <- probability_map
    bundle_file <- file.path(results_root, model_dirs[[model_name]], "acute_heat_ggplots.rds")
    bundle <- if (file.exists(bundle_file)) readRDS(bundle_file) else list()
    if (identical(x$has_fixed_heat, FALSE)) bundle$rr_plot <- NULL
    bundle$probability_map <- probability_map
    saveRDS(bundle, bundle_file)
  }
  # Match canvas dimensions AND panel dimensions, including legend/margin space.
  grobs <- lapply(plots, ggplotGrob)
  widths <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$widths))
  heights <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$heights))
  for (i in seq_along(grobs)) {
    grobs[[i]]$widths <- widths
    grobs[[i]]$heights <- heights
    for (extension in c("png", "pdf")) {
      save_project_map(file.path(map_dir, paste0(names(plots)[i], ".", extension)), grobs[[i]])
    }
  }
  
  saveRDS(irr_maps[[model_name]], file.path(results_root, model_dirs[[model_name]], "acute_heat_irr_map.rds"))
  if (isTRUE(map_inputs[[model_name]]$full_response)) {
    export_grouped_heat_irrs(file.path(results_root, model_dirs[[model_name]]),
      model_variant, map_dir, input_dir, model_name)
  }
  if (interactive()) print(irr_maps[[model_name]])
}

message("Finished. IRR maps saved in: ", paste(unname(model_dirs), collapse = ", "))
