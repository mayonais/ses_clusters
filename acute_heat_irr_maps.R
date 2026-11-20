# Run: Rscript acute_heat_irr_maps.R --nb-pc
load_nb_map_inputs <- function(results_dir) {
  map_style <- readRDS("Los_Angeles_discrete_none/heat_map_style.rds")
  stopifnot(identical(map_style$schema_version, 1L))
  reference_irr <- map_style$irr
  reference_probability <- map_style$probability
  reference_ar <- map_style$excess_rate
  if (!file.exists(file.path(results_dir, "summaries.rds")))
    stop("Fit the revised NB model with hri_ed_analysis.R before exporting maps.")
  prepared <- load_heat_data("Los_Angeles_discrete_none")
  results <- readRDS(file.path(results_dir, "summaries.rds"))
  has_linear_heat <- "acute_heat" %in% results$sampling_info$names.fixed
  interaction_count <- sum(grepl("cluster.*:acute_heat|acute_heat:cluster", results$sampling_info$names.fixed))
  stopifnot(interaction_count == if (has_linear_heat) 13L else 0L)
  heat_sd_f <- prepared$heat_sd
  heat_increment <- sprintf("1 SD = %.2f\u00b0F", heat_sd_f)
  heat_structure <- if (has_linear_heat) "SES heat slopes + shared/zone curves" else "Shared/zone heat curves"
  model_label <- paste("Negative binomial", heat_structure, "PC hyperpriors", sep = " | ")
  reference <- prepared$reference_cluster
  assignments <- st_drop_geometry(prepared$clustered_zctas) %>%
    transmute(GEOID = as.character(GEOID), cluster = as.character(cluster), cluster_label)
  zones <- read.csv("datasets/cz_zip.csv", colClasses = "character")
  zip_info <- prepared$daily_ed %>% distinct(zip, cluster) %>%
    transmute(GEOID = as.character(zip), cluster = as.character(cluster)) %>%
    left_join(zones %>% rename(GEOID = zip, climate_zone = cz), by = "GEOID")
  key <- function(cluster, zone) paste(cluster, ifelse(is.na(zone), "missing", zone), sep = ":")
  idx <- match(key(zip_info$cluster, zip_info$climate_zone),
               key(results$heat$cluster, results$heat$climate_zone))
  stopifnot(!anyNA(idx), !anyDuplicated(zip_info$GEOID))
  zip_info <- bind_cols(zip_info,
    results$heat[idx, c("IRR_median", "IRR_lower", "IRR_upper", "P_IRR_gt_1")]) %>%
    mutate(heat_sd_f = heat_sd_f)
  baseline <- results$clusters %>% mutate(cluster = as.character(cluster)) %>%
    select(-family) %>% left_join(distinct(select(assignments, cluster, cluster_label)), by = "cluster")

  list(prepared = prepared, results = results, reference_irr = reference_irr,
    reference_probability = reference_probability, reference_ar = reference_ar,
    has_linear_heat = has_linear_heat, heat_sd_f = heat_sd_f, heat_increment = heat_increment,
    heat_structure = heat_structure, model_label = model_label, reference = reference,
    assignments = assignments, zip_info = zip_info, baseline = baseline, results_dir = results_dir)
}

prepare_nb_map_geometry <- function(assignments) {
  geometry <- readRDS("Los_Angeles_discrete_none/los_angeles_zcta_2020.rds")
  mainland <- mainland_probability_geometry(geometry, assignments)
  boundary <- mainland$boundary
  polygons <- mainland$polygons
  stopifnot(setequal(assignments$GEOID, as.character(geometry$zctas$ZCTA5CE20)[
    as.character(geometry$zctas$ZCTA5CE20) %in% assignments$GEOID]))
  boundary_only <- setdiff(assignments$GEOID, polygons$GEOID)
  no_area <- as.character(geometry$zctas$ZCTA5CE20[st_dimension(geometry$zctas) < 2L])
  not_area <- data.frame(GEOID = boundary_only, reason = ifelse(boundary_only %in% no_area,
    "Boundary-only county intersection", "Outside mainland display extent"))
  clusters <- polygons %>% group_by(cluster) %>% summarise(.groups = "drop")
  # Label every disconnected polygon, so a cluster number can appear repeatedly.
  label_parts <- suppressWarnings(st_cast(st_cast(
    st_collection_extract(clusters, "POLYGON"), "MULTIPOLYGON"), "POLYGON"))
  # Maximal-inscribed-circle centers keep labels away from region boundaries.
  circles <- st_inscribed_circle(st_geometry(label_parts), nQuadSegs = 0)
  centers <- lapply(circles, function(x) st_point(unname(x[1, ])))
  radius <- vapply(circles, function(x) sqrt(sum((x[1, ] - x[2, ])^2)), numeric(1))
  label_points <- st_sf(st_drop_geometry(label_parts), geometry = st_sfc(centers, crs = st_crs(label_parts)))
  xy <- st_coordinates(label_points)
  label_data <- st_drop_geometry(label_points) %>% mutate(x = xy[, 1], y = xy[, 2],
    small_region = radius < 1500)
  units_per_mm <- unname(diff(st_bbox(boundary)[c("xmin", "xmax")])) / (6.5 * 25.4)
  label_data$label_size <- pmin(4.5, pmax(1.6,
    2 * radius / units_per_mm / pmax(1, nchar(label_data$cluster) * 0.65)))
  stopifnot(setequal(label_data$cluster, assignments$cluster), nrow(label_data) == nrow(label_parts))
  coverage <- st_drop_geometry(polygons) %>% distinct(GEOID, cluster)
  list(mainland = mainland, boundary = boundary, polygons = polygons, clusters = clusters,
    label_data = label_data, coverage = coverage, not_area = not_area)
}

build_nb_pc_figures <- function(inputs, geometry) {
  with(c(inputs, geometry), {
    map_caption <- "Black outlines and repeated numbers: SES clusters. Gray hatching: no modeled area result."
    plots <- list()
    export <- function(plot, name) plots[[name]] <<- plot
    map_plot <- function(estimates, value, title, subtitle, scale, by = "cluster", caption = map_caption) {
      map_data <- polygons %>% left_join(estimates, by = by)
      stopifnot(!anyNA(map_data[[value]]))
      if (identical(value, "IRR_median")) {
        stopifnot(identical(scale$limits, c(1, 2)), identical(scale$breaks, seq(1, 2, by = 0.2)),
          identical(toupper(unname(scale$map(scale$limits))), c("#FFFFFF", "#B2182B")),
          all(grDevices::col2rgb(scale$map(1)) >= 250))
      } else {
        reference_scale <- (if (identical(value, "P_AR_gt_0")) reference_ar else
          reference_probability)$fill_scale
        stopifnot(identical(scale$map(seq(0, 1, length.out = 101)),
          reference_scale$map(seq(0, 1, length.out = 101))))
      }
      ggplot() + geom_sf(data = boundary, fill = "grey92", colour = "grey60", linewidth = 0.2) +
        missing_result_hatch_layer(mainland) +
        geom_sf(data = map_data, aes(fill = .data[[value]]), colour = "white", linewidth = 0.10) +
        geom_sf(data = clusters, fill = NA, colour = "black", linewidth = 0.30) +
        geom_text_repel(data = label_data, aes(x = x, y = y, label = cluster, size = label_size), seed = 42,
          fontface = "bold", colour = "black", bg.color = "white", bg.r = 0.035,
          force = 0, force_pull = 0, max.iter = 0, max.time = 0,
          box.padding = 0, point.padding = 0, min.segment.length = Inf,
          max.overlaps = Inf) + scale_size_identity() + scale +
        coord_sf(datum = NA, expand = FALSE) +
        labs(title = title, subtitle = subtitle, caption = caption) +
        (if (identical(value, "P_AR_gt_0")) reference_ar$theme else
          if (identical(value, "P_IRR_gt_1")) reference_probability$theme else reference_irr$theme)
    }
    ratio_scale <- function(values, legend) {
      stopifnot(all(is.finite(values)), all(values > 0))
      scale_fill_gradientn(colours = c("white", "#B2182B"),
        limits = c(1, 2), breaks = seq(1, 2, by = 0.2),
        oob = scales::squish,
        labels = function(x) sprintf("%.1f", x), name = legend,
        guide = guide_colourbar(barheight = grid::unit(65, "mm")))
    }
    prob_scale <- function(legend, reference = reference_probability) {
      scale <- reference$fill_scale$clone()
      scale$name <- legend
      scale$labels <- scales::label_percent(accuracy = 1)
      scale
    }

    acute_map <- zip_info %>% mutate(log_irr = log(IRR_median)) %>% select(-cluster)
    heat_subtitle <- paste0("Negative-binomial model | Posterior median IRRs for +", sprintf("%.2f\u00b0F", heat_sd_f),
      " (1 SD)\n", heat_structure, " | ZIP iid + PC hyperpriors")
    heat_caption <- if (has_linear_heat)
      "SD is calculated county-wide. Missing-zone ZIPs retain the shared heat curve and SES-specific slope." else
      "SD is calculated county-wide. Missing-zone ZIPs use the shared heat curve; linear cluster heat slopes are excluded."
    if (any(zip_info$IRR_median < 1 | zip_info$IRR_median > 2))
      heat_caption <- paste(heat_caption, "IRRs outside 1.0-2.0 use the nearest legend end color; estimates are unchanged.")
    heat_caption <- paste(strwrap(heat_caption, width = 110L), collapse = "\n")
    export(map_plot(acute_map, "IRR_median", "Los Angeles acute-heat incidence rate ratios",
      heat_subtitle, ratio_scale(zip_info$IRR_median[zip_info$GEOID %in% polygons$GEOID],
        "Acute-heat IRR\n(1 = no change)"), by = "GEOID",
      caption = heat_caption), "acute_heat_irr_1sd_map")
    export(map_plot(acute_map, "P_IRR_gt_1", "Los Angeles posterior exceedance probabilities",
      paste0("Negative-binomial model | Posterior probability of an increased rate\nfor +", sprintf("%.2f\u00b0F", heat_sd_f),
        " (1 SD) from each ZIP's average | ", heat_structure, " + PC hyperpriors"),
      prob_scale(expression(Pr(IRR["1SD"] > 1 ~ "|" ~ Y))), by = "GEOID",
      caption = heat_caption), "acute_heat_probability_map")

    cluster_order <- as.character(prepared$cluster_vuln$cluster)
    forest_theme <- theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, size = 9), plot.margin = margin(12, 18, 12, 12))
    forest <- function(data, point, lower, upper, title, subtitle, xlab, caption, ratio = TRUE) {
      p <- ggplot(data, aes(y = row)) +
        geom_vline(xintercept = if (ratio) 1 else 0, linetype = "dashed", colour = "grey50") +
        geom_segment(aes(x = .data[[lower]], xend = .data[[upper]], yend = row),
          colour = "#2166AC", linewidth = 0.65) +
        geom_point(aes(x = .data[[point]]), colour = "#2166AC", size = 2.7) +
        labs(title = title, subtitle = subtitle, x = xlab, y = NULL, caption = caption) + forest_theme
      if (ratio) p <- p + scale_x_log10(breaks = scales::breaks_log(n = 6), labels = scales::label_number())
      p
    }
    baseline <- baseline %>% mutate(row = factor(cluster, levels = cluster_order,
      labels = paste0("Cluster ", cluster_order, ifelse(cluster_order == reference, " (reference)", ""))))
    export(forest(baseline, "IRR", "lower", "upper", "Baseline cluster IRRs",
      paste0(model_label, "\nConditional SES contrasts; reference: cluster ", reference),
      "Baseline cluster IRR (log scale)",
      "Points: exp(posterior mean log IRR); bars: 95% credible intervals. Reference IRR = 1 by definition.\nClusters ordered by vulnerability, least vulnerable at bottom."), "baseline_cluster_irr_forest")
    heat <- unique(zip_info[c("cluster", "climate_zone", "IRR_median", "IRR_lower", "IRR_upper", "P_IRR_gt_1", "heat_sd_f")])
    heat$zone_label <- ifelse(is.na(heat$climate_zone), "Climate zone unavailable", paste("Climate zone", heat$climate_zone))
    zone_levels <- c(paste("Climate zone", sort(unique(as.numeric(na.omit(heat$climate_zone))))), "Climate zone unavailable")
    heat$zone_label <- factor(heat$zone_label, levels = zone_levels)
    heat$row <- factor(heat$cluster, levels = cluster_order, labels = paste("Cluster", cluster_order))
    heat_table <- heat %>% select(-row)
    export(forest(heat, "IRR_median", "IRR_lower", "IRR_upper", "Acute-heat IRRs per one SD",
      heat_subtitle, sprintf("Acute-heat IRR for +%.2f\u00b0F (log scale)", heat_sd_f),
      paste0("ZIP mean heat versus one SD above that mean; remaining model terms held fixed.\nPoints: posterior medians; bars: 95% credible intervals. Heat responses vary by ",
        if (has_linear_heat) "SES cluster and climate zone." else "climate zone; clusters in the same zone share the heat curve.")) +
        facet_wrap(~zone_label, ncol = 3),
      "acute_heat_irr_1sd_forest")

    ar_file <- file.path(results_dir, "cluster_excess_rates.rds")
    if (file.exists(ar_file)) {
      ar <- readRDS(ar_file)
      provenance <- readRDS(file.path(results_dir, "provenance.rds"))
      stopifnot(identical(ar$provenance$model_md5, provenance$model_md5), all(ar$summary$n_samples == 2000))
      ar <- ar$summary %>% mutate(cluster = sub("^cluster", "", cluster))
      ar_subtitle <- paste0("Negative-binomial model | baseline PC priors\n",
        "Cluster-average excess ED rate: observed heat versus\neach SES cluster's monthly 2018 mean heat")
      export(map_plot(ar, "P_AR_gt_0", "Probability of Positive Cluster-Average Excess ED Rate",
        ar_subtitle, prob_scale(expression(Pr(AR[c](theta) > 0 ~ "|" ~ Y)), reference_ar), caption = NULL), "excess_rate_probability_map")
      ar$row <- factor(ar$cluster, levels = cluster_order, labels = paste("Cluster", cluster_order))
      export(forest(ar, "AR_mean", "AR_q025", "AR_q975", "Cluster-average signed excess rates",
        ar_subtitle, "Excess visits per 10,000 people per day",
        "Points: posterior means; bars: 95% credible intervals from 2,000 joint posterior draws.\nPositive and negative contributions included; each usable ZIP-day receives equal weight.", ratio = FALSE),
        "excess_rate_forest")
    } else message("Excess-rate sampling still pending; exporting two IRR maps and two forest plots.")
    # Keep the original three-map set; baseline cluster IRRs have their own forest.
    expected_maps <- c("acute_heat_irr_1sd_map", "acute_heat_probability_map", "excess_rate_probability_map")
    stopifnot(all(names(plots)[grepl("_map$", names(plots))] %in% expected_maps))
    if (file.exists(ar_file)) stopifnot(setequal(names(plots), c(expected_maps,
      "baseline_cluster_irr_forest", "acute_heat_irr_1sd_forest", "excess_rate_forest")))
    list(plots = plots, heat_table = heat_table)
  })
}

save_nb_pc_figures <- function(plots, output_dir, figure_names = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  grobs <- lapply(plots, ggplotGrob)
  map_names <- names(grobs)[grepl("_map$", names(grobs))]
  common_widths <- do.call(grid::unit.pmax, lapply(grobs[map_names], function(g) g$widths))
  common_heights <- do.call(grid::unit.pmax, lapply(grobs[map_names], function(g) g$heights))
  save_names <- if (is.null(figure_names)) names(grobs) else figure_names
  stopifnot(length(save_names) > 0L, all(save_names %in% names(grobs)))
  for (name in save_names) {
    if (name %in% map_names) {
      grobs[[name]]$widths <- common_widths
      grobs[[name]]$heights <- common_heights
    }
    for (ext in c("png", "pdf")) {
      target <- file.path(output_dir, paste0(name, ".", ext))
      if (identical(name, "acute_heat_irr_1sd_forest"))
        ggsave(target, grobs[[name]], width = 12, height = 9.5, dpi = map_dpi, bg = "white",
          device = if (ext == "pdf") grDevices::cairo_pdf else NULL) else
        save_project_map(target, grobs[[name]])
    }
    message("Saved ", name)
  }
}

export_nb_pc_figures <- function(output_dir = "maps/nb_pc_prior",
  results_dir = "Los_Angeles_discrete_none/climate_zone_added_nb_pc", figure_names = NULL) {
  .libPaths(c("C:/Users/yysma/AppData/Local/R/win-library/4.4", .libPaths()))
  library(dplyr)
  library(ggplot2)
  library(sf)
  library(ggrepel)
  source("shared_functions.R")
  if (!interactive() && grDevices::dev.cur() == 1L) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  inputs <- load_nb_map_inputs(results_dir)
  geometry <- prepare_nb_map_geometry(inputs$assignments)
  figures <- build_nb_pc_figures(inputs, geometry)
  tables <- list(acute_heat_irr_1sd_by_zip = inputs$zip_info, baseline_cluster_irrs = inputs$baseline,
    not_area_filled_on_mainland_map = geometry$not_area,
    cluster_region_label_positions = geometry$label_data[c("cluster", "x", "y")],
    mapped_zip_cluster_ids = geometry$coverage, acute_heat_irr_1sd_by_cluster_zone = figures$heat_table)
  for (name in names(tables))
    write.csv(tables[[name]], file.path(results_dir, paste0(name, ".csv")), row.names = FALSE)
  save_nb_pc_figures(figures$plots, output_dir, figure_names)
  with(inputs, {
    writeLines(c("BASELINE PC NEGATIVE-BINOMIAL FIGURES", model_label,
      paste("Source model:", file.path(results_dir, "nbinomial_model.rds")),
      paste0("Baseline cluster IRRs: conditional SES contrasts, relative to cluster ", reference, "."),
      paste0("Acute-heat IRRs: ", heat_increment, " above each ZIP's overall mean, versus that mean."),
      if (has_linear_heat) "Acute heat retains common/SES-specific slopes, a shared RW2 curve and correlated climate-zone deviations." else
        "Acute heat uses only the shared RW2 curve and correlated climate-zone deviations; no linear cluster heat slopes.",
      "Excess rates: observed heat versus the SES-cluster monthly mean; signed visits/10,000 people/day.",
      "Maps use the existing mainland extent, 10 by 8 inch canvases, SES outlines and haloed cluster numbers.",
      "Visual reference: colors and themes preserved in Los_Angeles_discrete_none/heat_map_style.rds.",
      "IRR legend is fixed at 1.0-2.0 with ticks every 0.2; probability legends remain fixed at 0-100%.",
      "IRRs outside the display limits use the nearest end color; saved estimates remain unmodified.",
      "IRR uses the positive side of the climate_zone_added palette: white at 1.0 and red (#B2182B) at 2.0.",
      "Probability maps retain their climate_zone_added blue-off-white-red palette with midpoint 50%.",
      "Only three maps are exported: acute-heat IRR, acute-heat IRR exceedance probability, and excess-rate exceedance probability.",
      "Forest plots: baseline cluster IRRs (separately labeled), acute-heat IRRs, and excess rates.",
      "Each PNG has a PDF counterpart. Full diagnostics remain in the source model's folder."),
      file.path(results_dir, "figure_notes.txt"))
  })
}

if (!exists(".map_functions_only", inherits = FALSE) || !isTRUE(.map_functions_only)) {
  if (exists("model_variant", inherits = FALSE) && !identical(model_variant, "nb_pc_prior"))
    stop("Historical plotting modes are retired in this script; use the current NB workflow.")
  export_nb_pc_figures()
}
