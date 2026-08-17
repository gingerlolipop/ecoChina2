# 11.2 Visualization for corrected climate models
# Threshold 0.4 branch
# ==============================================================================
#
# This supplemental visualization script combines:
#   - original plain_rf and plain_mf results;
#   - corrected rf_var and mf_var results;
#   - threshold-0.4 assigned-zone results;
#   - assigned-zone population/species projections from script 6.1;
#   - continuous dual-suitability projections from script 6.21;
#   - pixel-level ranking and uncertainty from script 8.1.
#
# Required prior scripts:
#   5.1 assessment threshold 0.4
#   6.1 future results threshold 0.4
#   6.21 dual-suitability population/species niches
#   8.1 rank dual suitability by pixel
#
# Outputs:
#   visualization var threshold0.4/
#
# Existing script 11 outputs are not overwritten.
# ==============================================================================

library(terra)
library(data.table)
library(ggplot2)

rm(list = ls())
gc()


# 0. Paths and settings ==========================================================

base_dir <- "H:/Jing/ecoChina2"

assessment_dir <- file.path(
  base_dir,
  "assessment_var"
)

assigned_result_root <- file.path(
  base_dir,
  "future tree niche var"
)

assigned_table_dir <- file.path(
  assigned_result_root,
  "tables"
)

dual_result_root <- file.path(
  base_dir,
  "future tree niche dual suitability var"
)

dual_table_dir <- file.path(
  dual_result_root,
  "tables"
)

ranking_root <- file.path(
  base_dir,
  "dual suit ranking var"
)

ranking_table_dir <- file.path(
  ranking_root,
  "tables"
)

result_map_root <- file.path(
  base_dir,
  "result maps"
)

reference_file <- file.path(
  base_dir,
  "raster/ecosys_ori.tif"
)

palette_file <- file.path(
  base_dir,
  "color_palette_China.csv"
)

output_root <- file.path(
  base_dir,
  "visualization var threshold0.4"
)

figure_dir <- file.path(
  output_root,
  "figures"
)

table_dir <- file.path(
  output_root,
  "tables"
)

assigned_species_page_dir <- file.path(
  figure_dir,
  "assigned species maps"
)

dual_species_page_dir <- file.path(
  figure_dir,
  "dual species maps"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  assigned_species_page_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dual_species_page_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

old_threshold <- 0.2
new_threshold <- 0.4
tie_tol <- 1e-4
novel_value <- 99

comparison_method_order <- c(
  "plain_rf",
  "rf_var",
  "plain_mf",
  "mf_var"
)

corrected_method_order <- c(
  "rf_var",
  "mf_var"
)

method_labels <- c(
  plain_rf = "Original Plain RF",
  rf_var = "Plain RF",
  plain_mf = "Original Plain MF RF",
  mf_var = "Plain MF RF"
)

method_colors <- c(
  plain_rf = "#8FB3D9",
  rf_var = "#2E5E8C",
  plain_mf = "#D5B27B",
  mf_var = "#8C6A2E",
  multiclass_rf = "#4F4F4F"
)

future_order <- c(
  "2011-2040SSP245",
  "2041-2070SSP245",
  "2071-2100SSP245",
  "2011-2040SSP585",
  "2041-2070SSP585",
  "2071-2100SSP585"
)

scenario_order <- c(
  "normal",
  future_order
)

period_levels <- c(
  "2011-2040",
  "2041-2070",
  "2071-2100"
)

ssp_levels <- c(
  "SSP245",
  "SSP585"
)

species_per_page <- 10L
species_page_columns <- 5L
species_page_rows <- 2L

model_zoneID <- c(
  1:7,
  9:50,
  52:55
)

multiclass_reference_map_file <- file.path(
  result_map_root,
  "multiclass_rf",
  "assigned_zone_normal_multiclass_rf.tif"
)

chord_dir <- file.path(
  figure_dir,
  "chord diagrams"
)

dir.create(
  chord_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 1. Helpers ====================================================================

safe_name <- function(x) {
  
  x <- gsub(
    "[^A-Za-z0-9_]+",
    "_",
    x
  )
  
  x <- gsub(
    "^_+|_+$",
    "",
    x
  )
  
  ifelse(
    nchar(x) == 0,
    "unnamed",
    x
  )
}


require_file <- function(file) {
  
  if (!file.exists(file)) {
    stop(
      "Missing required file: ",
      file
    )
  }
  
  file
}


save_plot <- function(
    plot_object,
    filename,
    width,
    height) {
  
  ggsave(
    filename = filename,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )
}


assigned_map_file <- function(
    method,
    scenario) {
  
  threshold <- if (
    method %in% corrected_method_order
  ) {
    new_threshold
  } else {
    old_threshold
  }
  
  file.path(
    result_map_root,
    method,
    paste0(
      "assigned_zone_",
      scenario,
      "_threshold",
      threshold,
      "_tol",
      tie_tol,
      "_novel",
      novel_value,
      "_maskNA8_noNovelNormal.tif"
    )
  )
}


scenario_fields <- function(scenario) {
  
  data.table(
    period = sub(
      "SSP.*$",
      "",
      scenario
    ),
    ssp = sub(
      "^.*(SSP[0-9]+)$",
      "\\1",
      scenario
    )
  )
}


plot_zone_map <- function(
    raster_file,
    title,
    palette) {
  
  if (!file.exists(raster_file)) {
    
    plot.new()
    
    title(
      main = paste0(
        title,
        "\nmissing"
      )
    )
    
    return(invisible(FALSE))
  }
  
  raster <- rast(
    raster_file
  )
  
  indexed <- subst(
    raster,
    from = palette$zoneID,
    to = seq_len(
      nrow(palette)
    ),
    others = NA
  )
  
  plot(
    indexed,
    col = palette$COLOR,
    breaks = seq(
      0.5,
      nrow(palette) + 0.5,
      by = 1
    ),
    legend = FALSE,
    axes = FALSE,
    box = FALSE,
    main = title
  )
  
  invisible(TRUE)
}


plot_binary_species_page <- function(
    species_names,
    raster_files,
    page_title,
    output_file,
    reference_mask,
    fill_color = "#2E8B57") {
  
  png(
    output_file,
    width = 3200,
    height = 1500,
    res = 250
  )
  
  old_par <- par(
    no.readonly = TRUE
  )
  
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  
  par(
    mfrow = c(
      species_page_rows,
      species_page_columns
    ),
    mar = c(
      0.5,
      0.5,
      2,
      0.5
    ),
    oma = c(
      0,
      0,
      2.5,
      0
    )
  )
  
  for (index in seq_along(
    species_names
  )) {
    
    plot(
      reference_mask,
      col = "grey94",
      legend = FALSE,
      axes = FALSE,
      box = FALSE,
      main = species_names[index]
    )
    
    if (file.exists(
      raster_files[index]
    )) {
      
      species_raster <- rast(
        raster_files[index]
      )
      
      plot(
        species_raster,
        add = TRUE,
        col = grDevices::adjustcolor(
          fill_color,
          alpha.f = 0.7
        ),
        legend = FALSE
      )
    }
  }
  
  if (length(species_names) <
      species_per_page) {
    
    for (unused in seq_len(
      species_per_page -
      length(species_names)
    )) {
      plot.new()
    }
  }
  
  mtext(
    page_title,
    outer = TRUE,
    line = 0.5,
    cex = 1.2
  )
  
  invisible(TRUE)
}


global_mean <- function(raster) {
  
  as.numeric(
    global(
      raster,
      fun = "mean",
      na.rm = TRUE
    )[1, 1]
  )
}



div <- function(a, b) {
  
  ifelse(
    is.finite(b) & b > 0,
    a / b,
    NA_real_
  )
}


zone_color_vector <- function(zone_values) {
  
  zone_values <- sort(unique(as.integer(zone_values)))
  zone_values <- zone_values[!is.na(zone_values)]
  
  color_map <- setNames(
    rep(NA_character_, length(zone_values)),
    as.character(zone_values)
  )
  
  matched <- match(
    zone_values,
    palette_map$zoneID
  )
  
  good <- !is.na(matched)
  
  color_map[as.character(zone_values[good])] <-
    palette_map$COLOR[matched[good]]
  
  color_map
}


calculate_zone_metrics_from_maps <- function(original_map, predicted_map) {
  
  original_modeled <- subst(
    original_map,
    from = model_zoneID,
    to = model_zoneID,
    others = NA
  )
  
  predicted_modeled <- subst(
    predicted_map,
    from = model_zoneID,
    to = model_zoneID,
    others = NA
  )
  
  names(original_modeled) <- "original_zone"
  names(predicted_modeled) <- "predicted_zone"
  
  confusion <- as.data.table(
    crosstab(
      c(
        original_modeled,
        predicted_modeled
      ),
      long = TRUE,
      useNA = FALSE
    )
  )
  
  setnames(
    confusion,
    names(confusion),
    c(
      "original_zone",
      "predicted_zone",
      "n"
    )
  )
  
  confusion[
    ,
    `:=`(
      original_zone = as.integer(original_zone),
      predicted_zone = as.integer(predicted_zone),
      n = as.numeric(n)
    )
  ]
  
  total_compared <- sum(
    confusion$n,
    na.rm = TRUE
  )
  
  rbindlist(
    lapply(
      model_zoneID,
      function(zone_value) {
        
        tp <- confusion[
          original_zone == zone_value & predicted_zone == zone_value,
          sum(n, na.rm = TRUE)
        ]
        
        fn <- confusion[
          original_zone == zone_value & predicted_zone != zone_value,
          sum(n, na.rm = TRUE)
        ]
        
        fp <- confusion[
          original_zone != zone_value & predicted_zone == zone_value,
          sum(n, na.rm = TRUE)
        ]
        
        tn <- total_compared - tp - fn - fp
        
        recall_value <- div(tp, tp + fn)
        specificity_value <- div(tn, tn + fp)
        precision_value <- div(tp, tp + fp)
        
        data.table(
          zone = zone_value,
          f1 = div(
            2 * precision_value * recall_value,
            precision_value + recall_value
          ),
          balanced_accuracy = div(
            recall_value + specificity_value,
            2
          )
        )
      }
    ),
    fill = TRUE
  )
}


plot_chord <- function(
    flow,
    item_order,
    item_color,
    item_label,
    out_file,
    main,
    label_cex = 0.65,
    label_track_height = 0.13) {
  
  flow <- as.data.table(copy(flow))
  
  flow[
    ,
    `:=`(
      from = as.character(from),
      to = as.character(to),
      n = as.numeric(n)
    )
  ]
  
  flow <- flow[n > 0]
  
  from_id <- item_order[item_order %in% unique(flow$from)]
  to_id <- rev(item_order[item_order %in% unique(flow$to)])
  
  from_sector <- paste0("O_", from_id)
  to_sector <- paste0("P_", to_id)
  sector_order <- c(from_sector, to_sector)
  
  flow[
    ,
    `:=`(
      from_sector = paste0("O_", from),
      to_sector = paste0("P_", to)
    )
  ]
  
  sector_color <- c(
    setNames(unname(item_color[from_id]), from_sector),
    setNames(unname(item_color[to_id]), to_sector)
  )
  
  sector_label <- c(
    setNames(unname(item_label[from_id]), from_sector),
    setNames(unname(item_label[to_id]), to_sector)
  )
  
  link_color <- unname(item_color[flow$from])
  
  circlize::circos.clear()
  pdf(out_file, width = 16, height = 16, useDingbats = FALSE)
  
  on.exit({
    circlize::circos.clear()
    dev.off()
  }, add = TRUE)
  
  par(mar = c(0.5, 0.5, 2.8, 0.5), xpd = NA)
  
  circlize::circos.par(
    start.degree = 90,
    clock.wise = FALSE,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.002, 0.002),
    canvas.xlim = c(-1.38, 1.38),
    canvas.ylim = c(-1.28, 1.28),
    points.overflow.warning = FALSE
  )
  
  circlize::chordDiagram(
    x = as.data.frame(flow[, .(from_sector, to_sector, n)]),
    order = sector_order,
    grid.col = sector_color,
    grid.border = NA,
    col = link_color,
    transparency = 0.72,
    directional = 1,
    direction.type = "diffHeight",
    diffHeight = circlize::mm_h(1),
    link.target.prop = FALSE,
    link.sort = "default",
    link.decreasing = TRUE,
    link.largest.ontop = TRUE,
    annotationTrack = "grid",
    annotationTrackHeight = circlize::mm_h(2),
    preAllocateTracks = list(track.height = label_track_height),
    big.gap = 14,
    small.gap = 0.15,
    reduce = -1
  )
  
  circlize::circos.trackPlotRegion(
    track.index = 1,
    bg.border = NA,
    panel.fun = function(x, y) {
      
      sector <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      ylim <- circlize::get.cell.meta.data("ylim")
      
      circlize::circos.text(
        x = mean(xlim),
        y = mean(ylim),
        labels = unname(sector_label[sector]),
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0.5, 0.5),
        cex = label_cex
      )
    }
  )
  
  mtext(main, side = 3, line = 0.5, font = 2, cex = 1.25)
  text(-1.27, 0, labels = "Original", srt = 90, font = 2, cex = 1.15)
  text(1.27, 0, labels = "Assigned", srt = 270, font = 2, cex = 1.15)
}


# 2. Read source tables ==========================================================

climate_summary <- fread(
  require_file(
    file.path(
      assessment_dir,
      "climate_test_model_summary_var.csv"
    )
  )
)

map_summary <- fread(
  require_file(
    file.path(
      assessment_dir,
      "normal_map_overall_metrics_var.csv"
    )
  )
)

comparison_delta <- fread(
  require_file(
    file.path(
      assessment_dir,
      "corrected_minus_original_summary_var.csv"
    )
  )
)

ecosystem_area <- fread(
  require_file(
    file.path(
      assigned_table_dir,
      "future_ecosystem_area_var.csv"
    )
  )
)

transition <- fread(
  require_file(
    file.path(
      assigned_table_dir,
      "future_ecosystem_transition_var.csv"
    )
  )
)

assigned_species_area <- fread(
  require_file(
    file.path(
      assigned_table_dir,
      "future_species_niche_area_var.csv"
    )
  )
)

assigned_population_area <- fread(
  require_file(
    file.path(
      assigned_table_dir,
      "future_population_niche_area_var.csv"
    )
  )
)

dual_species_area <- fread(
  require_file(
    file.path(
      dual_table_dir,
      "dual_species_niche_area_var.csv"
    )
  )
)

dual_population_area <- fread(
  require_file(
    file.path(
      dual_table_dir,
      "dual_population_niche_area_var.csv"
    )
  )
)

ranking_index <- fread(
  require_file(
    file.path(
      ranking_table_dir,
      "ranking_output_index_var.csv"
    )
  )
)

palette <- fread(
  require_file(
    palette_file
  )
)

reference_map <- rast(
  require_file(
    reference_file
  )
)

reference_mask <- ifel(
  !is.na(reference_map),
  1,
  NA
)

palette[
  ,
  zoneID := as.integer(zoneID)
]

palette_map <- palette[
  zoneID != 8
][order(zoneID)]

if (!(
  novel_value %in%
  palette_map$zoneID
)) {
  
  palette_map <- rbind(
    palette_map,
    data.table(
      zoneID = novel_value,
      zone = "Novel",
      COLOR = "#333333"
    ),
    fill = TRUE
  )
}

for (table_name in c(
  "ecosystem_area",
  "transition",
  "assigned_species_area",
  "assigned_population_area",
  "dual_species_area",
  "dual_population_area"
)) {
  
  table_object <- get(
    table_name
  )
  
  if ("period" %in%
      names(table_object)) {
    
    table_object[
      ,
      period := factor(
        period,
        levels = period_levels
      )
    ]
  }
  
  if ("ssp" %in%
      names(table_object)) {
    
    table_object[
      ,
      ssp := factor(
        ssp,
        levels = ssp_levels
      )
    ]
  }
  
  assign(
    table_name,
    table_object
  )
}


# 3. Independent climate assessment ============================================

climate_long <- melt(
  climate_summary,
  id.vars = c(
    "method",
    "model_family",
    "model_version"
  ),
  measure.vars = c(
    "mean_balanced_accuracy",
    "mean_f1",
    "mean_auc"
  ),
  variable.name = "metric",
  value.name = "value"
)

climate_long[
  ,
  `:=`(
    method = factor(method, levels = comparison_method_order),
    method_label = factor(
      method_labels[as.character(method)],
      levels = method_labels[comparison_method_order]
    ),
    metric_label = factor(
      metric,
      levels = c("mean_balanced_accuracy", "mean_f1", "mean_auc"),
      labels = c("Balanced accuracy", "F1", "AUC")
    )
  )
]

figure_1 <- ggplot(
  climate_long,
  aes(
    x = metric_label,
    y = value,
    group = method_label,
    color = method_label
  )
) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2.6) +
  scale_color_manual(values = method_colors[comparison_method_order]) +
  coord_cartesian(ylim = c(0.93, 1)) +
  labs(
    x = NULL,
    y = "Mean independent-test metric",
    color = NULL,
    title = "Climate-model assessment"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

save_plot(
  figure_1,
  file.path(figure_dir, "Figure_var_1_climate_assessment.png"),
  8.8,
  4.8
)


# 4. Normal-map assessment =======================================================

map_long <- melt(
  map_summary,
  id.vars = c("method", "model_family", "model_version"),
  measure.vars = c("exact_accuracy", "coverage"),
  variable.name = "metric",
  value.name = "value"
)

map_long[
  ,
  `:=`(
    method = factor(method, levels = comparison_method_order),
    method_label = factor(
      method_labels[as.character(method)],
      levels = method_labels[comparison_method_order]
    ),
    metric_label = factor(
      metric,
      levels = c("exact_accuracy", "coverage"),
      labels = c("Exact agreement", "Coverage")
    )
  )
]

figure_2 <- ggplot(
  map_long,
  aes(
    x = metric_label,
    y = value,
    group = method_label,
    color = method_label
  )
) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2.6) +
  scale_color_manual(values = method_colors[comparison_method_order]) +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(
    x = NULL,
    y = "Normal-map metric",
    color = NULL,
    title = "Normal-period reconstruction"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

save_plot(
  figure_2,
  file.path(figure_dir, "Figure_var_2_normal_map_assessment.png"),
  7.8,
  4.8
)


# 5. Corrected-minus-original deltas ============================================

comparison_long <- melt(
  comparison_delta,
  id.vars = c("model_family", "original_method", "corrected_method"),
  measure.vars = c(
    "delta_mean_balanced_accuracy",
    "delta_mean_f1",
    "delta_mean_auc",
    "delta_normal_exact_accuracy"
  ),
  variable.name = "metric",
  value.name = "delta"
)

comparison_long[
  ,
  metric_label := factor(
    metric,
    levels = c(
      "delta_mean_balanced_accuracy",
      "delta_mean_f1",
      "delta_mean_auc",
      "delta_normal_exact_accuracy"
    ),
    labels = c(
      "Balanced accuracy",
      "F1",
      "AUC",
      "Normal-map agreement"
    )
  )
]

comparison_long[
  ,
  model_family := factor(
    model_family,
    levels = c("single RF", "multi-Forest"),
    labels = c("Plain RF family", "Plain MF RF family")
  )
]

figure_3 <- ggplot(
  comparison_long,
  aes(
    x = metric_label,
    y = delta,
    color = model_family,
    group = model_family
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey40") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.8) +
  scale_color_manual(
    values = c(
      "Plain RF family" = method_colors[["rf_var"]],
      "Plain MF RF family" = method_colors[["mf_var"]]
    )
  ) +
  labs(
    x = NULL,
    y = "Final minus original",
    color = NULL,
    title = "Effect of fixing the plain-model climate variables"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 18, hjust = 1),
    panel.grid.minor = element_blank()
  )

save_plot(
  figure_3,
  file.path(figure_dir, "Figure_var_3_corrected_minus_original.png"),
  9.2,
  4.8
)


# 6. Normal assigned maps ========================================================

normal_files <- c(
  original = reference_file,
  plain_rf = assigned_map_file(
    "plain_rf",
    "normal"
  ),
  rf_var = assigned_map_file(
    "rf_var",
    "normal"
  ),
  plain_mf = assigned_map_file(
    "plain_mf",
    "normal"
  ),
  mf_var = assigned_map_file(
    "mf_var",
    "normal"
  )
)

normal_titles <- c(
  original = "Original vegetation map",
  plain_rf = "Original Plain RF",
  rf_var = "Plain RF",
  plain_mf = "Original Plain MF RF",
  mf_var = "Plain MF RF"
)

png(
  file.path(
    figure_dir,
    "Figure_var_4_normal_maps.png"
  ),
  width = 3000,
  height = 1900,
  res = 250
)

old_par <- par(
  no.readonly = TRUE
)

par(
  mfrow = c(
    2,
    3
  ),
  mar = c(
    1,
    1,
    2.2,
    1
  )
)

for (map_name in names(
  normal_files
)) {
  
  plot_zone_map(
    normal_files[map_name],
    normal_titles[map_name],
    palette_map
  )
}

plot.new()

legend(
  "center",
  legend = paste0(
    palette_map$zoneID,
    ": ",
    palette_map$zone
  ),
  fill = palette_map$COLOR,
  cex = 0.37,
  bty = "n",
  ncol = 2,
  title = "Vegetation zone"
)

par(old_par)
dev.off()


# 7. Future assigned maps ========================================================

future_titles <- c(
  `2011-2040SSP245` = "2011-2040",
  `2041-2070SSP245` = "2041-2070",
  `2071-2100SSP245` = "2071-2100",
  `2011-2040SSP585` = "2011-2040",
  `2041-2070SSP585` = "2041-2070",
  `2071-2100SSP585` = "2071-2100"
)

for (method in corrected_method_order) {
  
  png(
    file.path(
      figure_dir,
      paste0("Figure_var_5_future_maps_", method, ".png")
    ),
    width = 2800,
    height = 1900,
    res = 250
  )
  
  old_par <- par(no.readonly = TRUE)
  
  par(
    mfrow = c(2, 3),
    mar = c(1, 1, 2.4, 1),
    oma = c(0, 0, 2.5, 0)
  )
  
  for (scenario in future_order) {
    
    full_title <- paste(
      sub("SSP.*$", "", scenario),
      sub("^.*(SSP[0-9]+)$", "\\1", scenario)
    )
    
    plot_zone_map(
      assigned_map_file(method, scenario),
      full_title,
      palette_map
    )
  }
  
  mtext(
    paste0("Future ecosystem maps | ", method_labels[method]),
    outer = TRUE,
    line = 0.7,
    cex = 1.2
  )
  
  par(old_par)
  dev.off()
}


# 7b. Normal confusion Sankey and chord diagrams ================================

if (requireNamespace("ggalluvial", quietly = TRUE)) {
  
  normal_map_confusion <- fread(
    require_file(
      file.path(assessment_dir, "normal_map_confusion_long_var.csv")
    )
  )
  
  rf_flow <- normal_map_confusion[
    method == "rf_var" & original_zone != predicted_zone,
    .(count = sum(n)),
    by = .(from = original_zone, to = predicted_zone)
  ]
  
  setorder(rf_flow, -count, from, to)
  rf_flow <- rf_flow[seq_len(min(20L, nrow(rf_flow)))]
  
  if (nrow(rf_flow) > 0) {
    
    rf_flow[
      ,
      `:=`(
        from_chr = as.character(from),
        to_chr = as.character(to)
      )
    ]
    
    figure_4 <- ggplot(
      rf_flow,
      aes(y = count, axis1 = from_chr, axis2 = to_chr)
    ) +
      ggalluvial::geom_alluvium(
        aes(fill = from_chr),
        width = 1 / 12,
        alpha = 0.8,
        show.legend = FALSE
      ) +
      ggalluvial::geom_stratum(
        width = 1 / 8,
        fill = "grey95",
        colour = "grey45"
      ) +
      ggalluvial::stat_stratum(
        geom = "text",
        aes(label = after_stat(stratum)),
        size = 2.5
      ) +
      scale_fill_manual(values = zone_color_vector(rf_flow$from)) +
      scale_x_discrete(
        limits = c("Original", "Assigned"),
        expand = c(0.08, 0.08)
      ) +
      labs(
        x = NULL,
        y = "Pixel count",
        title = "Major normal-map confusion flows | Plain RF"
      ) +
      theme_bw(base_size = 10.5) +
      theme(panel.grid = element_blank())
    
    save_plot(
      figure_4,
      file.path(figure_dir, "Figure_var_4_major_ecotype_confusion_flows.png"),
      8.8,
      6.2
    )
  }
}

if (requireNamespace("circlize", quietly = TRUE)) {
  
  normal_map_confusion <- fread(
    require_file(
      file.path(assessment_dir, "normal_map_confusion_long_var.csv")
    )
  )
  
  category_lookup <- unique(
    palette[
      zoneID %in% model_zoneID,
      .(zoneID, category2 = as.character(category2))
    ],
    by = "zoneID"
  )
  
  zone_to_category <- setNames(
    as.character(category_lookup$category2),
    as.character(category_lookup$zoneID)
  )
  
  zone_order <- as.character(model_zoneID)
  zone_colors <- zone_color_vector(model_zoneID)
  zone_labels <- setNames(as.character(model_zoneID), as.character(model_zoneID))
  
  category_palette <- palette[
    zoneID %in% model_zoneID,
    .(
      first_zone = min(zoneID),
      COLOR = COLOR[which.min(zoneID)]
    ),
    by = category2
  ][order(first_zone)]
  
  category_order <- as.character(category_palette$category2)
  category_colors <- setNames(category_palette$COLOR, category_order)
  category_labels <- setNames(
    gsub("_", " ", category_order, fixed = TRUE),
    category_order
  )
  
  for (method_name in corrected_method_order) {
    
    method_dt <- normal_map_confusion[
      method == method_name &
        original_zone %in% model_zoneID &
        predicted_zone %in% model_zoneID &
        n > 0
    ]
    
    zone_flow <- method_dt[
      ,
      .(n = sum(n)),
      by = .(
        from = as.character(original_zone),
        to = as.character(predicted_zone)
      )
    ]
    
    plot_chord(
      zone_flow,
      zone_order,
      zone_colors,
      zone_labels,
      file.path(chord_dir, paste0("normal_map_zone_chord_", method_name, ".pdf")),
      paste0(method_labels[method_name], ": Original zone to assigned zone")
    )
    
    category_flow <- copy(zone_flow)
    category_flow[, from := unname(zone_to_category[from])]
    category_flow[, to := unname(zone_to_category[to])]
    category_flow <- category_flow[
      !is.na(from) & !is.na(to),
      .(n = sum(n)),
      by = .(from, to)
    ]
    
    plot_chord(
      category_flow,
      category_order,
      category_colors,
      category_labels,
      file.path(chord_dir, paste0("normal_map_category_chord_", method_name, ".pdf")),
      paste0(method_labels[method_name], ": Original category to assigned category"),
      label_cex = 0.88,
      label_track_height = 0.18
    )
  }
}


# 8. Novel area and ecosystem transitions =======================================

novel_area <- ecosystem_area[zoneID == novel_value]

novel_area[
  ,
  method_label := factor(
    method_labels[method],
    levels = method_labels[corrected_method_order]
  )
]

figure_6 <- ggplot(
  novel_area,
  aes(x = period, y = area_km2, group = method_label, color = method_label)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.1) +
  scale_color_manual(values = method_colors[corrected_method_order]) +
  facet_wrap(~ ssp, nrow = 1, scales = "free_y") +
  labs(
    x = "Future period",
    y = expression("Novel area (km"^2*")"),
    color = NULL,
    title = "Projected novel ecosystem area (threshold 0.4)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

save_plot(
  figure_6,
  file.path(figure_dir, "Figure_var_6_novel_area.png"),
  8.5,
  4.8
)

transition_summary <- transition[
  ,
  .(area_km2 = sum(area_km2, na.rm = TRUE)),
  by = .(method, scenario, period, ssp, transition_type)
]

transition_summary[
  ,
  total_area_km2 := sum(area_km2, na.rm = TRUE),
  by = .(method, scenario)
]

transition_summary[, area_share := area_km2 / total_area_km2]

transition_summary[
  ,
  `:=`(
    method_label = factor(
      method_labels[method],
      levels = method_labels[corrected_method_order]
    ),
    transition_type = factor(
      transition_type,
      levels = c("stable", "changed", "novel"),
      labels = c("Stable zone", "Changed zone", "Novel")
    )
  )
]

figure_7 <- ggplot(
  transition_summary,
  aes(x = period, y = area_share, fill = transition_type)
) +
  geom_col(width = 0.55, color = "white", linewidth = 0.15) +
  facet_grid(method_label ~ ssp) +
  scale_fill_manual(
    values = c(
      "Stable zone" = "#4C9F70",
      "Changed zone" = "#D19149",
      "Novel" = "#333333"
    )
  ) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  labs(
    x = "Future period",
    y = "Share of mapped area",
    fill = NULL,
    title = "Normal-to-future ecosystem transitions"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.spacing = grid::unit(1.1, "lines")
  )

save_plot(
  figure_7,
  file.path(figure_dir, "Figure_var_7_transition_share.png"),
  9.8,
  6.0
)


# 9. Assigned-zone species and population areas =================================

assigned_species_area[
  ,
  method_label := factor(
    method_labels[method],
    levels = method_labels[corrected_method_order]
  )
]

figure_8 <- ggplot(
  assigned_species_area,
  aes(
    x = period,
    y = future_area_km2,
    group = method_label,
    color = method_label
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.7) +
  scale_color_manual(values = method_colors[corrected_method_order]) +
  facet_grid(Species ~ ssp, scales = "free_y") +
  labs(
    x = "Future period",
    y = expression("Assigned-zone species area (km"^2*")"),
    color = NULL,
    title = "Species niches from assigned ecosystem zones"
  ) +
  theme_bw(base_size = 9.2) +
  theme(
    legend.position = "top",
    strip.text.y = element_text(angle = 0),
    panel.grid.minor = element_blank(),
    panel.spacing = grid::unit(0.7, "lines")
  )

save_plot(
  figure_8,
  file.path(figure_dir, "Figure_var_8_assigned_species_area.png"),
  10.4,
  max(8, 1.15 * length(unique(assigned_species_area$Species)))
)


# 10. Continuous dual-suitability niches ========================================

dual_species_long <- melt(
  dual_species_area,
  id.vars = c("Species", "method", "scenario"),
  measure.vars = c("suitable_area_km2", "suitability_weighted_area_km2"),
  variable.name = "area_metric",
  value.name = "area_value"
)

dual_species_long[, fields := scenario_fields(scenario)$period, by = scenario]

dual_species_long[
  ,
  period := factor(sub("SSP.*$", "", scenario), levels = period_levels)
]

dual_species_long[
  ,
  ssp := factor(sub("^.*(SSP[0-9]+)$", "\\1", scenario), levels = ssp_levels)
]

dual_species_long[
  ,
  `:=`(
    method_label = factor(
      method_labels[method],
      levels = method_labels[corrected_method_order]
    ),
    area_metric_label = factor(
      area_metric,
      levels = c("suitable_area_km2", "suitability_weighted_area_km2"),
      labels = c("Area with dual suitability >= 0.4", "Suitability-weighted area")
    )
  )
]

figure_9 <- ggplot(
  dual_species_long,
  aes(x = period, y = area_value, group = method_label, color = method_label)
) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1.5) +
  scale_color_manual(values = method_colors[corrected_method_order]) +
  facet_grid(Species + area_metric_label ~ ssp, scales = "free_y") +
  labs(
    x = "Future period",
    y = expression("Area metric (km"^2*")"),
    color = NULL,
    title = "Species niches from continuous dual suitability"
  ) +
  theme_bw(base_size = 8.8) +
  theme(
    legend.position = "top",
    strip.text.y = element_text(angle = 0),
    panel.grid.minor = element_blank(),
    panel.spacing = grid::unit(0.7, "lines")
  )

save_plot(
  figure_9,
  file.path(figure_dir, "Figure_var_9_dual_species_area.png"),
  11.2,
  max(10, 1.85 * length(unique(dual_species_long$Species)))
)

dual_population_area[
  ,
  period := factor(sub("SSP.*$", "", scenario), levels = period_levels)
]

dual_population_area[
  ,
  ssp := factor(sub("^.*(SSP[0-9]+)$", "\\1", scenario), levels = ssp_levels)
]

dual_population_area[
  ,
  method_label := factor(
    method_labels[method],
    levels = method_labels[corrected_method_order]
  )
]

figure_10a <- ggplot(
  dual_population_area,
  aes(
    x = period,
    y = suitable_area_km2,
    group = interaction(method_label, PopulationID),
    color = method_label
  )
) +
  geom_line(linewidth = 0.45, alpha = 0.45) +
  geom_point(size = 0.8, alpha = 0.45) +
  scale_color_manual(values = method_colors[corrected_method_order]) +
  facet_grid(Species ~ ssp, scales = "free_y") +
  labs(
    x = "Future period",
    y = expression("Population area with dual suitability >= 0.4 (km"^2*")"),
    color = NULL,
    title = "Population niches from source-zone dual suitability"
  ) +
  theme_bw(base_size = 8.8) +
  theme(
    legend.position = "top",
    strip.text.y = element_text(angle = 0),
    panel.grid.minor = element_blank(),
    panel.spacing = grid::unit(0.7, "lines")
  )

save_plot(
  figure_10a,
  file.path(figure_dir, "Figure_var_10a_dual_population_area.png"),
  11.2,
  max(8, 1.25 * length(unique(dual_population_area$Species)))
)

normal_map_zone_metrics <- fread(
  require_file(file.path(assessment_dir, "normal_map_zone_metrics_var.csv"))
)

bubble_dt <- normal_map_zone_metrics[
  method %in% corrected_method_order,
  .(method, zone, f1)
]

if (file.exists(multiclass_reference_map_file)) {
  
  multiclass_map <- rast(multiclass_reference_map_file)
  
  if (compareGeom(reference_map, multiclass_map, stopOnError = FALSE)) {
    
    multiclass_metrics <- calculate_zone_metrics_from_maps(
      reference_map,
      multiclass_map
    )
    
    multiclass_metrics[, method := "multiclass_rf"]
    
    bubble_dt <- rbind(
      bubble_dt,
      multiclass_metrics[, .(method, zone, f1)],
      fill = TRUE
    )
  }
}

ref_area <- cellSize(reference_map, unit = "km")
ref_modeled <- subst(reference_map, from = model_zoneID, to = model_zoneID, others = NA)
area_dt <- as.data.table(zonal(ref_area, ref_modeled, fun = "sum", na.rm = TRUE))
setnames(area_dt, names(area_dt)[1:2], c("zone", "area_km2"))
area_dt[, `:=`(zone = as.integer(zone), area_km2 = as.numeric(area_km2))]

bubble_dt <- merge(bubble_dt, area_dt, by = "zone", all.x = TRUE, sort = FALSE)

bubble_method_levels <- c(
  corrected_method_order,
  if ("multiclass_rf" %in% bubble_dt$method) "multiclass_rf" else character(0)
)

bubble_method_labels <- c(method_labels, multiclass_rf = "Multiclass RF")

bubble_dt[
  ,
  `:=`(
    method_label = factor(
      bubble_method_labels[method],
      levels = bubble_method_labels[bubble_method_levels]
    ),
    zone_factor = factor(
      as.character(zone),
      levels = rev(as.character(sort(unique(zone))))
    )
  )
]

bubble_fill_values <- c(
  setNames(
    unname(method_colors[corrected_method_order]),
    method_labels[corrected_method_order]
  ),
  "Multiclass RF" = method_colors[["multiclass_rf"]]
)

figure_10b <- ggplot(
  bubble_dt,
  aes(x = f1, y = zone_factor, size = area_km2, fill = method_label)
) +
  geom_point(shape = 21, color = "grey20", alpha = 0.8) +
  scale_fill_manual(values = bubble_fill_values) +
  facet_wrap(~ method_label, nrow = 1) +
  labs(
    x = "Zone-level F1",
    y = "Vegetation zone",
    size = expression("Original area (km"^2*")"),
    fill = NULL,
    title = "Reference-map F1 by vegetation zone"
  ) +
  theme_bw(base_size = 9.5) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.spacing = grid::unit(1.7, "lines"),
    axis.text.y = element_text(size = 6.4)
  )

save_plot(
  figure_10b,
  file.path(figure_dir, "Figure_var_10b_reference_map_F1_bubble.png"),
  13.5,
  7.0
)


# 11. Pixel-level ranking and uncertainty =======================================

ranking_summary_results <- list()

valid_ranking_jobs <- ranking_index[
  status %in% c("created", "reused")
]

for (row_index in seq_len(nrow(valid_ranking_jobs))) {
  
  method <- valid_ranking_jobs$method[row_index]
  scenario <- valid_ranking_jobs$scenario[row_index]
  summary_file <- valid_ranking_jobs$ranked_summary_file[row_index]
  
  if (!file.exists(summary_file)) {
    next
  }
  
  summary_raster <- rast(summary_file)
  
  if (nlyr(summary_raster) != 6L) {
    next
  }
  
  names(summary_raster) <- c(
    "n_zone_ranked",
    "n_zone_above_threshold",
    "top1_minus_top2",
    "top1_suit",
    "top2_suit",
    "novel_by_threshold"
  )
  
  fields <- if (scenario == "normal") {
    data.table(period = "normal", ssp = "normal")
  } else {
    scenario_fields(scenario)
  }
  
  ranking_summary_results[[length(ranking_summary_results) + 1L]] <- data.table(
    method = method,
    scenario = scenario,
    period = fields$period,
    ssp = fields$ssp,
    mean_n_zone_ranked = global_mean(summary_raster[["n_zone_ranked"]]),
    mean_n_zone_above_threshold = global_mean(summary_raster[["n_zone_above_threshold"]]),
    mean_top1_minus_top2 = global_mean(summary_raster[["top1_minus_top2"]]),
    mean_top1_suit = global_mean(summary_raster[["top1_suit"]]),
    novel_share = global_mean(summary_raster[["novel_by_threshold"]]),
    ranked_summary_file = summary_file
  )
}

ranking_summary <- rbindlist(ranking_summary_results, fill = TRUE)

fwrite(ranking_summary, file.path(table_dir, "pixel_ranking_summary_var.csv"))

ranking_future <- ranking_summary[scenario != "normal"]

ranking_future[
  ,
  `:=`(
    period = factor(period, levels = period_levels),
    ssp = factor(ssp, levels = ssp_levels),
    method_label = factor(
      method_labels[method],
      levels = method_labels[corrected_method_order]
    )
  )
]

figure_11a <- ggplot(
  ranking_future,
  aes(x = period, y = novel_share, group = method_label, color = method_label)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.0) +
  scale_color_manual(values = method_colors[corrected_method_order]) +
  facet_wrap(~ ssp, nrow = 1) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100, 1), "%")) +
  labs(
    x = "Future period",
    y = "Pixels with all zones < 0.4",
    color = NULL,
    title = "Pixel-level novel share from ranked dual suitability"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

save_plot(
  figure_11a,
  file.path(figure_dir, "Figure_var_11a_ranking_novel_share.png"),
  8.8,
  4.8
)

ranking_margin_long <- melt(
  ranking_future,
  id.vars = c("method", "method_label", "scenario", "period", "ssp"),
  measure.vars = c(
    "mean_top1_suit",
    "mean_top1_minus_top2",
    "mean_n_zone_above_threshold"
  ),
  variable.name = "metric",
  value.name = "value"
)

ranking_margin_long[
  ,
  metric_label := factor(
    metric,
    levels = c(
      "mean_top1_suit",
      "mean_top1_minus_top2",
      "mean_n_zone_above_threshold"
    ),
    labels = c(
      "Mean top suitability",
      "Mean top1 - top2",
      "Mean zones >= 0.4"
    )
  )
]

figure_11b <- ggplot(
  ranking_margin_long,
  aes(x = period, y = value, group = method_label, color = method_label)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_color_manual(values = method_colors[corrected_method_order]) +
  facet_grid(metric_label ~ ssp, scales = "free_y") +
  labs(
    x = "Future period",
    y = NULL,
    color = NULL,
    title = "Pixel-level ranking strength and uncertainty"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

save_plot(
  figure_11b,
  file.path(figure_dir, "Figure_var_11b_ranking_uncertainty.png"),
  9.8,
  7.0
)


# 12. Multi-page species maps ====================================================

species_names <- sort(
  unique(
    assigned_species_area$Species
  )
)

number_of_pages <- ceiling(
  length(species_names) /
    species_per_page
)

for (method in corrected_method_order) {
  for (scenario in future_order) {
    for (page_index in seq_len(
      number_of_pages
    )) {
      
      start_index <-
        (page_index - 1L) *
        species_per_page +
        1L
      
      end_index <- min(
        page_index *
          species_per_page,
        length(species_names)
      )
      
      page_species <- species_names[
        start_index:end_index
      ]
      
      assigned_files <- file.path(
        assigned_result_root,
        method,
        scenario,
        "species niche",
        paste0(
          page_species,
          "_species_niche.tif"
        )
      )
      
      dual_files <- file.path(
        dual_result_root,
        method,
        scenario,
        "species niche binary",
        paste0(
          safe_name(page_species),
          "_species_dual_binary_threshold",
          new_threshold,
          ".tif"
        )
      )
      
      plot_binary_species_page(
        page_species,
        assigned_files,
        paste0(
          "Assigned-zone species niches | ",
          method_labels[method],
          " | ",
          scenario,
          " | page ",
          page_index
        ),
        file.path(
          assigned_species_page_dir,
          paste0(
            "assigned_species_",
            method,
            "_",
            scenario,
            "_page",
            page_index,
            ".png"
          )
        ),
        reference_mask,
        fill_color = "#2E8B57"
      )
      
      plot_binary_species_page(
        page_species,
        dual_files,
        paste0(
          "Dual-suitability species niches >= 0.4 | ",
          method_labels[method],
          " | ",
          scenario,
          " | page ",
          page_index
        ),
        file.path(
          dual_species_page_dir,
          paste0(
            "dual_species_",
            method,
            "_",
            scenario,
            "_page",
            page_index,
            ".png"
          )
        ),
        reference_mask,
        fill_color = "#2E5E8C"
      )
    }
  }
}


# 13. Save figure-source tables ==================================================

fwrite(
  climate_long,
  file.path(
    table_dir,
    "figure_climate_assessment_data_var.csv"
  )
)

fwrite(
  map_long,
  file.path(
    table_dir,
    "figure_normal_map_assessment_data_var.csv"
  )
)

fwrite(
  comparison_long,
  file.path(
    table_dir,
    "figure_corrected_minus_original_data_var.csv"
  )
)

fwrite(
  novel_area,
  file.path(
    table_dir,
    "figure_novel_area_data_var.csv"
  )
)

fwrite(
  transition_summary,
  file.path(
    table_dir,
    "figure_transition_share_data_var.csv"
  )
)

fwrite(
  dual_species_long,
  file.path(
    table_dir,
    "figure_dual_species_area_data_var.csv"
  )
)

fwrite(
  ranking_margin_long,
  file.path(
    table_dir,
    "figure_ranking_uncertainty_data_var.csv"
  )
)

cat(
  "\nCOMPLETE\n",
  "Threshold: ",
  new_threshold,
  "\n",
  "Methods: ",
  paste(
    corrected_method_order,
    collapse = ", "
  ),
  "\n",
  "Figures: ",
  figure_dir,
  "\n",
  "Tables: ",
  table_dir,
  "\n",
  sep = ""
)
