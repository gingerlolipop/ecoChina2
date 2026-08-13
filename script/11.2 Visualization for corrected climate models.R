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
  plain_rf = "Original single RF",
  rf_var = "Corrected single RF",
  plain_mf = "Original multi-Forest",
  mf_var = "Corrected multi-Forest"
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
      COLOR = "#BDBDBD"
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
    method = factor(
      method,
      levels =
        comparison_method_order
    ),
    method_label = factor(
      method_labels[
        as.character(method)
      ],
      levels = method_labels[
        comparison_method_order
      ]
    ),
    metric_label = factor(
      metric,
      levels = c(
        "mean_balanced_accuracy",
        "mean_f1",
        "mean_auc"
      ),
      labels = c(
        "Balanced accuracy",
        "F1",
        "AUC"
      )
    )
  )
]

figure_1 <- ggplot(
  climate_long,
  aes(
    x = method_label,
    y = value,
    fill = model_version
  )
) +
  geom_col(
    width = 0.72
  ) +
  facet_wrap(
    ~ metric_label,
    nrow = 1
  ) +
  coord_cartesian(
    ylim = c(
      0,
      1
    )
  ) +
  labs(
    x = NULL,
    y = "Mean independent-test metric",
    fill = "Version",
    title = "Climate-model assessment"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    legend.position = "top"
  )

save_plot(
  figure_1,
  file.path(
    figure_dir,
    "Figure_var_1_climate_assessment.png"
  ),
  11,
  4.5
)


# 4. Normal-map assessment =======================================================

map_long <- melt(
  map_summary,
  id.vars = c(
    "method",
    "model_family",
    "model_version"
  ),
  measure.vars = c(
    "exact_accuracy",
    "coverage"
  ),
  variable.name = "metric",
  value.name = "value"
)

map_long[
  ,
  `:=`(
    method = factor(
      method,
      levels =
        comparison_method_order
    ),
    method_label = factor(
      method_labels[
        as.character(method)
      ],
      levels = method_labels[
        comparison_method_order
      ]
    ),
    metric_label = factor(
      metric,
      levels = c(
        "exact_accuracy",
        "coverage"
      ),
      labels = c(
        "Exact agreement",
        "Coverage"
      )
    )
  )
]

figure_2 <- ggplot(
  map_long,
  aes(
    x = method_label,
    y = value,
    fill = model_version
  )
) +
  geom_col(
    width = 0.72
  ) +
  facet_wrap(
    ~ metric_label,
    nrow = 1
  ) +
  coord_cartesian(
    ylim = c(
      0,
      1
    )
  ) +
  labs(
    x = NULL,
    y = "Normal-map metric",
    fill = "Version",
    title = "Normal-period reconstruction"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    legend.position = "top"
  )

save_plot(
  figure_2,
  file.path(
    figure_dir,
    "Figure_var_2_normal_map_assessment.png"
  ),
  8.5,
  4.5
)


# 5. Corrected-minus-original deltas ============================================

comparison_long <- melt(
  comparison_delta,
  id.vars = c(
    "model_family",
    "original_method",
    "corrected_method"
  ),
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

figure_3 <- ggplot(
  comparison_long,
  aes(
    x = metric_label,
    y = delta,
    fill = model_family
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  geom_col(
    position = position_dodge(
      width = 0.75
    ),
    width = 0.68
  ) +
  labs(
    x = NULL,
    y = "Corrected minus original",
    fill = "Model family",
    title = "Effect of correcting the plain-model variable sets"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    ),
    legend.position = "top"
  )

save_plot(
  figure_3,
  file.path(
    figure_dir,
    "Figure_var_3_corrected_minus_original.png"
  ),
  8.5,
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
  plain_rf = "Original single RF",
  rf_var = "Corrected single RF",
  plain_mf = "Original multi-Forest",
  mf_var = "Corrected multi-Forest"
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
  normal = "Normal",
  `2011-2040SSP245` =
    "2011-2040 SSP245",
  `2041-2070SSP245` =
    "2041-2070 SSP245",
  `2071-2100SSP245` =
    "2071-2100 SSP245",
  `2011-2040SSP585` =
    "2011-2040 SSP585",
  `2041-2070SSP585` =
    "2041-2070 SSP585",
  `2071-2100SSP585` =
    "2071-2100 SSP585"
)

for (method in corrected_method_order) {
  
  png(
    file.path(
      figure_dir,
      paste0(
        "Figure_var_5_future_maps_",
        method,
        ".png"
      )
    ),
    width = 3200,
    height = 1900,
    res = 250
  )
  
  old_par <- par(
    no.readonly = TRUE
  )
  
  par(
    mfrow = c(
      2,
      4
    ),
    mar = c(
      1,
      1,
      2.2,
      1
    )
  )
  
  for (scenario in scenario_order) {
    
    plot_zone_map(
      assigned_map_file(
        method,
        scenario
      ),
      future_titles[scenario],
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
    cex = 0.34,
    bty = "n",
    ncol = 2,
    title = method_labels[method]
  )
  
  par(old_par)
  dev.off()
}


# 8. Novel area and ecosystem transitions =======================================

novel_area <- ecosystem_area[
  zoneID == novel_value
]

novel_area[
  ,
  method_label := factor(
    method_labels[method],
    levels = method_labels[
      corrected_method_order
    ]
  )
]

figure_6 <- ggplot(
  novel_area,
  aes(
    x = period,
    y = area_km2,
    group = method_label,
    linetype = method_label
  )
) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_point(
    size = 2
  ) +
  facet_wrap(
    ~ ssp,
    nrow = 1,
    scales = "free_y"
  ) +
  labs(
    x = "Future period",
    y = expression(
      "Novel area (km"^2*")"
    ),
    linetype = "Method",
    title =
      "Projected novel ecosystem area (threshold 0.4)"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    legend.position = "top"
  )

save_plot(
  figure_6,
  file.path(
    figure_dir,
    "Figure_var_6_novel_area.png"
  ),
  8.5,
  4.5
)

transition_summary <- transition[
  ,
  .(
    area_km2 = sum(
      area_km2,
      na.rm = TRUE
    )
  ),
  by = .(
    method,
    scenario,
    period,
    ssp,
    transition_type
  )
]

transition_summary[
  ,
  total_area_km2 := sum(
    area_km2,
    na.rm = TRUE
  ),
  by = .(
    method,
    scenario
  )
]

transition_summary[
  ,
  area_share :=
    area_km2 /
    total_area_km2
]

transition_summary[
  ,
  `:=`(
    method_label = factor(
      method_labels[method],
      levels = method_labels[
        corrected_method_order
      ]
    ),
    transition_type = factor(
      transition_type,
      levels = c(
        "stable",
        "changed",
        "novel"
      ),
      labels = c(
        "Stable zone",
        "Changed zone",
        "Novel"
      )
    )
  )
]

figure_7 <- ggplot(
  transition_summary,
  aes(
    x = period,
    y = area_share,
    fill = transition_type
  )
) +
  geom_col(
    width = 0.72
  ) +
  facet_grid(
    method_label ~ ssp
  ) +
  scale_y_continuous(
    labels = function(x) {
      paste0(
        round(
          x * 100
        ),
        "%"
      )
    }
  ) +
  labs(
    x = "Future period",
    y = "Share of mapped area",
    fill = "Transition",
    title =
      "Normal-to-future ecosystem transitions"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    legend.position = "top"
  )

save_plot(
  figure_7,
  file.path(
    figure_dir,
    "Figure_var_7_transition_share.png"
  ),
  9.5,
  6
)


# 9. Assigned-zone species and population areas =================================

assigned_species_area[
  ,
  method_label := factor(
    method_labels[method],
    levels = method_labels[
      corrected_method_order
    ]
  )
]

figure_8 <- ggplot(
  assigned_species_area,
  aes(
    x = period,
    y = future_area_km2,
    group = method_label,
    linetype = method_label
  )
) +
  geom_line(
    linewidth = 0.65
  ) +
  geom_point(
    size = 1.4
  ) +
  facet_grid(
    Species ~ ssp,
    scales = "free_y"
  ) +
  labs(
    x = "Future period",
    y = expression(
      "Assigned-zone species area (km"^2*")"
    ),
    linetype = "Method",
    title =
      "Species niches from assigned ecosystem zones"
  ) +
  theme_bw(
    base_size = 9
  ) +
  theme(
    legend.position = "top",
    strip.text.y = element_text(
      angle = 0
    )
  )

save_plot(
  figure_8,
  file.path(
    figure_dir,
    "Figure_var_8_assigned_species_area.png"
  ),
  10,
  max(
    8,
    1.2 *
      length(
        unique(
          assigned_species_area$Species
        )
      )
  )
)


# 10. Continuous dual-suitability niches ========================================

dual_species_long <- melt(
  dual_species_area,
  id.vars = c(
    "Species",
    "method",
    "scenario"
  ),
  measure.vars = c(
    "suitable_area_km2",
    "suitability_weighted_area_km2"
  ),
  variable.name = "area_metric",
  value.name = "area_value"
)

dual_species_long[
  ,
  fields := scenario_fields(
    scenario
  )$period,
  by = scenario
]

dual_species_long[
  ,
  period := factor(
    sub(
      "SSP.*$",
      "",
      scenario
    ),
    levels = period_levels
  )
]

dual_species_long[
  ,
  ssp := factor(
    sub(
      "^.*(SSP[0-9]+)$",
      "\\1",
      scenario
    ),
    levels = ssp_levels
  )
]

dual_species_long[
  ,
  `:=`(
    method_label = factor(
      method_labels[method],
      levels = method_labels[
        corrected_method_order
      ]
    ),
    area_metric_label = factor(
      area_metric,
      levels = c(
        "suitable_area_km2",
        "suitability_weighted_area_km2"
      ),
      labels = c(
        "Area with dual suitability >= 0.4",
        "Suitability-weighted area"
      )
    )
  )
]

figure_9 <- ggplot(
  dual_species_long,
  aes(
    x = period,
    y = area_value,
    group = method_label,
    linetype = method_label
  )
) +
  geom_line(
    linewidth = 0.65
  ) +
  geom_point(
    size = 1.3
  ) +
  facet_grid(
    Species + area_metric_label ~ ssp,
    scales = "free_y"
  ) +
  labs(
    x = "Future period",
    y = expression(
      "Area metric (km"^2*")"
    ),
    linetype = "Method",
    title =
      "Species niches from continuous dual suitability"
  ) +
  theme_bw(
    base_size = 8.5
  ) +
  theme(
    legend.position = "top",
    strip.text.y = element_text(
      angle = 0
    )
  )

save_plot(
  figure_9,
  file.path(
    figure_dir,
    "Figure_var_9_dual_species_area.png"
  ),
  11,
  max(
    10,
    2.0 *
      length(
        unique(
          dual_species_long$Species
        )
      )
  )
)

dual_population_area[
  ,
  period := factor(
    sub(
      "SSP.*$",
      "",
      scenario
    ),
    levels = period_levels
  )
]

dual_population_area[
  ,
  ssp := factor(
    sub(
      "^.*(SSP[0-9]+)$",
      "\\1",
      scenario
    ),
    levels = ssp_levels
  )
]

dual_population_area[
  ,
  method_label := factor(
    method_labels[method],
    levels = method_labels[
      corrected_method_order
    ]
  )
]

figure_10 <- ggplot(
  dual_population_area,
  aes(
    x = period,
    y = suitable_area_km2,
    group = interaction(
      method_label,
      PopulationID
    ),
    linetype = method_label
  )
) +
  geom_line(
    linewidth = 0.55,
    alpha = 0.8
  ) +
  facet_grid(
    Species ~ ssp,
    scales = "free_y"
  ) +
  labs(
    x = "Future period",
    y = expression(
      "Population area with dual suitability >= 0.4 (km"^2*")"
    ),
    linetype = "Method",
    title =
      "Population niches from source-zone dual suitability"
  ) +
  theme_bw(
    base_size = 8.5
  ) +
  theme(
    legend.position = "top",
    strip.text.y = element_text(
      angle = 0
    )
  )

save_plot(
  figure_10,
  file.path(
    figure_dir,
    "Figure_var_10_dual_population_area.png"
  ),
  11,
  max(
    8,
    1.3 *
      length(
        unique(
          dual_population_area$Species
        )
      )
  )
)


# 11. Pixel-level ranking and uncertainty =======================================

ranking_summary_results <- list()

valid_ranking_jobs <- ranking_index[
  status %in% c(
    "created",
    "reused"
  )
]

for (row_index in seq_len(
  nrow(valid_ranking_jobs)
)) {
  
  method <- valid_ranking_jobs$method[
    row_index
  ]
  
  scenario <- valid_ranking_jobs$scenario[
    row_index
  ]
  
  summary_file <-
    valid_ranking_jobs$ranked_summary_file[
      row_index
    ]
  
  if (!file.exists(summary_file)) {
    next
  }
  
  summary_raster <- rast(
    summary_file
  )
  
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
  
  fields <- if (
    scenario == "normal"
  ) {
    data.table(
      period = "normal",
      ssp = "normal"
    )
  } else {
    scenario_fields(
      scenario
    )
  }
  
  ranking_summary_results[[length(ranking_summary_results) + 1L]] <- data.table(
    method = method,
    scenario = scenario,
    period = fields$period,
    ssp = fields$ssp,
    mean_n_zone_ranked =
      global_mean(
        summary_raster[["n_zone_ranked"]]
      ),
    mean_n_zone_above_threshold =
      global_mean(
        summary_raster[["n_zone_above_threshold"]]
      ),
    mean_top1_minus_top2 =
      global_mean(
        summary_raster[["top1_minus_top2"]]
      ),
    mean_top1_suit =
      global_mean(
        summary_raster[["top1_suit"]]
      ),
    novel_share =
      global_mean(
        summary_raster[["novel_by_threshold"]]
      ),
    ranked_summary_file =
      summary_file
  )
}

ranking_summary <- rbindlist(
  ranking_summary_results,
  fill = TRUE
)

fwrite(
  ranking_summary,
  file.path(
    table_dir,
    "pixel_ranking_summary_var.csv"
  )
)

ranking_future <- ranking_summary[
  scenario != "normal"
]

ranking_future[
  ,
  `:=`(
    period = factor(
      period,
      levels = period_levels
    ),
    ssp = factor(
      ssp,
      levels = ssp_levels
    ),
    method_label = factor(
      method_labels[method],
      levels = method_labels[
        corrected_method_order
      ]
    )
  )
]

figure_11a <- ggplot(
  ranking_future,
  aes(
    x = period,
    y = novel_share,
    group = method_label,
    linetype = method_label
  )
) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_point(
    size = 2
  ) +
  facet_wrap(
    ~ ssp,
    nrow = 1
  ) +
  scale_y_continuous(
    labels = function(x) {
      paste0(
        round(
          x * 100,
          1
        ),
        "%"
      )
    }
  ) +
  labs(
    x = "Future period",
    y = "Pixels with all zones < 0.4",
    linetype = "Method",
    title =
      "Pixel-level novel share from ranked dual suitability"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    legend.position = "top"
  )

save_plot(
  figure_11a,
  file.path(
    figure_dir,
    "Figure_var_11a_ranking_novel_share.png"
  ),
  8.5,
  4.5
)

ranking_margin_long <- melt(
  ranking_future,
  id.vars = c(
    "method",
    "method_label",
    "scenario",
    "period",
    "ssp"
  ),
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
  aes(
    x = period,
    y = value,
    group = method_label,
    linetype = method_label
  )
) +
  geom_line(
    linewidth = 0.7
  ) +
  geom_point(
    size = 1.7
  ) +
  facet_grid(
    metric_label ~ ssp,
    scales = "free_y"
  ) +
  labs(
    x = "Future period",
    y = NULL,
    linetype = "Method",
    title =
      "Pixel-level ranking strength and uncertainty"
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
    legend.position = "top"
  )

save_plot(
  figure_11b,
  file.path(
    figure_dir,
    "Figure_var_11b_ranking_uncertainty.png"
  ),
  9.5,
  7
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
