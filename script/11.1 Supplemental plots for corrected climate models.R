# 11.1 Supplemental plots for corrected climate models
# ==============================================================================
# Produces comparison and future-result figures for:
#
#   original models:
#     plain_rf
#     plain_mf
#
#   corrected models:
#     rf_var
#     mf_var
#
# Required prior scripts:
#   5.1 Assessment for corrected climate models
#   6.1 Future results for rf_var and mf_var
#
# New outputs are written to:
#   H:/Jing/ecoChina2/visualization_var
#
# Existing visualization outputs are not overwritten.
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

future_root <- file.path(
  base_dir,
  "future tree niche var"
)

future_table_dir <- file.path(
  future_root,
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
  "visualization_var"
)

figure_dir <- file.path(
  output_root,
  "figures"
)

table_dir <- file.path(
  output_root,
  "tables"
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

old_threshold <- 0.2
new_threshold <- 0.4
tie_tol <- 1e-4
novel_value <- 99

method_order <- c(
  "plain_rf",
  "rf_var",
  "plain_mf",
  "mf_var"
)

method_labels <- c(
  plain_rf = "Original single RF",
  rf_var = "Corrected single RF",
  plain_mf = "Original multi-Forest",
  mf_var = "Corrected multi-Forest"
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


# 1. Helpers ====================================================================

require_file <- function(file) {
  
  if (!file.exists(file)) {
    stop(
      "Missing required file: ",
      file
    )
  }
  
  file
}


assigned_map_file <- function(
    method,
    scenario) {
  
  threshold <- if (
    method %in% c(
      "rf_var",
      "mf_var"
    )
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


recode_zone_raster <- function(
    raster,
    palette) {
  
  zone_values <- palette$zoneID
  
  subst(
    raster,
    from = zone_values,
    to = seq_along(zone_values),
    others = NA
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
    
    return(invisible(NULL))
  }
  
  x <- rast(
    raster_file
  )
  
  x_index <- recode_zone_raster(
    x,
    palette
  )
  
  plot(
    x_index,
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
  
  invisible(NULL)
}


# 2. Read tables ================================================================

climate_summary_file <- require_file(
  file.path(
    assessment_dir,
    "climate_test_model_summary_var.csv"
  )
)

map_summary_file <- require_file(
  file.path(
    assessment_dir,
    "normal_map_overall_metrics_var.csv"
  )
)

comparison_file <- require_file(
  file.path(
    assessment_dir,
    "corrected_minus_original_summary_var.csv"
  )
)

ecosystem_area_file <- require_file(
  file.path(
    future_table_dir,
    "future_ecosystem_area_var.csv"
  )
)

transition_file <- require_file(
  file.path(
    future_table_dir,
    "future_ecosystem_transition_var.csv"
  )
)

species_area_file <- require_file(
  file.path(
    future_table_dir,
    "future_species_niche_area_var.csv"
  )
)

palette_file <- require_file(
  palette_file
)

reference_file <- require_file(
  reference_file
)

climate_summary <- fread(
  climate_summary_file
)

map_summary <- fread(
  map_summary_file
)

comparison_delta <- fread(
  comparison_file
)

ecosystem_area <- fread(
  ecosystem_area_file
)

transition <- fread(
  transition_file
)

species_area <- fread(
  species_area_file
)

palette <- fread(
  palette_file
)

palette[
  ,
  zoneID := as.integer(zoneID)
]

# Hide Zone 8 from legends and maps, consistent with the assignment mask.
palette_map <- palette[
  zoneID != 8
][order(zoneID)]

for (table_name in c(
  "ecosystem_area",
  "transition",
  "species_area"
)) {
  
  table_object <- get(
    table_name
  )
  
  if ("period" %in% names(table_object)) {
    table_object[
      ,
      period := factor(
        period,
        levels = period_levels
      )
    ]
  }
  
  if ("ssp" %in% names(table_object)) {
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


# 3. Figure 1: independent climate assessment ==================================

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
      levels = method_order
    ),
    method_label = factor(
      method_labels[
        as.character(method)
      ],
      levels = method_labels[
        method_order
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
    y = "Mean test-set metric",
    fill = "Model version",
    title = "Independent climate-model assessment"
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
  width = 11,
  height = 4.5
)


# 4. Figure 2: normal-map assessment ===========================================

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
      levels = method_order
    ),
    method_label = factor(
      method_labels[
        as.character(method)
      ],
      levels = method_labels[
        method_order
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
    fill = "Model version",
    title = "Normal-period map reconstruction"
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
  width = 8.5,
  height = 4.5
)


# 5. Figure 3: corrected-minus-original changes =================================

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
    title = "Effect of using mcRFop variables in plain climate models"
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
  width = 8.5,
  height = 4.8
)


# 6. Figure 4: normal assigned maps =============================================

normal_map_files <- c(
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

normal_map_titles <- c(
  original = "Original vegetation map",
  plain_rf = "Original single RF",
  rf_var = "Corrected single RF",
  plain_mf = "Original multi-Forest",
  mf_var = "Corrected multi-Forest"
)

normal_map_figure_file <- file.path(
  figure_dir,
  "Figure_var_4_normal_maps.png"
)

png(
  normal_map_figure_file,
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
  normal_map_files
)) {
  
  plot_zone_map(
    normal_map_files[
      map_name
    ],
    normal_map_titles[
      map_name
    ],
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
  cex = 0.38,
  bty = "n",
  ncol = 2,
  title = "Vegetation zone"
)

par(old_par)
dev.off()


# 7. Future map panels for each corrected method ================================

future_scenarios <- c(
  "normal",
  "2011-2040SSP245",
  "2041-2070SSP245",
  "2071-2100SSP245",
  "2011-2040SSP585",
  "2041-2070SSP585",
  "2071-2100SSP585"
)

future_titles <- c(
  normal = "Normal",
  `2011-2040SSP245` = "2011-2040 SSP245",
  `2041-2070SSP245` = "2041-2070 SSP245",
  `2071-2100SSP245` = "2071-2100 SSP245",
  `2011-2040SSP585` = "2011-2040 SSP585",
  `2041-2070SSP585` = "2041-2070 SSP585",
  `2071-2100SSP585` = "2071-2100 SSP585"
)

for (method in c(
  "rf_var",
  "mf_var"
)) {
  
  output_file <- file.path(
    figure_dir,
    paste0(
      "Figure_var_5_future_maps_",
      method,
      ".png"
    )
  )
  
  png(
    output_file,
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
  
  for (scenario in future_scenarios) {
    
    plot_zone_map(
      assigned_map_file(
        method,
        scenario
      ),
      future_titles[
        scenario
      ],
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
    title = method_labels[
      method
    ]
  )
  
  par(old_par)
  dev.off()
}


# 8. Future novel-zone area ======================================================

novel_area <- ecosystem_area[
  zoneID == novel_value
]

if (nrow(novel_area) > 0) {
  
  novel_area[
    ,
    method_label := factor(
      method_labels[
        method
      ],
      levels = method_labels[
        c(
          "rf_var",
          "mf_var"
        )
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
      title = "Projected novel ecosystem area (threshold 0.4)"
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
      "Figure_var_6_future_novel_area.png"
    ),
    width = 8.5,
    height = 4.5
  )
}


# 9. Stable, changed and novel transition shares ================================

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
  area_share := area_km2 /
    total_area_km2
]

transition_summary[
  ,
  `:=`(
    method_label = factor(
      method_labels[
        method
      ],
      levels = method_labels[
        c(
          "rf_var",
          "mf_var"
        )
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

fwrite(
  transition_summary,
  file.path(
    table_dir,
    "future_transition_summary_var.csv"
  )
)

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
    title = "Normal-to-future ecosystem transitions"
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
    "Figure_var_7_future_transition_share.png"
  ),
  width = 9.5,
  height = 6
)


# 10. Species future suitable area ==============================================

species_area[
  ,
  method_label := factor(
    method_labels[
      method
    ],
    levels = method_labels[
      c(
        "rf_var",
        "mf_var"
      )
    ]
  )
]

figure_8 <- ggplot(
  species_area,
  aes(
    x = period,
    y = future_area_km2,
    group = method_label,
    linetype = method_label
  )
) +
  geom_line(
    linewidth = 0.7
  ) +
  geom_point(
    size = 1.5
  ) +
  facet_grid(
    Species ~ ssp,
    scales = "free_y"
  ) +
  labs(
    x = "Future period",
    y = expression(
      "Projected species niche area (km"^2*")"
    ),
    linetype = "Method",
    title = "Projected species niche area"
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
    "Figure_var_8_future_species_area.png"
  ),
  width = 10,
  height = max(
    8,
    1.2 * length(
      unique(
        species_area$Species
      )
    )
  )
)


# 11. Save plotting-source tables ===============================================

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
  novel_area,
  file.path(
    table_dir,
    "figure_future_novel_area_data_var.csv"
  )
)

fwrite(
  species_area,
  file.path(
    table_dir,
    "figure_future_species_area_data_var.csv"
  )
)

cat(
  "\nCOMPLETE\n",
  "Figures: ",
  figure_dir,
  "\n",
  "Figure-source tables: ",
  table_dir,
  "\n",
  sep = ""
)
