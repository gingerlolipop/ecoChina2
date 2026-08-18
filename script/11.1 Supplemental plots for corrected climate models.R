# 11.1 Supplemental visualization: final Plain RF and Plain MF RF
# ==============================================================================
# Uses ONLY the final selected-variable models evaluated by script 5.1.
# Old incorrectly coded models are not used here.
#
# Supplemental outputs:
#   Figure S1 climate/soil zone-level metrics incl. TSS
#   Figure S2 selected predictor counts by zone
#   Figure S3 normal-map zone metrics incl. TSS
#   Figure S4 future maps for Plain RF
#   Figure S5 future maps for Plain MF RF
# ============================================================================== 

library(terra)
library(data.table)
library(ggplot2)

rm(list = ls())
gc()


# 0. Paths and settings ==========================================================

base_dir <- "H:/Jing/ecoChina2"
assessment_dir <- file.path(base_dir, "assessment_var")
result_map_root <- file.path(base_dir, "result maps")
palette_file <- file.path(base_dir, "color_palette_China.csv")
output_root <- file.path(base_dir, "visualization var threshold0.4", "supplement")
figure_dir <- file.path(output_root, "figures")
table_dir <- file.path(output_root, "tables")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

method_order <- c("rf_var", "mf_var")
method_labels <- c(
  rf_var = "Plain RF",
  mf_var = "Plain MF RF"
)

model_zoneID <- c(1:7, 9:50, 52:55)
map_threshold <- 0.4
tie_tol <- 1e-4
novel_value <- 99L

future_order <- c(
  "2011-2040SSP245",
  "2041-2070SSP245",
  "2071-2100SSP245",
  "2011-2040SSP585",
  "2041-2070SSP585",
  "2071-2100SSP585"
)

metric_columns <- c(
  "balanced_accuracy",
  "f1",
  "auc",
  "precision",
  "recall",
  "specificity",
  "tss"
)

metric_labels <- c(
  balanced_accuracy = "Balanced accuracy",
  f1 = "F1",
  auc = "AUC",
  precision = "Precision",
  recall = "Recall",
  specificity = "Specificity",
  tss = "TSS"
)


# 1. Helpers ====================================================================

require_file <- function(file) {
  if (!file.exists(file)) {
    stop("Missing required file: ", file)
  }
  file
}

save_plot <- function(plot_object, filename, width, height) {
  ggsave(
    filename = filename,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

assigned_map_file <- function(method, scenario) {
  file.path(
    result_map_root,
    method,
    paste0(
      "assigned_zone_",
      scenario,
      "_threshold",
      map_threshold,
      "_tol",
      tie_tol,
      "_novel",
      novel_value,
      "_maskNA8_noNovelNormal.tif"
    )
  )
}

plot_zone_map <- function(raster_file, title, palette) {
  require_file(raster_file)
  
  x <- rast(raster_file)
  x_index <- subst(
    x,
    from = palette$zoneID,
    to = seq_len(nrow(palette)),
    others = NA
  )
  
  plot(
    x_index,
    col = palette$COLOR,
    breaks = seq(0.5, nrow(palette) + 0.5, by = 1),
    legend = FALSE,
    axes = FALSE,
    box = FALSE,
    main = title
  )
}


# 2. Read final assessment =======================================================

model_zone_metrics <- fread(
  require_file(
    file.path(
      assessment_dir,
      "model_test_zone_metrics_var.csv"
    )
  )
)

normal_map_zone_metrics <- fread(
  require_file(
    file.path(
      assessment_dir,
      "normal_map_zone_metrics_var.csv"
    )
  )
)

palette <- fread(require_file(palette_file))
palette[, zoneID := as.integer(zoneID)]
palette_map <- palette[zoneID != 8][order(zoneID)]

if (!(novel_value %in% palette_map$zoneID)) {
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

zone_colors <- setNames(
  palette_map[zoneID %in% model_zoneID]$COLOR,
  as.character(palette_map[zoneID %in% model_zoneID]$zoneID)
)


# 3. Figure S1: zone-level climate and soil metrics =============================

metric_long <- melt(
  model_zone_metrics[
    method %in% method_order & zone %in% model_zoneID
  ],
  id.vars = c("niche", "method", "zone"),
  measure.vars = metric_columns,
  variable.name = "metric",
  value.name = "value"
)

metric_long[
  ,
  `:=`(
    niche_label = factor(
      niche,
      levels = c("climate", "soil"),
      labels = c("Climate", "Soil")
    ),
    method_label = factor(
      method_labels[method],
      levels = method_labels[method_order]
    ),
    metric_label = factor(
      metric_labels[metric],
      levels = unname(metric_labels[metric_columns])
    ),
    zone_factor = factor(
      as.character(zone),
      levels = as.character(model_zoneID)
    )
  )
]

figure_s1 <- ggplot(
  metric_long,
  aes(
    x = metric_label,
    y = value,
    color = zone_factor,
    shape = method_label
  )
) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.10,
      dodge.width = 0.5,
      seed = 11
    ),
    size = 1.45,
    alpha = 0.82
  ) +
  facet_wrap(~ niche_label, nrow = 1) +
  scale_color_manual(values = zone_colors, guide = "none") +
  scale_shape_manual(values = c(16, 17)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = NULL,
    y = "Independent-test metric",
    shape = NULL,
    title = "Final selected-variable RF performance by vegetation zone"
  ) +
  theme_bw(base_size = 10.5) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

save_plot(
  figure_s1,
  file.path(figure_dir, "Figure_S1_zone_level_model_metrics_TSS.png"),
  12.5,
  5.2
)


# 4. Figure S2: selected predictor counts =======================================

predictor_dt <- unique(
  model_zone_metrics[
    ,
    .(
      niche,
      method,
      zone,
      n_variables
    )
  ]
)

predictor_dt[
  ,
  `:=`(
    method_label = factor(
      method_labels[method],
      levels = method_labels[method_order]
    ),
    zone_factor = factor(
      as.character(zone),
      levels = as.character(model_zoneID)
    )
  )
]

figure_s2 <- ggplot(
  predictor_dt,
  aes(
    x = zone_factor,
    y = n_variables,
    shape = method_label
  )
) +
  geom_point(
    position = position_dodge(width = 0.45),
    size = 1.8,
    color = "black"
  ) +
  facet_wrap(~ niche, nrow = 2, scales = "free_y") +
  scale_shape_manual(values = c(16, 17)) +
  labs(
    x = "Vegetation zone",
    y = "Selected predictors",
    shape = NULL,
    title = "Number of selected predictors by vegetation zone"
  ) +
  theme_bw(base_size = 9.5) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 90, size = 6),
    panel.grid.minor = element_blank()
  )

save_plot(
  figure_s2,
  file.path(figure_dir, "Figure_S2_selected_predictor_counts.png"),
  12.5,
  7.0
)


# 5. Figure S3: normal-map zone metrics =========================================

map_metric_columns <- c(
  "balanced_accuracy",
  "f1",
  "precision",
  "recall",
  "specificity",
  "tss"
)

map_metric_long <- melt(
  normal_map_zone_metrics[
    method %in% method_order & zone %in% model_zoneID
  ],
  id.vars = c("method", "zone"),
  measure.vars = map_metric_columns,
  variable.name = "metric",
  value.name = "value"
)

map_metric_long[
  ,
  `:=`(
    method_label = factor(
      method_labels[method],
      levels = method_labels[method_order]
    ),
    metric_label = factor(
      metric_labels[metric],
      levels = unname(metric_labels[map_metric_columns])
    ),
    zone_factor = factor(
      as.character(zone),
      levels = as.character(model_zoneID)
    )
  )
]

figure_s3 <- ggplot(
  map_metric_long,
  aes(
    x = value,
    y = zone_factor,
    color = zone_factor,
    shape = method_label
  )
) +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 1.35,
    alpha = 0.85
  ) +
  facet_wrap(~ metric_label, ncol = 3, scales = "free_x") +
  scale_color_manual(values = zone_colors, guide = "none") +
  scale_shape_manual(values = c(16, 17)) +
  labs(
    x = "Normal-map metric",
    y = "Vegetation zone",
    shape = NULL,
    title = "Final normal-map reconstruction metrics by zone"
  ) +
  theme_bw(base_size = 9.2) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 5.8),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

save_plot(
  figure_s3,
  file.path(figure_dir, "Figure_S3_normal_map_zone_metrics_TSS.png"),
  10.5,
  10.0
)


# 6. Figure S4-S5: future maps ===================================================

future_titles <- c(
  "2011-2040 SSP245",
  "2041-2070 SSP245",
  "2071-2100 SSP245",
  "2011-2040 SSP585",
  "2041-2070 SSP585",
  "2071-2100 SSP585"
)

for (method_index in seq_along(method_order)) {
  method_key <- method_order[[method_index]]
  
  output_file <- file.path(
    figure_dir,
    paste0(
      "Figure_S",
      method_index + 3L,
      "_future_maps_",
      method_key,
      ".png"
    )
  )
  
  png(output_file, width = 2800, height = 1900, res = 250)
  
  old_par <- par(no.readonly = TRUE)
  par(
    mfrow = c(2, 3),
    mar = c(1, 1, 2.2, 1),
    oma = c(0, 0, 2.5, 0)
  )
  
  for (scenario_index in seq_along(future_order)) {
    plot_zone_map(
      assigned_map_file(method_key, future_order[[scenario_index]]),
      future_titles[[scenario_index]],
      palette_map
    )
  }
  
  mtext(
    paste0("Future ecosystem maps | ", method_labels[[method_key]]),
    outer = TRUE,
    line = 0.6,
    cex = 1.2
  )
  
  par(old_par)
  dev.off()
}

fwrite(
  model_zone_metrics,
  file.path(table_dir, "supplement_model_zone_metrics_TSS.csv")
)

fwrite(
  normal_map_zone_metrics,
  file.path(table_dir, "supplement_normal_map_zone_metrics_TSS.csv")
)

cat("\nCOMPLETE\n")
