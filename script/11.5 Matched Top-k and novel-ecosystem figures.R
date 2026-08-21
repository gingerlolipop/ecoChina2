# 11.5 Matched Top-k and novel-ecosystem figures
# ==============================================================================
# Run after script 8.3. This script reads existing CSV files only; it does not
# read rasters or repeat any model, projection, suitability, or area analysis.
#
# Main figures
# ------------
# 1. Matched Top-k change composition. Reference and future maps use the same
#    Top-k assignment rule in each row.
# 2. Novel ecosystem area. This is the Zone 99 result (all 53 current ecotypes
#    have dual suitability < 0.4), so it is independent of Top-k and is shown
#    once rather than repeated as Top-1, Top-3, and Top-5.
# ==============================================================================

library(data.table)
library(ggplot2)

rm(list = ls())
gc()


# 0. Paths and settings =========================================================

base_dir <- "H:/Jing/ecoChina2"

matched_file <- file.path(
  base_dir,
  "assessment_var",
  "future_topk_matched",
  "matched_topk_change_summary.csv"
)

existing_novel_file <- file.path(
  base_dir,
  "visualization var threshold0.4",
  "tables",
  "figure_novel_area_data_var.csv"
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

for (directory in c(
  figure_dir,
  table_dir
)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

method_order <- c(
  "rf_var",
  "mf_var"
)

method_labels <- c(
  rf_var = "Plain RF",
  mf_var = "Plain MF RF"
)

period_order <- c(
  "2011-2040",
  "2041-2070",
  "2071-2100"
)

ssp_order <- c(
  "SSP245",
  "SSP585"
)

rank_order <- c(
  "Top-1",
  "Top-3",
  "Top-5"
)

change_order <- c(
  "Stable ecotype niche",
  "Changed ecotype niche",
  "Novel ecosystem"
)

change_colors <- c(
  "Stable ecotype niche" = "#4C78A8",
  "Changed ecotype niche" = "#F2B134",
  "Novel ecosystem" = "#333333"
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


percent_labels <- function(x) {
  paste0(
    round(100 * x),
    "%"
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
    dpi = 300,
    bg = "white"
  )
}


theme_article <- function(base_size = 11) {
  theme_bw(
    base_size = base_size
  ) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(
        fill = "grey92",
        color = "grey35",
        linewidth = 0.45
      ),
      strip.text = element_text(
        face = "bold"
      ),
      plot.title = element_text(
        face = "bold",
        size = rel(1.15)
      ),
      plot.subtitle = element_text(
        size = rel(0.95),
        margin = margin(
          b = 8
        )
      ),
      plot.caption = element_text(
        hjust = 0,
        color = "grey30",
        size = rel(0.82),
        margin = margin(
          t = 7
        )
      ),
      legend.position = "top",
      legend.title = element_blank()
    )
}


# 2. Matched Top-k ecotype-niche change ========================================

matched <- fread(
  require_file(matched_file)
)

required_matched_columns <- c(
  "method",
  "method_label",
  "period",
  "ssp",
  "k",
  "rank_label",
  "change_code",
  "area_km2",
  "area_share"
)

missing_matched_columns <- setdiff(
  required_matched_columns,
  names(matched)
)

if (length(missing_matched_columns) > 0L) {
  stop(
    "Matched Top-k table is missing columns: ",
    paste(
      missing_matched_columns,
      collapse = ", "
    )
  )
}

matched[
  ,
  change_label := fcase(
    change_code == 1L, "Stable ecotype niche",
    change_code == 2L, "Changed ecotype niche",
    change_code == 3L, "Novel ecosystem"
  )
]

matched[
  ,
  `:=`(
    method = factor(
      method,
      levels = method_order
    ),
    period = factor(
      period,
      levels = period_order
    ),
    ssp = factor(
      ssp,
      levels = ssp_order
    ),
    rank_label = factor(
      rank_label,
      levels = rank_order
    ),
    change_label = factor(
      change_label,
      levels = change_order
    )
  )
]

fwrite(
  matched,
  file.path(
    table_dir,
    "Figure_var_7_matched_topk_ecotype_niche_change.csv"
  )
)

for (method_value in method_order) {
  plot_data <- matched[
    method == method_value
  ]
  
  method_title <- unname(
    method_labels[method_value]
  )
  
  matched_plot <- ggplot(
    plot_data,
    aes(
      x = period,
      y = area_share,
      fill = change_label
    )
  ) +
    geom_col(
      width = 0.64,
      color = "white",
      linewidth = 0.25
    ) +
    facet_grid(
      rows = vars(rank_label),
      cols = vars(ssp),
      switch = "y"
    ) +
    scale_fill_manual(
      values = change_colors,
      drop = FALSE
    ) +
    scale_y_continuous(
      breaks = seq(
        0,
        1,
        by = 0.2
      ),
      labels = percent_labels,
      expand = expansion(
        mult = c(
          0,
          0
        )
      )
    ) +
    coord_cartesian(
      ylim = c(
        0,
        1
      ),
      expand = FALSE
    ) +
    labs(
      title = paste0(
        "Projected changes in ecotype-niche assignment | ",
        method_title
      ),
      subtitle = paste0(
        "Reference and future maps use the same Top-k assignment rule; ",
        "areas are cell-size weighted."
      ),
      x = "Future period",
      y = "Mapped area (%)",
      fill = NULL,
      caption = paste0(
        "Top-3 and Top-5 restore the observed ecotype when it occurs within ",
        "the first k ranks. Novel ecosystem is Zone 99 and is independent of k."
      )
    ) +
    theme_article(
      base_size = 11
    ) +
    theme(
      strip.placement = "outside",
      axis.text.x = element_text(
        angle = 25,
        hjust = 1
      )
    )
  
  save_plot(
    matched_plot,
    file.path(
      figure_dir,
      paste0(
        "Figure_var_7_matched_topk_ecotype_niche_change_",
        method_value,
        ".png"
      )
    ),
    width = 10.5,
    height = 8.2
  )
}


# 3. Actual novel-ecosystem area ===============================================

novel <- fread(
  require_file(existing_novel_file)
)

required_novel_columns <- c(
  "method",
  "period",
  "ssp",
  "area_km2"
)

missing_novel_columns <- setdiff(
  required_novel_columns,
  names(novel)
)

if (length(missing_novel_columns) > 0L) {
  stop(
    "Novel-ecosystem table is missing columns: ",
    paste(
      missing_novel_columns,
      collapse = ", "
    )
  )
}

novel[
  ,
  `:=`(
    method_label = unname(
      method_labels[method]
    ),
    period = factor(
      period,
      levels = period_order
    ),
    ssp = factor(
      ssp,
      levels = ssp_order
    ),
    area_million_km2 = area_km2 / 1e6,
    result_class = "Novel ecosystem"
  )
]

novel <- unique(
  novel[
    ,
    .(
      method,
      method_label,
      period,
      ssp,
      area_km2,
      area_million_km2,
      result_class
    )
  ]
)

stopifnot(
  nrow(novel) ==
    length(method_order) *
    length(period_order) *
    length(ssp_order)
)

fwrite(
  novel,
  file.path(
    table_dir,
    "Figure_var_6_novel_ecosystem_area.csv"
  )
)

novel_plot <- ggplot(
  novel,
  aes(
    x = period,
    y = area_million_km2,
    group = method_label,
    linetype = method_label,
    shape = method_label
  )
) +
  geom_line(
    color = "black",
    linewidth = 0.7
  ) +
  geom_point(
    color = "black",
    fill = "white",
    size = 2.2,
    stroke = 0.65
  ) +
  facet_wrap(
    vars(ssp),
    nrow = 1
  ) +
  scale_linetype_manual(
    values = c(
      "Plain RF" = "solid",
      "Plain MF RF" = "dashed"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Plain RF" = 21,
      "Plain MF RF" = 24
    )
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.04,
        0.08
      )
    )
  ) +
  labs(
    title = "Projected novel ecosystem area",
    subtitle = paste0(
      "All 53 current ecotypes have dual suitability < 0.4; ",
      "this definition is independent of Top-k."
    ),
    x = "Future period",
    y = expression(
      "Novel ecosystem area (million " * km^2 * ")"
    ),
    linetype = NULL,
    shape = NULL
  ) +
  theme_article(
    base_size = 11
  ) +
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    )
  )

save_plot(
  novel_plot,
  file.path(
    figure_dir,
    "Figure_var_6_novel_ecosystem_area.png"
  ),
  width = 9.4,
  height = 4.8
)


cat(
  "Matched Top-k and novel-ecosystem figures completed.\n",
  "No raster or model analysis was rerun.\n",
  "Use Figure_var_7_matched_topk_ecotype_niche_change_*.png for Top-k.\n",
  "Use Figure_var_6_novel_ecosystem_area.png for Zone 99.\n",
  sep = ""
)
