# 5.4 Reader-facing feature-importance outputs
# ==============================================================================
# Run after script 5.3. This script reads the existing FI CSV files only. It
# does not load RF objects, extract importance again, or recalculate any raster
# quantity.
#
# The complete 53-ecotype audit remains unchanged in the original tables. The
# reader-facing tables and figures exclude the "no vegetation" class because it
# is not used in the ecological interpretation of vegetation-category niches.
#
# Outputs
# -------
# assessment_var/feature_importance_analysis/
#   reader_tables/
#   reader_figures/
#   reader_figures/category representative zones/
# ==============================================================================

library(data.table)
library(ggplot2)

rm(list = ls())
gc()


# 0. Paths and settings =========================================================

base_dir <- "H:/Jing/ecoChina2"

feature_importance_root <- file.path(
  base_dir,
  "assessment_var",
  "feature_importance_analysis"
)

source_table_dir <- file.path(
  feature_importance_root,
  "tables"
)

reader_table_dir <- file.path(
  feature_importance_root,
  "reader_tables"
)

reader_figure_dir <- file.path(
  feature_importance_root,
  "reader_figures"
)

representative_figure_dir <- file.path(
  reader_figure_dir,
  "category representative zones"
)

for (directory in c(
  reader_table_dir,
  reader_figure_dir,
  representative_figure_dir
)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

top_n_category_variables <- 6L
top_n_zone_variables <- 10L

niche_order <- c(
  "climate",
  "soil"
)

climate_period_levels <- c(
  "Annual",
  "Winter",
  "Spring",
  "Summer",
  "Autumn"
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


read_fi_table <- function(filename) {
  fread(
    require_file(
      file.path(
        source_table_dir,
        filename
      )
    )
  )
}


is_no_vegetation <- function(x) {
  normalized <- tolower(
    gsub(
      "[ -]+",
      "_",
      trimws(
        as.character(x)
      )
    )
  )
  
  normalized == "no_vegetation"
}


filter_reader_categories <- function(table) {
  result <- copy(table)
  
  if ("category1" %in% names(result)) {
    result <- result[
      is.na(category1) |
        !is_no_vegetation(category1)
    ]
  }
  
  if ("category2" %in% names(result)) {
    result <- result[
      is.na(category2) |
        !is_no_vegetation(category2)
    ]
  }
  
  result[]
}


percent_label <- function(x) {
  paste0(
    round(100 * x),
    "%"
  )
}


safe_name <- function(x) {
  x <- gsub(
    "[^A-Za-z0-9_]+",
    "_",
    x
  )
  
  gsub(
    "^_+|_+$",
    "",
    x
  )
}


wrap_text <- function(
    x,
    width = 42L) {
  vapply(
    x,
    function(value) {
      paste(
        strwrap(
          value,
          width = width
        ),
        collapse = "\n"
      )
    },
    character(1)
  )
}


category_roles_only <- function(x) {
  vapply(
    strsplit(
      as.character(x),
      "; ",
      fixed = TRUE
    ),
    function(roles) {
      roles <- roles[
        grepl(
          "^category ",
          roles
        )
      ]
      
      paste(
        gsub(
          "^category ",
          "",
          roles
        ),
        collapse = "; "
      )
    },
    character(1)
  )
}


niche_title <- function(niche) {
  if (
    identical(
      niche,
      "climate"
    )
  ) {
    "Climatic"
  } else {
    "Topsoil"
  }
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


theme_article <- function(base_size = 9.5) {
  theme_bw(
    base_size = base_size
  ) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      strip.text = element_text(
        face = "bold"
      ),
      plot.title = element_text(
        face = "bold"
      )
    )
}


# 2. Read and filter the existing summaries ====================================

method_agreement <- filter_reader_categories(
  read_fi_table(
    "FI_04_rf_vs_mf_importance_agreement.csv"
  )
)

category_variables <- filter_reader_categories(
  read_fi_table(
    "FI_07_category_variable_summary.csv"
  )
)

category_factor_groups <- filter_reader_categories(
  read_fi_table(
    "FI_10_category_ecological_factor_groups.csv"
  )
)

category_climate_period <- filter_reader_categories(
  read_fi_table(
    "FI_11_category_climate_period_summary.csv"
  )
)

representative_zones <- filter_reader_categories(
  read_fi_table(
    "FI_14_representative_and_area_extreme_zones.csv"
  )
)

selected_zone_factors <- filter_reader_categories(
  read_fi_table(
    "FI_15_selected_zone_top_factors.csv"
  )
)

area_concentration <- filter_reader_categories(
  read_fi_table(
    "FI_16_zone_area_and_importance_concentration.csv"
  )
)


# 3. Reader-facing tables =======================================================

reader_source_files <- c(
  "FI_04_rf_vs_mf_importance_agreement.csv",
  "FI_07_category_variable_summary.csv",
  "FI_08_category_common_and_top_factors.csv",
  "FI_09_category_contrasts_vs_all_zones.csv",
  "FI_10_category_ecological_factor_groups.csv",
  "FI_11_category_climate_period_summary.csv",
  "FI_12_theory_review_shortlist.csv",
  "FI_13_category_importance_profile_similarity.csv",
  "FI_14_representative_and_area_extreme_zones.csv",
  "FI_15_selected_zone_top_factors.csv",
  "FI_16_zone_area_and_importance_concentration.csv",
  "FI_17_area_importance_concentration_correlations.csv"
)

for (filename in reader_source_files) {
  reader_table <- filter_reader_categories(
    read_fi_table(filename)
  )
  
  fwrite(
    reader_table,
    file.path(
      reader_table_dir,
      filename
    )
  )
}


# 4. Category-level climatic and topsoil variables =============================

category_top_variables <- category_variables[
  category_variable_rank <=
    top_n_category_variables
]

category_top_variables[
  ,
  `:=`(
    lower_share = pmax(
      0,
      mean_importance_share -
        se_importance_share
    ),
    upper_share =
      mean_importance_share +
      se_importance_share
  )
]

for (niche_value in niche_order) {
  figure_data <- copy(
    category_top_variables[
      niche == niche_value
    ]
  )
  
  setorder(
    figure_data,
    category_label,
    mean_importance_share,
    feature
  )
  
  figure_data[
    ,
    feature_panel := factor(
      paste(
        category2,
        feature,
        sep = "|||"
      ),
      levels = unique(
        paste(
          category2,
          feature,
          sep = "|||"
        )
      )
    )
  ]
  
  variable_figure <- ggplot(
    figure_data,
    aes(
      x = mean_importance_share,
      y = feature_panel
    )
  ) +
    geom_errorbarh(
      aes(
        xmin = lower_share,
        xmax = upper_share
      ),
      height = 0.16,
      linewidth = 0.42,
      color = "grey45"
    ) +
    geom_point(
      aes(
        size = selection_frequency
      ),
      shape = 21,
      fill = "white",
      color = "black",
      stroke = 0.45
    ) +
    facet_wrap(
      vars(category_label),
      scales = "free_y",
      ncol = 2
    ) +
    scale_y_discrete(
      labels = function(x) {
        sub(
          "^.*\\|\\|\\|",
          "",
          x
        )
      }
    ) +
    scale_x_continuous(
      labels = percent_label,
      expand = expansion(
        mult = c(
          0,
          0.08
        )
      )
    ) +
    scale_size_continuous(
      range = c(
        1.7,
        3.8
      ),
      limits = c(
        0,
        1
      ),
      labels = percent_label
    ) +
    labs(
      x = "Mean normalized permutation importance",
      y = NULL,
      size = "Ecotypes selecting\nthe variable",
      title = paste0(
        niche_title(niche_value),
        " variables contributing to ecotype-niche models"
      ),
      subtitle =
        "Category means across ecotypes; error bars show +/-1 SE."
    ) +
    theme_article(
      base_size = 9.5
    ) +
    theme(
      legend.position = "top",
      panel.spacing = grid::unit(
        0.75,
        "lines"
      )
    )
  
  save_plot(
    variable_figure,
    file.path(
      reader_figure_dir,
      paste0(
        "Figure_FI_reader_category_variables_",
        niche_value,
        ".png"
      )
    ),
    width = 10.5,
    height = 10.0
  )
}


# 5. Ecological factor groups by vegetation category ===========================

for (niche_value in niche_order) {
  group_data <- copy(
    category_factor_groups[
      niche == niche_value
    ]
  )
  
  group_levels <- sort(
    unique(
      group_data$feature_group
    )
  )
  
  group_colors <- setNames(
    hcl.colors(
      length(group_levels),
      palette = "Dark 3"
    ),
    group_levels
  )
  
  group_data[
    ,
    feature_group := factor(
      feature_group,
      levels = group_levels
    )
  ]
  
  group_figure <- ggplot(
    group_data,
    aes(
      x = category_label,
      y = mean_group_importance_share,
      fill = feature_group
    )
  ) +
    geom_col(
      width = 0.68,
      color = "white",
      linewidth = 0.15
    ) +
    scale_fill_manual(
      values = group_colors
    ) +
    scale_y_continuous(
      breaks = seq(
        0,
        1,
        by = 0.2
      ),
      labels = percent_label,
      expand = expansion(
        mult = c(
          0,
          0
        )
      )
    ) +
    coord_flip(
      ylim = c(
        0,
        1
      ),
      expand = FALSE
    ) +
    labs(
      x = NULL,
      y = "Mean share of normalized importance",
      fill = NULL,
      title = paste0(
        niche_title(niche_value),
        " feature groups across vegetation categories"
      ),
      subtitle =
        "RF and MF RF importance shares were averaged within each ecotype before category summaries."
    ) +
    theme_article(
      base_size = 10
    ) +
    theme(
      legend.position = "right"
    )
  
  save_plot(
    group_figure,
    file.path(
      reader_figure_dir,
      paste0(
        "Figure_FI_reader_category_factor_groups_",
        niche_value,
        ".png"
      )
    ),
    width = 10.5,
    height = 6.2
  )
}


# 6. Annual and seasonal climate contribution ==================================

climate_period_data <- copy(
  category_climate_period
)

climate_period_data[
  ,
  climate_period := factor(
    climate_period,
    levels = climate_period_levels
  )
]

climate_period_colors <- c(
  Annual = "#4D4D4D",
  Winter = "#5B8DB8",
  Spring = "#78A96B",
  Summer = "#D58A45",
  Autumn = "#9A6FB0"
)

climate_period_figure <- ggplot(
  climate_period_data,
  aes(
    x = category_label,
    y = mean_period_importance_share,
    fill = climate_period
  )
) +
  geom_col(
    width = 0.68,
    color = "white",
    linewidth = 0.15
  ) +
  scale_fill_manual(
    values = climate_period_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    breaks = seq(
      0,
      1,
      by = 0.2
    ),
    labels = percent_label,
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  coord_flip(
    ylim = c(
      0,
      1
    ),
    expand = FALSE
  ) +
  labs(
    x = NULL,
    y = "Mean share of normalized climate importance",
    fill = NULL,
    title = "Annual and seasonal contributions to climatic niche models",
    subtitle =
      "Season refers to the ClimateNA variable period."
  ) +
  theme_article(
    base_size = 10
  ) +
  theme(
    legend.position = "top"
  )

save_plot(
  climate_period_figure,
  file.path(
    reader_figure_dir,
    "Figure_FI_reader_climate_period_composition.png"
  ),
  width = 10.0,
  height = 5.8
)


# 7. Agreement between Plain RF and Plain MF RF ================================

method_agreement_long <- rbindlist(
  list(
    method_agreement[
      ,
      .(
        zoneID,
        zone_name,
        category_label,
        niche,
        agreement_metric = "Spearman correlation",
        agreement_value = importance_spearman
      )
    ],
    method_agreement[
      ,
      .(
        zoneID,
        zone_name,
        category_label,
        niche,
        agreement_metric = "Top-5 overlap",
        agreement_value = top5_jaccard
      )
    ]
  ),
  use.names = TRUE,
  fill = TRUE
)

method_agreement_figure <- ggplot(
  method_agreement_long,
  aes(
    x = category_label,
    y = agreement_value
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.5
  ) +
  geom_point(
    position = position_jitter(
      width = 0.12,
      height = 0,
      seed = 49
    ),
    shape = 21,
    fill = "grey75",
    color = "black",
    size = 1.5,
    stroke = 0.35,
    alpha = 0.72
  ) +
  facet_grid(
    agreement_metric ~ niche
  ) +
  scale_y_continuous(
    breaks = seq(
      0,
      1,
      by = 0.2
    )
  ) +
  coord_cartesian(
    ylim = c(
      0,
      1
    )
  ) +
  labs(
    x = NULL,
    y = "Agreement",
    title = "Agreement of feature importance between Plain RF and Plain MF RF",
    subtitle = "Each point represents one ecotype model."
  ) +
  theme_article(
    base_size = 9
  ) +
  theme(
    axis.text.x = element_text(
      angle = 40,
      hjust = 1
    )
  )

save_plot(
  method_agreement_figure,
  file.path(
    reader_figure_dir,
    "Figure_FI_reader_rf_vs_mf_importance_agreement.png"
  ),
  width = 12.0,
  height = 7.0
)


# 8. Observed ecotype area and importance concentration ========================

category_values <- sort(
  unique(
    area_concentration$category2
  )
)

category_colors <- setNames(
  hcl.colors(
    length(category_values),
    palette = "Dark 3"
  ),
  category_values
)

area_concentration_figure <- ggplot(
  area_concentration,
  aes(
    x = reference_area_km2,
    y = effective_number_of_features,
    color = category2,
    shape = selected_for_detail
  )
) +
  geom_point(
    size = 2.0,
    alpha = 0.78,
    stroke = 0.42
  ) +
  facet_wrap(
    vars(niche),
    nrow = 1,
    scales = "free_y"
  ) +
  scale_x_log10() +
  scale_color_manual(
    values = category_colors
  ) +
  scale_shape_manual(
    values = c(
      `FALSE` = 16,
      `TRUE` = 17
    )
  ) +
  labs(
    x = expression(
      "Observed reference area (" * km^2 * ", log scale)"
    ),
    y = "Effective number of important variables",
    color = "Vegetation category",
    shape = "Selected for\ndetailed comparison",
    title = "Observed ecotype area and feature-importance concentration",
    subtitle =
      "The effective number is exp(Shannon entropy) of normalized positive permutation importance."
  ) +
  theme_article(
    base_size = 9.5
  ) +
  theme(
    legend.position = "right"
  )

save_plot(
  area_concentration_figure,
  file.path(
    reader_figure_dir,
    "Figure_FI_reader_area_vs_importance_concentration.png"
  ),
  width = 11.2,
  height = 6.0
)


# 9. Representative and observed-area-extreme ecotypes =========================

category_role_pattern <- paste(
  c(
    "category profile medoid",
    "category largest observed area",
    "category smallest observed area"
  ),
  collapse = "|"
)

representative_zones <- representative_zones[
  !is.na(selection_roles) &
    grepl(
      category_role_pattern,
      selection_roles
    )
]

representative_figure_index <- list()

for (category_value in sort(
  unique(
    representative_zones$category2
  )
)) {
  selected_zones <- representative_zones[
    category2 == category_value
  ]
  
  figure_data <- copy(
    selected_zone_factors[
      category2 == category_value &
        zoneID %in% selected_zones$zoneID &
        consensus_rank <= top_n_zone_variables
    ]
  )
  
  if (nrow(figure_data) == 0L) {
    next
  }
  
  figure_data <- merge(
    figure_data,
    selected_zones[
      ,
      .(
        zoneID,
        reader_selection_roles = selection_roles
      )
    ],
    by = "zoneID",
    all.x = TRUE,
    sort = FALSE
  )
  
  figure_data[
    ,
    role_short := category_roles_only(
      reader_selection_roles
    )
  ]
  
  figure_data[
    ,
    panel_label := paste0(
      "Zone ",
      zoneID,
      " | ",
      niche,
      "\n",
      wrap_text(
        zone_name,
        width = 38L
      ),
      "\n",
      wrap_text(
        role_short,
        width = 48L
      )
    )
  ]
  
  setorder(
    figure_data,
    panel_label,
    consensus_importance_share,
    feature
  )
  
  figure_data[
    ,
    feature_panel := factor(
      paste(
        panel_label,
        feature,
        sep = "|||"
      ),
      levels = unique(
        paste(
          panel_label,
          feature,
          sep = "|||"
        )
      )
    )
  ]
  
  representative_figure <- ggplot(
    figure_data,
    aes(
      x = consensus_importance_share,
      y = feature_panel
    )
  ) +
    geom_col(
      width = 0.68,
      fill = "grey35"
    ) +
    facet_wrap(
      vars(panel_label),
      scales = "free_y",
      ncol = 2
    ) +
    scale_y_discrete(
      labels = function(x) {
        sub(
          "^.*\\|\\|\\|",
          "",
          x
        )
      }
    ) +
    scale_x_continuous(
      labels = percent_label,
      expand = expansion(
        mult = c(
          0,
          0.06
        )
      )
    ) +
    labs(
      x = "Consensus normalized permutation importance",
      y = NULL,
      title = paste0(
        "Representative and observed-area-extreme ecotypes: ",
        gsub(
          "_",
          " ",
          category_value,
          fixed = TRUE
        )
      ),
      subtitle = paste0(
        "The profile medoid summarizes the category importance profile; ",
        "area extremes use the reference vegetation raster."
      )
    ) +
    theme_article(
      base_size = 8.8
    ) +
    theme(
      strip.text = element_text(
        face = "bold",
        size = 7.6
      ),
      panel.spacing = grid::unit(
        0.8,
        "lines"
      )
    )
  
  figure_file <- file.path(
    representative_figure_dir,
    paste0(
      "Figure_FI_reader_category_zones_",
      safe_name(category_value),
      ".png"
    )
  )
  
  save_plot(
    representative_figure,
    figure_file,
    width = 12.5,
    height = 10.0
  )
  
  representative_figure_index[[
    length(representative_figure_index) + 1L
  ]] <- data.table(
    category2 = category_value,
    n_unique_zones = uniqueN(
      figure_data$zoneID
    ),
    figure_file = figure_file
  )
}

representative_figure_index <- rbindlist(
  representative_figure_index,
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  representative_figure_index,
  file.path(
    reader_table_dir,
    "FI_18_reader_representative_figure_index.csv"
  )
)


# 10. Final checks and report ===================================================

reader_csv_files <- list.files(
  reader_table_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

for (reader_file in reader_csv_files) {
  reader_table <- fread(
    reader_file
  )
  
  if (
    "category1" %in% names(reader_table) &&
    any(
      is_no_vegetation(
        reader_table$category1
      ),
      na.rm = TRUE
    )
  ) {
    stop(
      "Excluded category remains in ",
      reader_file
    )
  }
  
  if (
    "category2" %in% names(reader_table) &&
    any(
      is_no_vegetation(
        reader_table$category2
      ),
      na.rm = TRUE
    )
  ) {
    stop(
      "Excluded category remains in ",
      reader_file
    )
  }
}

cat(
  "Reader-facing feature-importance outputs completed.\n",
  "No RF model or raster analysis was rerun.\n",
  "Complete audit retained in: ",
  source_table_dir,
  "\n",
  "Filtered reader tables: ",
  reader_table_dir,
  "\n",
  "Filtered reader figures: ",
  reader_figure_dir,
  "\n",
  sep = ""
)
