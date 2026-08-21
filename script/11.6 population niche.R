# 11.6 Main-text population-niche change map
# ==============================================================================
# Purpose
# -------
# Create one reader-facing spatial figure for the source-ecotype population
# niches of Chamaecyparis formosensis.  The script reads existing Plain-RF dual
# suitability rasters only; it does not refit models or rerun projections.
#
# Each panel compares 2011--2040 with 2071--2100 under SSP585 for one source
# ecotype.  Cells are classified as retained suitable area, loss of suitability,
# or gain of suitability at the final dual-suitability threshold of 0.4.
#
# Output
# ------
# visualization var threshold0.4/figures/
#   Figure_var_10e_population_niche_change_maps_chamFor.png
# visualization var threshold0.4/tables/
#   Figure_var_10e_population_niche_change_maps_chamFor.csv
# ==============================================================================

library(terra)
library(data.table)
library(ggplot2)

rm(list = ls())
gc()


# 0. Paths and settings =========================================================

base_dir <- "H:/Jing/ecoChina2"

population_table_file <- file.path(
  base_dir,
  "visualization var threshold0.4",
  "tables",
  "Figure_var_10a_population_detail_long.csv"
)

figure_dir <- file.path(
  base_dir,
  "visualization var threshold0.4",
  "figures"
)

table_dir <- file.path(
  base_dir,
  "visualization var threshold0.4",
  "tables"
)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

species_code <- "chamFor"
method_value <- "rf_var"
early_scenario <- "2011-2040SSP585"
late_scenario <- "2071-2100SSP585"
dual_threshold <- 0.4

# Raster aggregation is used only for plotting.  Exact areas are calculated
# from the full-resolution categorical change raster.
target_plot_cells_per_population <- 300000L

output_figure <- file.path(
  figure_dir,
  "Figure_var_10e_population_niche_change_maps_chamFor.png"
)

output_table <- file.path(
  table_dir,
  "Figure_var_10e_population_niche_change_maps_chamFor.csv"
)


# 1. Helpers ====================================================================

require_file <- function(file) {
  if (!file.exists(file)) {
    stop("Missing required file: ", file)
  }
  file
}

short_zone_label <- function(zone, zone_name) {
  wrapped_name <- paste(
    strwrap(
      as.character(zone_name),
      width = 38
    ),
    collapse = "\n"
  )
  
  paste0(
    "Zone ",
    as.integer(zone),
    "\n",
    wrapped_name
  )
}

resolve_dual_file <- function(
    population_table,
    source_zone_value,
    scenario_value) {
  candidates <- unique(
    population_table[
      source_zone == as.integer(source_zone_value) &
        scenario == as.character(scenario_value),
      source_dual_raster
    ]
  )
  
  candidates <- candidates[
    !is.na(candidates) &
      nzchar(candidates)
  ]
  
  existing_candidates <- candidates[
    file.exists(candidates)
  ]
  
  if (length(existing_candidates) > 0L) {
    return(existing_candidates[[1L]])
  }
  
  fallback_candidates <- c(
    file.path(
      base_dir,
      "dual suit",
      method_value,
      scenario_value,
      paste0("dual_suitability_zone", source_zone_value, ".tif")
    ),
    file.path(
      base_dir,
      "dual suitability",
      method_value,
      scenario_value,
      paste0("dual_suitability_zone", source_zone_value, ".tif")
    )
  )
  
  existing_fallback <- fallback_candidates[
    file.exists(fallback_candidates)
  ]
  
  if (length(existing_fallback) == 0L) {
    stop(
      "No dual-suitability raster was found for source zone ",
      source_zone_value,
      " and scenario ",
      scenario_value,
      ". Checked the source_dual_raster column and: ",
      paste(fallback_candidates, collapse = "; ")
    )
  }
  
  existing_fallback[[1L]]
}

build_change_raster <- function(early_file, late_file) {
  early <- rast(require_file(early_file))[[1]]
  late <- rast(require_file(late_file))[[1]]
  
  if (!compareGeom(early, late, stopOnError = FALSE)) {
    stop(
      "Early- and late-century rasters have different geometry: ",
      early_file,
      " versus ",
      late_file
    )
  }
  
  change <- ifel(
    is.na(early) | is.na(late),
    NA,
    ifel(
      early >= dual_threshold & late >= dual_threshold,
      1L,
      ifel(
        early >= dual_threshold & late < dual_threshold,
        2L,
        ifel(
          early < dual_threshold & late >= dual_threshold,
          3L,
          0L
        )
      )
    )
  )
  
  names(change) <- "change_code"
  change
}

summarize_exact_area <- function(change, source_zone, zone_name) {
  cell_area <- cellSize(change, unit = "km")
  area_table <- as.data.table(
    zonal(
      cell_area,
      change,
      fun = "sum",
      na.rm = TRUE
    )
  )
  
  setnames(
    area_table,
    names(area_table),
    c("change_code", "area_km2")
  )
  
  complete_codes <- data.table(change_code = 0:3)
  area_table <- merge(
    complete_codes,
    area_table,
    by = "change_code",
    all.x = TRUE
  )
  area_table[is.na(area_km2), area_km2 := 0]
  
  area_table[
    ,
    `:=`(
      source_zone = as.integer(source_zone),
      zone_name = as.character(zone_name),
      Species = species_code,
      method = method_value,
      early_scenario = early_scenario,
      late_scenario = late_scenario,
      dual_threshold = dual_threshold,
      change_class = factor(
        change_code,
        levels = 0:3,
        labels = c(
          "Below threshold in both periods",
          "Retained suitable area",
          "Loss of suitability",
          "Gain of suitability"
        )
      )
    )
  ]
  
  area_table[]
}

prepare_plot_data <- function(change, source_zone, panel_label) {
  aggregation_factor <- max(
    1L,
    as.integer(
      ceiling(
        sqrt(
          ncell(change) /
            target_plot_cells_per_population
        )
      )
    )
  )
  
  plot_raster <- if (aggregation_factor > 1L) {
    aggregate(
      change,
      fact = aggregation_factor,
      fun = "modal",
      na.rm = TRUE
    )
  } else {
    change
  }
  
  plot_table <- as.data.table(
    as.data.frame(
      plot_raster,
      xy = TRUE,
      na.rm = TRUE
    )
  )
  
  plot_table[
    ,
    `:=`(
      source_zone = as.integer(source_zone),
      panel_label = as.character(panel_label),
      change_class = factor(
        change_code,
        levels = 0:3,
        labels = c(
          "Below threshold in both periods",
          "Retained suitable area",
          "Loss of suitability",
          "Gain of suitability"
        )
      )
    )
  ]
  
  plot_table[]
}


# 2. Population selection =======================================================

population_table <- fread(
  require_file(population_table_file)
)

required_columns <- c(
  "Species",
  "source_zone",
  "zone_name",
  "method",
  "scenario",
  "source_dual_raster"
)

missing_columns <- setdiff(
  required_columns,
  names(population_table)
)

if (length(missing_columns) > 0L) {
  stop(
    "Population table is missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

population_table <- population_table[
  Species == species_code &
    method == method_value
]

if (nrow(population_table) == 0L) {
  stop(
    "No ",
    species_code,
    " rows were found for ",
    method_value,
    "."
  )
}

population_lookup <- unique(
  population_table[
    ,
    .(
      source_zone = as.integer(source_zone),
      zone_name = as.character(zone_name)
    )
  ]
)

setorder(population_lookup, source_zone)


# 3. Build exact change classes and plotting data ===============================

plot_list <- vector("list", nrow(population_lookup))
area_list <- vector("list", nrow(population_lookup))

for (row_index in seq_len(nrow(population_lookup))) {
  source_zone_value <- population_lookup$source_zone[[row_index]]
  zone_name_value <- population_lookup$zone_name[[row_index]]
  
  early_file <- resolve_dual_file(
    population_table,
    source_zone_value,
    early_scenario
  )
  
  late_file <- resolve_dual_file(
    population_table,
    source_zone_value,
    late_scenario
  )
  
  change <- build_change_raster(
    early_file,
    late_file
  )
  
  panel_label <- short_zone_label(
    source_zone_value,
    zone_name_value
  )
  
  area_list[[row_index]] <- summarize_exact_area(
    change,
    source_zone_value,
    zone_name_value
  )
  
  plot_list[[row_index]] <- prepare_plot_data(
    change,
    source_zone_value,
    panel_label
  )
  
  rm(change)
  gc()
}

plot_data <- rbindlist(plot_list, use.names = TRUE)
area_summary <- rbindlist(area_list, use.names = TRUE)

panel_order <- vapply(
  seq_len(nrow(population_lookup)),
  function(row_index) {
    short_zone_label(
      population_lookup$source_zone[[row_index]],
      population_lookup$zone_name[[row_index]]
    )
  },
  character(1)
)

plot_data[
  ,
  panel_label := factor(
    panel_label,
    levels = panel_order
  )
]

fwrite(
  area_summary,
  output_table
)


# 4. Reader-facing figure =======================================================

change_colors <- c(
  "Below threshold in both periods" = "#F1F1F1",
  "Retained suitable area" = "#4C78A8",
  "Loss of suitability" = "#D95F02",
  "Gain of suitability" = "#1B9E77"
)

population_map <- ggplot(
  plot_data,
  aes(
    x = x,
    y = y,
    fill = change_class
  )
) +
  geom_raster() +
  facet_wrap(
    ~ panel_label,
    ncol = 2,
    drop = FALSE
  ) +
  coord_equal(expand = FALSE) +
  scale_fill_manual(
    values = change_colors,
    breaks = c(
      "Retained suitable area",
      "Loss of suitability",
      "Gain of suitability"
    ),
    drop = FALSE
  ) +
  labs(
    title = expression(
      italic("Chamaecyparis formosensis") *
        " source-ecotype population niches"
    ),
    subtitle = paste0(
      "Change in dual suitability >= ",
      dual_threshold,
      " from ",
      early_scenario,
      " to ",
      late_scenario,
      " | Plain RF"
    ),
    fill = NULL
  ) +
  theme_void(base_size = 10.5) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0,
      size = 14
    ),
    plot.subtitle = element_text(
      hjust = 0,
      size = 10
    ),
    strip.text = element_text(
      face = "bold",
      size = 8.5,
      lineheight = 0.95
    ),
    legend.position = "bottom",
    legend.key.width = grid::unit(1.1, "cm"),
    panel.spacing = grid::unit(0.35, "lines"),
    plot.margin = margin(8, 10, 8, 10)
  )

ggsave(
  output_figure,
  population_map,
  width = 10.5,
  height = 13.0,
  dpi = 320,
  bg = "white"
)

cat(
  "COMPLETE\n",
  "Figure: ", output_figure, "\n",
  "Area table: ", output_table, "\n",
  "No model fitting or projection was rerun.\n",
  sep = ""
)
