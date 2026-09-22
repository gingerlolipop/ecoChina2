# Centralized manuscript visualization for the clean Multi-Forest workflow
# ==============================================================================
# Run after scripts 5, 5.3, 6, 7 and 8. Analysis tables are never recalculated
# here; this script turns their canonical outputs into manuscript-ready figures.
# Set ECOCHINA2_STACK_SCENARIOS to a comma-separated subset of scenarios when
# choosing the detailed Top-k raster-stack figures. Set ECOCHINA2_FIGURE_IDS
# to a comma-separated list of figure-family IDs to redraw only selected
# outputs (for example, importance_category_variables,topk_exploded_stack).
# Bivariate suitability figures are opt-in: set ECOCHINA2_BIVARIATE_SCENARIOS
# and, for zone-level pages, ECOCHINA2_BIVARIATE_ZONES. They read existing
# climate and soil suitability rasters and never refit or predict models.

library(data.table)
library(terra)
library(ggplot2)
library(patchwork)

rm(list = ls())
gc()

script_revision <- "2026-09-15-assessment-compat-v2"


# 0. Paths, constants and palette ==============================================

find_project_root <- function(path = getwd()) {
  configured <- Sys.getenv("ECOCHINA2_DIR", unset = "")
  if (nzchar(configured)) {
    return(normalizePath(configured, winslash = "/", mustWork = TRUE))
  }
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "script", "11. visualization.R"))) return(path)
    parent <- dirname(path)
    if (parent == path) stop("Run inside the repository or set ECOCHINA2_DIR.")
    path <- parent
  }
}

base_dir <- find_project_root()
figure_dir <- file.path(base_dir, "figures")
figure_subdirs <- file.path(
  figure_dir,
  c("assessment", "importance", "maps", "niche", "topk")
)
for (directory in c(figure_dir, figure_subdirs)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}
cat("[SCRIPT 11 REVISION]", script_revision, "\n")
cat("[PROJECT ROOT]", base_dir, "\n")
cat("[FIGURE ROOT]", normalizePath(
  figure_dir, winslash = "/", mustWork = TRUE
), "\n")

modeled_zones <- c(1:7, 9:50, 52:55)
scenarios <- c(
  "normal",
  "2011-2040SSP245", "2041-2070SSP245", "2071-2100SSP245",
  "2011-2040SSP585", "2041-2070SSP585", "2071-2100SSP585"
)
future_scenarios <- setdiff(scenarios, "normal")
k_values <- 1:5
dual_threshold <- 0.4
maximum_map_cells <- 250000L

parse_scenario_setting <- function(name, default) {
  setting <- trimws(Sys.getenv(name, unset = default))
  selected <- if (tolower(setting) == "all") {
    scenarios
  } else {
    trimws(strsplit(setting, ",", fixed = TRUE)[[1]])
  }
  selected <- selected[nzchar(selected)]
  unknown <- setdiff(selected, scenarios)
  if (length(unknown)) {
    stop("Unknown ", name, " value(s): ", paste(unknown, collapse = ", "))
  }
  unique(selected)
}

bivariate_scenarios <- parse_scenario_setting(
  "ECOCHINA2_BIVARIATE_SCENARIOS",
  ""
)
bivariate_zone_setting <- trimws(Sys.getenv(
  "ECOCHINA2_BIVARIATE_ZONES", unset = ""
))
bivariate_zones <- if (!nzchar(bivariate_zone_setting)) {
  integer()
} else if (tolower(bivariate_zone_setting) == "all") {
  modeled_zones
} else {
  suppressWarnings(as.integer(trimws(strsplit(
    bivariate_zone_setting, ",", fixed = TRUE
  )[[1]])))
}
if (anyNA(bivariate_zones) || length(setdiff(bivariate_zones, modeled_zones))) {
  stop(
    "ECOCHINA2_BIVARIATE_ZONES must contain modeled zone IDs or 'all'."
  )
}
bivariate_zones <- unique(bivariate_zones)

first_existing <- function(paths, label, required = TRUE) {
  hit <- paths[file.exists(paths)]
  if (length(hit)) return(hit[[1]])
  if (required) stop("Missing ", label, ":\n", paste(paths, collapse = "\n"))
  NA_character_
}

reference_file <- first_existing(
  c(
    file.path(base_dir, "raster", "ecosys_ori.tif"),
    file.path(base_dir, "data", "ecosys_ori.tif")
  ),
  "reference raster"
)
palette_file <- first_existing(
  c(
    file.path(base_dir, "color_palette_China.csv"),
    file.path(base_dir, "data", "zone_palette.csv")
  ),
  "zone palette"
)
if (!file.exists(palette_file)) {
  stop("Missing zone palette. Run script/color_palette.R first.")
}
zone_palette <- fread(palette_file)
required_palette <- c("zoneID", "zone", "category2", "COLOR")
if (!all(required_palette %in% names(zone_palette))) {
  stop("zone_palette.csv needs: ", paste(required_palette, collapse = ", "))
}
zone_palette[, zoneID := as.integer(zoneID)]
zone_colors <- setNames(zone_palette$COLOR, as.character(zone_palette$zoneID))
zone_labels <- setNames(
  paste0("Zone ", zone_palette$zoneID),
  as.character(zone_palette$zoneID)
)

theme_manuscript <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey94", colour = "grey75"),
    legend.key.height = grid::unit(0.35, "cm"),
    plot.title.position = "plot"
  )

scenario_label <- function(x) {
  fifelse(
    x == "normal",
    "Reference (1961–1990)",
    paste0(sub("SSP.*$", "", x), " · ", sub("^.*(SSP[0-9]+)$", "\\1", x))
  )
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_-]+", "_", x)
  gsub("^_+|_+$", "", x)
}

script6_safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  gsub("^_+|_+$", "", x)
}

script6_legacy_safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "unnamed")
}

required_file <- function(path, label = basename(path)) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path)
  path
}

files_current <- function(outputs, inputs) {
  outputs <- outputs[!is.na(outputs) & nzchar(outputs)]
  if (!length(outputs) || !all(file.exists(outputs))) return(FALSE)
  inputs <- inputs[!is.na(inputs) & nzchar(inputs) & file.exists(inputs)]
  if (!length(inputs)) return(TRUE)
  isTRUE(min(file.info(outputs)$mtime) >= max(file.info(inputs)$mtime))
}

first_current_raster <- function(paths, template_file, inputs, label) {
  template <- rast(template_file)
  valid <- vapply(paths, function(path) {
    if (!file.exists(path) || !files_current(path, inputs)) return(FALSE)
    tryCatch({
      candidate <- rast(path)
      nlyr(candidate) >= 1L &&
        compareGeom(candidate, template, stopOnError = FALSE)
    }, error = function(e) FALSE)
  }, logical(1))
  if (any(valid)) return(paths[which(valid)[1]])
  stop("Missing or stale ", label, "; rerun script/8. future top k.R.")
}

read_required <- function(path, label = basename(path)) {
  fread(required_file(path, label))
}

binary_map_file <- function(scenario) {
  first_existing(
    c(
      file.path(
        base_dir, "result maps", "mf_var",
        paste0(
          "assigned_zone_", scenario,
          "_threshold0.4_tol1e-04_novel99_maskNA8_noNovelNormal.tif"
        )
      ),
      file.path(
        base_dir, "result maps", "mf_var",
        paste0("assigned_zone_", scenario, ".tif")
      )
    ),
    paste("binary Multi-Forest map", scenario)
  )
}

multiclass_map_file <- function(scenario) {
  established <- if (scenario == "normal") {
    "assigned_zone_normal_multiclass_rf.tif"
  } else {
    paste0("assigned_zone_", scenario, "_multiclass_rf.tif")
  }
  first_existing(
    c(
      file.path(base_dir, "result maps", "multiclass_rf", established),
      file.path(base_dir, "result maps", "multiclass_rf", paste0("assigned_zone_", scenario, ".tif")),
      file.path(base_dir, "multiclass", "maps", paste0("assigned_zone_", scenario, ".tif"))
    ),
    paste("multiclass map", scenario)
  )
}

multiclass_table_file <- function(name) {
  first_existing(
    c(
      file.path(base_dir, "assessment", "multiclass_rf", name),
      file.path(base_dir, "multiclass", "tables", name)
    ),
    paste("multiclass table", name)
  )
}

niche_table_file <- function(kind) {
  legacy_name <- c(
    ecosystem = "future_ecosystem_area_var.csv",
    transition = "future_ecosystem_transition_var.csv",
    population = "dual_population_niche_area_var.csv",
    species = "dual_species_niche_area_var.csv"
  )[[kind]]
  current_name <- c(
    ecosystem = "ecosystem_area.csv",
    transition = "assigned_transitions.csv",
    population = "population_area.csv",
    species = "species_area.csv"
  )[[kind]]
  legacy_root <- if (kind %in% c("ecosystem", "transition")) {
    file.path(base_dir, "future tree niche var", "tables")
  } else {
    file.path(base_dir, "future tree niche dual suitability var", "tables")
  }
  first_existing(
    c(
      file.path(legacy_root, legacy_name),
      file.path(base_dir, "population_species", "tables", current_name)
    ),
    paste(kind, "niche table")
  )
}

keep_multiforest <- function(x) {
  if ("method" %in% names(x)) x <- x[method == "mf_var"]
  x
}

read_zone_size <- function() {
  file <- file.path(base_dir, "data", "zone_lookup.csv")
  if (file.exists(file)) return(fread(file))
  if (!"count" %in% names(zone_palette)) {
    stop("Neither data/zone_lookup.csv nor palette counts are available.")
  }
  zone_palette[, .(zoneID, n_cells = as.numeric(count))]
}

save_gg <- function(plot, filename, width = 9, height = 6, dpi = 320) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = filename, plot = plot, width = width, height = height,
    units = "in", dpi = dpi, bg = "white"
  )
  info <- file.info(filename)
  if (!file.exists(filename) || is.na(info$size) || info$size <= 0) {
    stop("Figure was not written correctly: ", filename)
  }
  saved <- normalizePath(filename, winslash = "/", mustWork = TRUE)
  cat("[SAVED]", saved, "\n")
  invisible(saved)
}

figure_log <- list()
allowed_sections <- c("assessment", "importance", "maps", "niche", "topk")
section_setting <- tolower(Sys.getenv("ECOCHINA2_FIGURE_SECTIONS", "all"))
selected_sections <- trimws(strsplit(section_setting, ",", fixed = TRUE)[[1]])
if (identical(selected_sections, "all")) selected_sections <- allowed_sections
unknown_sections <- setdiff(selected_sections, allowed_sections)
if (length(unknown_sections)) {
  stop("Unknown ECOCHINA2_FIGURE_SECTIONS: ", paste(unknown_sections, collapse = ", "))
}

figure_id_setting <- tolower(Sys.getenv("ECOCHINA2_FIGURE_IDS", "all"))
selected_figure_ids <- trimws(strsplit(figure_id_setting, ",", fixed = TRUE)[[1]])
selected_figure_ids <- selected_figure_ids[nzchar(selected_figure_ids)]
if (!length(selected_figure_ids)) selected_figure_ids <- "all"
if ("all" %in% selected_figure_ids) selected_figure_ids <- "all"

figure_section <- function(id) {
  if (grepl("^importance", id)) return("importance")
  if (grepl("^(species|population)", id)) return("niche")
  if (grepl("^topk", id)) return("topk")
  if (grepl("^(assessment|soil_gate)", id)) return("assessment")
  "maps"
}

run_figure <- function(id, expression) {
  section <- figure_section(id)
  if (!(section %in% selected_sections)) {
    figure_log[[length(figure_log) + 1L]] <<- data.table(
      figure = id, section = section, status = "skipped",
      message = "section not selected", seconds = 0
    )
    return(invisible(NULL))
  }
  if (!identical(selected_figure_ids, "all") && !(id %in% selected_figure_ids)) {
    figure_log[[length(figure_log) + 1L]] <<- data.table(
      figure = id, section = section, status = "skipped",
      message = "figure ID not selected", seconds = 0
    )
    return(invisible(NULL))
  }
  started <- Sys.time()
  result <- tryCatch({
    value <- force(expression)
    if (inherits(value, "figure_skip")) {
      data.table(
        figure = id, section = section, status = "skipped",
        message = as.character(value$message),
        seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
      )
    } else {
      data.table(
        figure = id, section = section, status = "complete", message = NA_character_,
        seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
      )
    }
  }, figure_input_missing = function(e) {
    cat("[MISSING INPUT]", id, "|", conditionMessage(e), "\n")
    data.table(
      figure = id, section = section, status = "missing_input",
      message = conditionMessage(e),
      seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
    )
  }, error = function(e) {
    warning("[", id, "] ", conditionMessage(e), call. = FALSE)
    data.table(
      figure = id, section = section, status = "failed", message = conditionMessage(e),
      seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
    )
  })
  figure_log[[length(figure_log) + 1L]] <<- result
  invisible(result)
}

skip_figure <- function(message) {
  structure(list(message = message), class = "figure_skip")
}

missing_figure_input <- function(message) {
  stop(structure(
    list(message = message, call = NULL),
    class = c("figure_input_missing", "error", "condition")
  ))
}


# 1. Plotting helpers ===========================================================

display_factor <- function(r, max_cells = maximum_map_cells) {
  if (ncell(r) <= max_cells) return(1L)
  max(1L, as.integer(ceiling(sqrt(ncell(r) / max_cells))))
}

raster_points <- function(r, categorical = FALSE, max_cells = maximum_map_cells,
                          aggregate_fun = NULL) {
  r <- r[[1]]
  factor <- display_factor(r, max_cells)
  if (factor > 1L) {
    fun <- aggregate_fun
    if (is.null(fun)) fun <- if (categorical) "modal" else "mean"
    r <- aggregate(r, fact = factor, fun = fun, na.rm = TRUE)
  }
  values <- as.data.table(as.data.frame(r, xy = TRUE, na.rm = TRUE))
  if (ncol(values) != 3L) stop("Expected one raster layer.")
  setnames(values, names(values)[3], "value")
  values
}

zone_map_plot <- function(file, title, legend = FALSE, max_cells = maximum_map_cells) {
  raster <- if (inherits(file, "SpatRaster")) file else rast(required_file(file))
  x <- raster_points(raster, TRUE, max_cells)
  x[, zone := factor(as.integer(value), levels = zone_palette$zoneID)]
  ggplot(x, aes(x, y, fill = zone)) +
    geom_raster() +
    scale_fill_manual(
      values = zone_colors, labels = zone_labels, drop = TRUE,
      na.value = "transparent", name = "Ecotype"
    ) +
    coord_equal(expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      legend.position = if (legend) "right" else "none"
    )
}

continuous_map_plot <- function(file, title, limits = NULL,
                                max_cells = maximum_map_cells) {
  raster <- if (inherits(file, "SpatRaster")) file else rast(required_file(file))
  x <- raster_points(raster, FALSE, max_cells)
  ggplot(x, aes(x, y, fill = value)) +
    geom_raster() +
    scale_fill_viridis_c(
      option = "C", limits = limits, oob = scales::squish,
      na.value = "transparent", name = "Dual\nsuitability"
    ) +
    coord_equal(expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
}

zone_suitability_plot <- function(file, title, max_cells = 100000L) {
  x <- if (inherits(file, "SpatRaster")) file else rast(required_file(file))
  if (nlyr(x) < 2L) stop("Expected zone and suitability layers: ", file)
  zone <- raster_points(x[[1]], TRUE, max_cells)
  suitability <- raster_points(x[[2]], FALSE, max_cells)
  setnames(zone, "value", "zoneID")
  setnames(suitability, "value", "suitability")
  plot_data <- merge(zone, suitability, by = c("x", "y"), all = FALSE)
  plot_data[, zone := factor(as.integer(zoneID), levels = zone_palette$zoneID)]
  ggplot(plot_data, aes(x, y, fill = zone, alpha = suitability)) +
    geom_raster() +
    scale_fill_manual(values = zone_colors, labels = zone_labels, drop = TRUE) +
    scale_alpha_continuous(range = c(0.12, 1), limits = c(0, 1)) +
    coord_equal(expand = FALSE) +
    labs(title = title, x = NULL, y = NULL, fill = "Source ecotype", alpha = "Dual suitability") +
    theme_void(base_size = 8) +
    theme(plot.title = element_text(hjust = 0.5, size = 9))
}

# Climate x soil bivariate suitability ----------------------------------------
# Dual suitability is climate suitability retained only where soil suitability
# exceeds the gate. A bivariate map therefore reads the two component rasters;
# recolouring the one-dimensional dual raster cannot recover this information.

climate_suitability_file <- function(scenario, zoneID) {
  first_existing(
    c(
      file.path(
        base_dir, "clim suitability", "mf_var", scenario,
        paste0("clim_suit_zone", zoneID, ".tif")
      ),
      file.path(
        base_dir, "suitability", "climate", scenario,
        paste0("clim_suit_zone", zoneID, ".tif")
      )
    ),
    paste("climate suitability", scenario, "zone", zoneID)
  )
}

soil_suitability_file <- function(zoneID) {
  first_existing(
    c(
      file.path(
        base_dir, "soil suitability", "plain_mf", "normal",
        paste0("soil_suit_zone", zoneID, ".tif")
      ),
      file.path(
        base_dir, "suitability", "soil", "normal",
        paste0("soil_suit_zone", zoneID, ".tif")
      )
    ),
    paste("soil suitability zone", zoneID)
  )
}

bivariate_anchors <- list(
  low_low = "#F3F1EA",
  climate_high = "#2477B3",
  soil_high = "#B56A36",
  joint_high = "#155F4A"
)

bivariate_colour <- function(climate, soil) {
  climate <- pmin(1, pmax(0, as.numeric(climate)))
  soil <- pmin(1, pmax(0, as.numeric(soil)))
  answer <- rep(NA_character_, length(climate))
  valid <- is.finite(climate) & is.finite(soil)
  if (!any(valid)) return(answer)
  
  corner <- lapply(bivariate_anchors, function(colour) {
    as.numeric(grDevices::col2rgb(colour)) / 255
  })
  c_value <- climate[valid]
  s_value <- soil[valid]
  weights <- list(
    (1 - c_value) * (1 - s_value),
    c_value * (1 - s_value),
    (1 - c_value) * s_value,
    c_value * s_value
  )
  channels <- do.call(cbind, lapply(seq_len(3L), function(channel) {
    weights[[1L]] * corner$low_low[[channel]] +
      weights[[2L]] * corner$climate_high[[channel]] +
      weights[[3L]] * corner$soil_high[[channel]] +
      weights[[4L]] * corner$joint_high[[channel]]
  }))
  answer[valid] <- grDevices::rgb(
    channels[, 1L], channels[, 2L], channels[, 3L]
  )
  answer
}

component_pair <- function(scenario, zoneID, fact = 1L, cache = NULL) {
  zoneID <- as.integer(zoneID)[[1]]
  fact <- as.integer(fact)[[1]]
  cache_key <- paste(scenario, zoneID, fact, sep = "|")
  if (!is.null(cache) && exists(cache_key, envir = cache, inherits = FALSE)) {
    return(get(cache_key, envir = cache, inherits = FALSE))
  }
  
  climate <- rast(climate_suitability_file(scenario, zoneID))[[1]]
  soil <- rast(soil_suitability_file(zoneID))[[1]]
  if (!compareGeom(climate, soil, stopOnError = FALSE)) {
    message(
      "[DISPLAY RESAMPLE] soil -> climate geometry | ",
      scenario, " | zone ", zoneID
    )
    soil <- resample(soil, climate, method = "bilinear")
  }
  pair <- c(climate, soil)
  names(pair) <- c("climate_suitability", "soil_suitability")
  if (fact > 1L) {
    pair <- aggregate(pair, fact = fact, fun = "mean", na.rm = TRUE)
  }
  if (!is.null(cache)) assign(cache_key, pair, envir = cache)
  pair
}

bivariate_pair_data <- function(scenario, zoneID, max_cells = 70000L,
                                cache = NULL) {
  climate <- rast(climate_suitability_file(scenario, zoneID))[[1]]
  fact <- display_factor(climate, max_cells)
  pair <- component_pair(scenario, zoneID, fact, cache)
  values <- as.data.table(as.data.frame(pair, xy = TRUE, na.rm = TRUE))
  values[, colour := bivariate_colour(
    climate_suitability, soil_suitability
  )]
  values[!is.na(colour)]
}

bivariate_selector_data <- function(selector, scenario, mask_file = NULL,
                                    max_cells = 70000L, cache = NULL) {
  selector <- if (inherits(selector, "SpatRaster")) {
    selector[[1]]
  } else {
    rast(required_file(selector))[[1]]
  }
  if (!is.null(mask_file)) {
    include <- if (inherits(mask_file, "SpatRaster")) {
      mask_file[[1]]
    } else {
      rast(required_file(mask_file))[[1]]
    }
    if (!compareGeom(selector, include, stopOnError = FALSE)) {
      stop("Selector and bivariate display mask have incompatible geometry.")
    }
    # Mask before modal display aggregation, so the selected source population
    # is derived only from cells that are inside the thresholded species niche.
    selector <- ifel(!is.na(include) & include > 0, selector, NA)
  }
  fact <- display_factor(selector, max_cells)
  selector_display <- if (fact > 1L) {
    aggregate(selector, fact = fact, fun = "modal", na.rm = TRUE)
  } else {
    selector
  }
  names(selector_display) <- "zoneID"
  values <- as.data.table(as.data.frame(
    selector_display, xy = TRUE, na.rm = TRUE
  ))
  values[, zoneID := as.integer(zoneID)]
  values <- values[zoneID %in% modeled_zones]
  
  if (!nrow(values)) {
    values[, `:=`(
      climate_suitability = numeric(),
      soil_suitability = numeric(),
      colour = character()
    )]
    return(values)
  }
  
  values[, `:=`(
    climate_suitability = NA_real_,
    soil_suitability = NA_real_
  )]
  for (zone in sort(unique(values$zoneID))) {
    index <- which(values$zoneID == zone)
    pair <- component_pair(scenario, zone, fact, cache)
    extracted <- terra::extract(
      pair, as.matrix(values[index, .(x, y)]), ID = FALSE
    )
    values[index, `:=`(
      climate_suitability = as.numeric(extracted[[1]]),
      soil_suitability = as.numeric(extracted[[2]])
    )]
  }
  values[, colour := bivariate_colour(
    climate_suitability, soil_suitability
  )]
  values[!is.na(colour)]
}

bivariate_map_plot <- function(data, title) {
  if (!nrow(data)) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No mapped niche") +
        labs(title = title) +
        theme_void(base_size = 9) +
        theme(plot.title = element_text(hjust = 0.5, size = 9))
    )
  }
  ggplot(data, aes(x, y, fill = colour)) +
    geom_raster() +
    scale_fill_identity() +
    coord_equal(expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(hjust = 0.5, size = 9))
}

bivariate_key_plot <- function() {
  key <- CJ(
    climate_suitability = seq(0, 1, length.out = 101L),
    soil_suitability = seq(0, 1, length.out = 101L)
  )
  key[, colour := bivariate_colour(
    climate_suitability, soil_suitability
  )]
  ggplot(key, aes(climate_suitability, soil_suitability, fill = colour)) +
    geom_raster() +
    scale_fill_identity() +
    geom_vline(
      xintercept = dual_threshold, linetype = 2,
      colour = "white", linewidth = 0.45
    ) +
    geom_hline(
      yintercept = 0.2, linetype = 2,
      colour = "white", linewidth = 0.45
    ) +
    scale_x_continuous(
      breaks = c(0, dual_threshold, 1), expand = c(0, 0)
    ) +
    scale_y_continuous(breaks = c(0, 0.2, 1), expand = c(0, 0)) +
    coord_equal() +
    labs(
      title = "Climate \u00d7 soil suitability",
      subtitle = "Dashed: climate 0.4; soil gate 0.2",
      x = "Climate suitability", y = "Soil suitability"
    ) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 9, face = "bold"),
      plot.subtitle = element_text(size = 7)
    )
}

row_proportion <- function(x, group, value) {
  x[, proportion := get(value) / sum(get(value)), by = group]
  x
}

make_chord <- function(data, from, to, value, filename, grid_colors) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required for chord diagrams.")
  }
  flow <- as.data.frame(data[, .(
    from = as.character(get(from)),
    to = as.character(get(to)),
    value = as.numeric(get(value))
  )][value > 0])
  if (!nrow(flow)) stop("No positive flows for chord diagram.")
  pdf(filename, width = 10, height = 10, useDingbats = FALSE)
  on.exit({
    circlize::circos.clear()
    dev.off()
  }, add = TRUE)
  circlize::circos.clear()
  circlize::chordDiagram(
    flow, grid.col = grid_colors, transparency = 0.65,
    annotationTrack = c("grid", "name"), directional = 1,
    direction.type = c("arrows", "diffHeight"), diffHeight = -0.03
  )
}

run_figure("zone_palette_key", {
  key <- copy(zone_palette[zoneID %in% c(modeled_zones, 99L)])
  key[, category_label := tools::toTitleCase(gsub("_", " ", category2))]
  setorder(key, category_label, zoneID)
  key[, position := seq_len(.N), by = category_label]
  key[, zone_factor := factor(zoneID, levels = zone_palette$zoneID)]
  p <- ggplot(key, aes(position, 1, fill = zone_factor)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = zoneID), size = 3) +
    facet_wrap(~ category_label, ncol = 2) +
    scale_fill_manual(values = zone_colors, guide = "none") +
    coord_equal() +
    labs(title = "Ecotype-zone colour key", x = NULL, y = NULL) +
    theme_void(base_size = 10) +
    theme(
      strip.text = element_text(size = 9, face = "bold"),
      plot.title = element_text(size = 12, face = "bold")
    )
  save_gg(p, file.path(figure_dir, "maps", "zone_key.png"), 10, 8)
})

# 2. Model and reference-map assessment ========================================

assessment_dir <- file.path(base_dir, "assessment")

# Read the existing selected-variable Multi-Forest CSVs directly. The legacy
# names/column aliases are an input interface, not a reason to rerun script 5.
# Conversions below only summarize small cached tables; they never read raster
# values, fit models, or write replacement assessment data.
assessment_csv_cache <- new.env(parent = emptyenv())

assessment_csv_candidates <- function(filename) {
  stem <- sub("\\.csv$", "", filename)
  legacy_names <- paste0(stem, "_var.csv")
  if (filename == "rf_test_zone_metrics.csv") {
    legacy_names <- c("model_test_zone_metrics_var.csv", legacy_names)
  }
  c(
    file.path(assessment_dir, filename),
    file.path(base_dir, "assessment_var", legacy_names)
  )
}

normalize_assessment_csv <- function(x, filename) {
  if ("method" %in% names(x)) {
    x <- x[which(x[["method"]] == "mf_var")]
    if (!nrow(x)) stop("No mf_var rows in ", filename)
  }
  aliases <- c(
    predicted_zone = "assigned_zone", n = "pixels",
    recall = "sensitivity", zoneID = "zone",
    exact_accuracy = "exact_zone_agreement",
    original_category2 = "original_category",
    assigned_category2 = "assigned_category"
  )
  for (old in names(aliases)) {
    new <- aliases[[old]]
    if (old %in% names(x) && !new %in% names(x)) setnames(x, old, new)
  }
  x
}

read_assessment_csv <- function(filename) {
  if (exists(filename, envir = assessment_csv_cache, inherits = FALSE)) {
    return(copy(get(filename, envir = assessment_csv_cache, inherits = FALSE)))
  }
  candidates <- assessment_csv_candidates(filename)
  available <- candidates[file.exists(candidates)]
  x <- NULL
  for (path in available) {
    candidate <- fread(path)
    # Some pre-correction main-branch tables used these same short paths for
    # plain_mf/optimized_mf and a different map threshold. Never relabel those
    # methods as mf_var; continue to the corrected assessment_var source.
    if ("method" %in% names(candidate) &&
        !any(candidate[["method"]] == "mf_var", na.rm = TRUE)) {
      cat("[IGNORE NON-MF_VAR CSV]", path, "\n")
      next
    }
    x <- normalize_assessment_csv(candidate, basename(path))
    cat("[ASSESSMENT CSV]", path, "\n")
    break
  }
  if (is.null(x) && filename == "normal_map_category_confusion_long.csv") {
    # Use exactly the same pixel domain as the cached zone confusion table.
    x <- read_assessment_csv("normal_map_confusion_long.csv")
    category_lookup <- setNames(
      as.character(zone_palette$category2), as.character(zone_palette$zoneID)
    )
    x[, `:=`(
      original_category = unname(category_lookup[as.character(original_zone)]),
      assigned_category = unname(category_lookup[as.character(assigned_zone)])
    )]
    if (anyNA(x$original_category) || anyNA(x$assigned_category)) {
      stop("Zone confusion contains zones without a category2 palette entry.")
    }
    x <- x[, .(pixels = sum(pixels)), by = .(original_category, assigned_category)]
    cat("[ASSESSMENT CSV] Category totals from existing zone confusion CSV.\n")
  } else if (is.null(x) && filename == "rf_test_zone_metrics.csv") {
    parts <- lapply(c("climate", "soil"), function(niche_name) {
      part <- read_assessment_csv(paste0(niche_name, "_test_zone_metrics.csv"))
      part[, niche := niche_name]
      part
    })
    x <- rbindlist(parts, use.names = TRUE, fill = TRUE)
  } else if (is.null(x)) {
    missing_figure_input(paste0(
      "No existing assessment CSV for ", filename, ". Checked: ",
      paste(candidates, collapse = "; "),
      ". Producer: script/5. assessment.R; visualization does not run it."
    ))
  }
  
  required <- switch(filename,
                     "rf_test_zone_metrics.csv" = c("zone", "niche", "balanced_accuracy", "f1", "tss", "auc"),
                     "normal_map_confusion_long.csv" = c("original_zone", "assigned_zone", "pixels"),
                     "normal_map_category_confusion_long.csv" = c("original_category", "assigned_category", "pixels"),
                     "normal_map_zone_metrics.csv" = c("zone", "f1", "sensitivity", "specificity", "tss"),
                     "normal_map_overall_metrics.csv" = c("coverage", "exact_zone_agreement"),
                     character()
  )
  if (!all(required %in% names(x))) {
    stop(filename, " is missing columns: ", paste(setdiff(required, names(x)), collapse = ", "))
  }
  if (!nrow(x)) stop("Empty assessment table: ", filename)
  if ("pixels" %in% names(x) &&
      (any(!is.finite(x$pixels) | x$pixels < 0) || sum(x$pixels) <= 0)) {
    stop("Invalid pixel counts in ", filename)
  }
  if (filename == "normal_map_overall_metrics.csv") {
    if (nrow(x) != 1L) stop("Expected one mf_var overall assessment row.")
    if (!"category_agreement" %in% names(x)) {
      counts <- read_assessment_csv("normal_map_confusion_long.csv")
      if ("compared_pixels" %in% names(x) &&
          !isTRUE(all.equal(sum(counts$pixels), as.numeric(x$compared_pixels[[1L]])))) {
        stop("Cached overall and zone-confusion tables describe different pixel domains.")
      }
      categories <- read_assessment_csv("normal_map_category_confusion_long.csv")
      if (!isTRUE(all.equal(sum(counts$pixels), sum(categories$pixels)))) {
        stop("Cached category and zone-confusion tables have different totals.")
      }
      x[, category_agreement := categories[
        original_category == assigned_category, sum(pixels)
      ] / sum(categories$pixels)]
    }
  }
  assign(filename, copy(x), envir = assessment_csv_cache)
  copy(x)
}

run_figure("assessment_model_performance", {
  metrics <- read_assessment_csv("rf_test_zone_metrics.csv")
  selected <- c("balanced_accuracy", "f1", "tss", "auc")
  if (!all(c("zone", "niche", selected) %in% names(metrics))) {
    stop("rf_test_zone_metrics.csv has an unexpected schema.")
  }
  long <- melt(
    metrics, id.vars = c("zone", "niche"), measure.vars = selected,
    variable.name = "metric", value.name = "value"
  )
  long[, metric := factor(
    metric, levels = selected,
    labels = c("Balanced accuracy", "F1", "TSS", "AUC")
  )]
  p <- ggplot(long, aes(niche, value, fill = niche)) +
    geom_violin(width = 0.85, alpha = 0.35, colour = NA) +
    geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.75) +
    geom_jitter(width = 0.09, size = 0.8, alpha = 0.55) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c(climate = "#2E6F9E", soil = "#9A7B32")) +
    labs(x = NULL, y = "Held-out value", fill = "Niche") +
    theme_manuscript + theme(legend.position = "top")
  save_gg(p, file.path(figure_dir, "assessment", "model_performance.png"), 11, 4.4)
  
  p_zone <- ggplot(long, aes(factor(zone), value, colour = niche, group = niche)) +
    geom_point(size = 0.8, alpha = 0.75) +
    facet_grid(metric ~ niche, scales = "free_y") +
    scale_colour_manual(values = c(climate = "#2E6F9E", soil = "#9A7B32")) +
    labs(x = "Ecotype zone", y = "Held-out value", colour = "Niche") +
    theme_manuscript +
    theme(axis.text.x = element_text(angle = 90, size = 5), legend.position = "none")
  save_gg(p_zone, file.path(figure_dir, "assessment", "model_performance_by_zone.png"), 12, 8)
})

run_figure("assessment_reference_confusion", {
  confusion <- read_assessment_csv("normal_map_confusion_long.csv")
  confusion <- row_proportion(confusion, "original_zone", "pixels")
  p <- ggplot(confusion, aes(factor(assigned_zone), factor(original_zone), fill = proportion)) +
    geom_tile() +
    scale_fill_viridis_c(option = "C", trans = "sqrt", labels = scales::percent) +
    coord_equal() +
    labs(x = "Assigned zone", y = "Observed zone", fill = "Row share") +
    theme_manuscript + theme(axis.text = element_text(size = 5))
  save_gg(p, file.path(figure_dir, "assessment", "reference_confusion.png"), 9, 8)
  
  category <- read_assessment_csv("normal_map_category_confusion_long.csv")
  category <- row_proportion(category, "original_category", "pixels")
  p_category <- ggplot(
    category, aes(assigned_category, original_category, fill = proportion)
  ) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_viridis_c(option = "C", labels = scales::percent) +
    coord_equal() +
    labs(x = "Assigned category", y = "Observed category", fill = "Row share") +
    theme_manuscript +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_gg(
    p_category,
    file.path(figure_dir, "assessment", "reference_category_confusion.png"),
    9, 7
  )
})

run_figure("assessment_reference_metrics", {
  binary <- read_assessment_csv("normal_map_zone_metrics.csv")
  selected <- c("f1", "sensitivity", "specificity", "tss")
  long <- melt(
    binary, id.vars = "zone", measure.vars = selected,
    variable.name = "metric", value.name = "value"
  )
  p <- ggplot(long, aes(factor(zone), value, colour = factor(zone))) +
    geom_point(size = 1.4) +
    facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    scale_colour_manual(values = zone_colors, guide = "none") +
    labs(x = "Ecotype zone", y = "Reference-map metric") +
    theme_manuscript + theme(axis.text.x = element_text(angle = 90, size = 6))
  save_gg(p, file.path(figure_dir, "assessment", "reference_zone_metrics.png"), 11, 9)
  
  overall <- read_assessment_csv("normal_map_overall_metrics.csv")
  values <- melt(
    overall,
    measure.vars = c("coverage", "exact_zone_agreement", "category_agreement"),
    variable.name = "metric", value.name = "value"
  )
  values[, metric := factor(
    metric,
    levels = c("coverage", "exact_zone_agreement", "category_agreement"),
    labels = c("Coverage", "Exact zone", "Broad category")
  )]
  p_overall <- ggplot(values, aes(metric, value, fill = metric)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = scales::percent(value, accuracy = 0.1)), vjust = -0.3) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.06)) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(x = NULL, y = "Agreement") + theme_manuscript
  save_gg(p_overall, file.path(figure_dir, "assessment", "reference_overall.png"), 6.5, 4.5)
  
  multiclass_file <- multiclass_table_file("reference_map_by_zone.csv")
  if (file.exists(multiclass_file)) {
    multiclass <- fread(multiclass_file)
    if ("zoneID" %in% names(multiclass)) setnames(multiclass, "zoneID", "zone")
    binary_compare <- copy(binary)[, workflow := "Binary Multi-Forest"]
    multiclass[, workflow := "Multiclass RF"]
    comparison <- rbindlist(
      list(
        binary_compare[, c("zone", selected, "workflow"), with = FALSE],
        multiclass[, c("zone", selected, "workflow"), with = FALSE]
      ),
      use.names = TRUE
    )
    comparison <- melt(
      comparison, id.vars = c("zone", "workflow"), measure.vars = selected,
      variable.name = "metric", value.name = "value"
    )
    p_compare <- ggplot(
      comparison,
      aes(factor(zone), value, colour = workflow, group = workflow)
    ) +
      geom_point(size = 0.9, alpha = 0.8) +
      facet_wrap(~ metric, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = c(
        "Binary Multi-Forest" = "#2E6F9E", "Multiclass RF" = "#D95F02"
      )) +
      labs(x = "Ecotype zone", y = "Reference-map metric", colour = "Workflow") +
      theme_manuscript +
      theme(axis.text.x = element_text(angle = 90, size = 5), legend.position = "top")
    save_gg(
      p_compare,
      file.path(figure_dir, "assessment", "binary_multiclass_reference_metrics.png"),
      11, 9
    )
  }
})

run_figure("assessment_chords", {
  zone_flow <- read_assessment_csv("normal_map_confusion_long.csv")
  make_chord(
    zone_flow, "original_zone", "assigned_zone", "pixels",
    file.path(figure_dir, "assessment", "reference_zone_chord.pdf"),
    zone_colors
  )
  category_flow <- read_assessment_csv("normal_map_category_confusion_long.csv")
  categories <- sort(unique(c(
    category_flow$original_category, category_flow$assigned_category
  )))
  category_colors <- setNames(grDevices::hcl.colors(length(categories), "Dark 3"), categories)
  make_chord(
    category_flow, "original_category", "assigned_category", "pixels",
    file.path(figure_dir, "assessment", "reference_category_chord.pdf"),
    category_colors
  )
})

run_figure("soil_gate_sensitivity", {
  gate_file <- file.path(assessment_dir, "soil_gate_sensitivity.csv")
  if (!file.exists(gate_file)) {
    missing_figure_input(paste(
      "No soil_gate_sensitivity.csv. This separate sensitivity analysis is",
      "produced by script 3.2 and may require a full raster scan.",
      "It is not a script 7 projection output. No analysis was started."
    ))
  }
  gate <- fread(gate_file)
  if (all(c("below_0.4_area_km2", "valid_area_km2") %in% names(gate))) {
    gate[, below_threshold_share := fifelse(
      is.finite(valid_area_km2) & valid_area_km2 > 0,
      below_0.4_area_km2 / valid_area_km2, NA_real_
    )]
    gate_y_label <- "Area share below dual-suitability 0.4"
  } else if ("below_0.4_share" %in% names(gate)) {
    gate[, below_threshold_share := below_0.4_share]
    gate_y_label <- "Pixel share below dual-suitability 0.4"
  } else {
    stop("soil_gate_sensitivity.csv has neither area nor pixel shares.")
  }
  gate[, scenario_label := factor(scenario_label(scenario), levels = scenario_label(scenarios))]
  p <- ggplot(gate, aes(soil_gate, below_threshold_share, colour = scenario_label)) +
    geom_vline(xintercept = 0.2, linetype = 2, colour = "grey35") +
    geom_line(linewidth = 0.75) + geom_point(size = 1.5) +
    scale_x_continuous(breaks = unique(gate$soil_gate)) +
    scale_y_continuous(labels = scales::percent) +
    scale_colour_viridis_d(option = "D", end = 0.9) +
    labs(
      x = "Soil gate", y = gate_y_label,
      colour = "Scenario",
      caption = "Gate 0 is climate-only; the dashed line is the primary gate (0.2)."
    ) + theme_manuscript + theme(legend.position = "right")
  save_gg(p, file.path(figure_dir, "assessment", "soil_gate_sensitivity.png"), 9.5, 6)
})


# 3. Feature importance =========================================================

importance_table_dir <- file.path(assessment_dir, "feature_importance", "tables")

# Importance is displayed at four resolutions, always after normalization
# within niche and with binary Multi-Forest and multiclass RF equally weighted:
#   1. zone x variable;              2. vegetation category x variable;
#   3. zone x variable group;        4. vegetation category x variable group.
importance_category_order <- c(
  "forest", "scrub", "grassland_meadow_steppe", "desert",
  "cropland", "wetland", "alpine_vegetation"
)
importance_group_order <- c(
  "Temperature", "Precipitation", "Climatic moisture balance",
  "Growing season and degree-days", "Radiation and evaporative demand",
  "Humidity", "Snowfall",
  "Nutrient retention and base status", "Bulk density", "Coarse fragments",
  "Texture", "Organic carbon", "Soil reaction", "Carbonates", "Gypsum",
  "Sodicity", "Salinity", "Other"
)
importance_group_colours <- c(
  "Temperature" = "#D55E00",
  "Precipitation" = "#0072B2",
  "Climatic moisture balance" = "#009E73",
  "Growing season and degree-days" = "#6A994E",
  "Radiation and evaporative demand" = "#E69F00",
  "Humidity" = "#56B4E9",
  "Snowfall" = "#9FBAD0",
  "Nutrient retention and base status" = "#2E8B57",
  "Bulk density" = "#6B4C3B",
  "Coarse fragments" = "#7A7F87",
  "Texture" = "#C17C2C",
  "Organic carbon" = "#4E342E",
  "Soil reaction" = "#7B3294",
  "Carbonates" = "#C9A227",
  "Gypsum" = "#CC79A7",
  "Sodicity" = "#E76F51",
  "Salinity" = "#008C95",
  "Other" = "#BDBDBD"
)

importance_category_factor <- function(x, reverse = FALSE) {
  values <- as.character(x)
  present <- unique(values)
  ordered <- c(
    intersect(importance_category_order, present),
    sort(setdiff(present, importance_category_order))
  )
  labels <- tools::toTitleCase(gsub("_", " ", ordered))
  if (reverse) labels <- rev(labels)
  factor(tools::toTitleCase(gsub("_", " ", values)), levels = labels)
}

importance_group_style <- function(x) {
  present <- unique(as.character(x))
  levels <- c(
    intersect(importance_group_order, present),
    sort(setdiff(present, importance_group_order))
  )
  missing_colours <- setdiff(levels, names(importance_group_colours))
  if (length(missing_colours)) {
    stop(
      "No semantic colour is defined for variable category: ",
      paste(missing_colours, collapse = ", ")
    )
  }
  list(
    levels = levels,
    colours = importance_group_colours[levels]
  )
}

run_figure("importance_multiclass_global", {
  global_importance <- read_required(
    file.path(importance_table_dir, "FI_03_multiclass_global_importance.csv")
  )
  top <- global_importance[
    order(niche, -importance_share_within_niche),
    head(.SD, 12L),
    by = niche
  ]
  top[, variable_label := factor(
    variable, levels = rev(unique(top[order(importance_share_within_niche), variable]))
  )]
  p <- ggplot(top, aes(importance_share_within_niche, variable_label, fill = niche)) +
    geom_col(width = 0.72) +
    facet_wrap(~ niche, scales = "free_y", nrow = 1) +
    scale_x_continuous(labels = scales::percent) +
    scale_fill_manual(values = c(climate = "#2E6F9E", soil = "#9A7B32"), guide = "none") +
    labs(
      x = "Within-niche global multiclass importance", y = NULL,
      title = "Global multiclass RF permutation importance"
    ) + theme_manuscript
  save_gg(p, file.path(figure_dir, "importance", "multiclass_global.png"), 10, 6.5)
})

run_figure("importance_agreement", {
  agreement <- read_required(file.path(importance_table_dir, "FI_06_zone_agreement.csv"))
  agreement_long <- melt(
    agreement,
    id.vars = c("zoneID", "category2", "niche"),
    measure.vars = c("spearman", "cosine", "top5_jaccard"),
    variable.name = "measure", value.name = "agreement"
  )
  agreement_long[, measure := factor(
    measure, levels = c("spearman", "cosine", "top5_jaccard"),
    labels = c("Spearman", "Cosine similarity", "Top-5 Jaccard")
  )]
  p <- ggplot(
    agreement_long,
    aes(factor(zoneID), agreement, colour = category2)
  ) +
    geom_point(size = 1.25) + facet_grid(measure ~ niche, scales = "free_y") +
    labs(
      x = "Ecotype zone", y = "Binary Multi-Forest vs multiclass agreement",
      colour = "Vegetation category"
    ) + theme_manuscript +
    theme(axis.text.x = element_text(angle = 90, size = 5), legend.position = "bottom")
  save_gg(p, file.path(figure_dir, "importance", "binary_multiclass_agreement.png"), 12, 8)
})

# Level 1: every reader-facing vegetation zone x every individual variable.
# The no-vegetation category is excluded; variables absent from a zone-specific
# binary model remain explicit zeros in the aligned table.
run_figure("importance_zone_variables", {
  zone_variables <- read_required(
    file.path(importance_table_dir, "FI_reader_03_aligned_zone_profiles.csv")
  )
  maximum_importance <- max(zone_variables$consensus_share, na.rm = TRUE)
  if (!is.finite(maximum_importance) || maximum_importance <= 0) {
    maximum_importance <- 1
  }
  for (niche_value in c("climate", "soil")) {
    plot_data <- copy(zone_variables[niche == niche_value])
    plot_data[, category_label := importance_category_factor(category2)]
    plot_data[, zone_label := factor(
      sprintf("Zone %02d", zoneID),
      levels = rev(sprintf("Zone %02d", sort(unique(zoneID))))
    )]
    
    variable_order <- plot_data[, .(
      mean_importance = mean(consensus_share, na.rm = TRUE)
    ), by = .(variable_group, variable)][
      order(variable_group, -mean_importance, variable), variable
    ]
    group_style <- importance_group_style(plot_data$variable_group)
    plot_data[, variable_group := factor(
      variable_group, levels = group_style$levels
    )]
    plot_data[, variable := factor(variable, levels = variable_order)]
    
    p <- ggplot(plot_data, aes(variable, zone_label, fill = consensus_share)) +
      geom_tile() +
      facet_grid(
        category_label ~ variable_group,
        scales = "free", space = "free", switch = "y",
        labeller = labeller(
          category_label = label_wrap_gen(18),
          variable_group = label_wrap_gen(18)
        )
      ) +
      scale_fill_viridis_c(
        option = "C", trans = "sqrt",
        limits = c(0, maximum_importance), oob = scales::squish,
        labels = scales::percent,
        name = "Normalized MDA\nconsensus share"
      ) +
      labs(
        title = paste(
          "Zone-level",
          if (niche_value == "climate") "climatic" else "topsoil",
          "variable importance"
        ),
        subtitle = paste(
          "Rows retain every vegetation zone (no-vegetation excluded);",
          "columns retain every variable.",
          "Facets show vegetation and variable groups without aggregation."
        ),
        caption = paste(
          "A square-root colour transform reveals small non-zero values;",
          "legend labels remain normalized importance shares."
        ),
        x = NULL, y = NULL
      ) +
      theme_bw(base_size = 9) +
      theme(
        panel.grid = element_blank(),
        panel.spacing = grid::unit(0.8, "mm"),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 6.5),
        axis.text.y = element_text(size = 6),
        strip.text.x = element_text(size = 7, face = "bold"),
        strip.text.y.left = element_text(angle = 0, size = 7, face = "bold"),
        plot.title = element_text(face = "bold")
      )
    
    figure_width <- max(15, 5 + 0.30 * uniqueN(plot_data$variable))
    figure_height <- max(13, 4 + 0.19 * uniqueN(plot_data$zoneID))
    save_gg(
      p,
      file.path(
        figure_dir, "importance",
        paste0("zone_variables_", niche_value, ".png")
      ),
      figure_width, figure_height
    )
  }
})

# Level 3: individual variables are summed within group separately for every
# reader-facing vegetation zone.
run_figure("importance_zone_groups", {
  zone_groups <- read_required(
    file.path(importance_table_dir, "FI_reader_05_zone_grouped_importance.csv")
  )
  for (niche_value in c("climate", "soil")) {
    plot_data <- copy(zone_groups[niche == niche_value])
    plot_data[, category_label := importance_category_factor(category2)]
    plot_data[, zone_label := factor(
      sprintf("Zone %02d", zoneID),
      levels = rev(sprintf("Zone %02d", sort(unique(zoneID))))
    )]
    group_style <- importance_group_style(plot_data$variable_group)
    plot_data[, variable_group := factor(
      variable_group, levels = group_style$levels
    )]
    
    p <- ggplot(
      plot_data,
      aes(x = consensus_share, y = zone_label, fill = variable_group)
    ) +
      geom_col(
        position = position_stack(reverse = TRUE),
        width = 0.74, colour = "white", linewidth = 0.08
      ) +
      facet_grid(
        category_label ~ ., scales = "free_y", space = "free_y", switch = "y",
        labeller = labeller(category_label = label_wrap_gen(18))
      ) +
      scale_x_continuous(
        breaks = seq(0, 1, by = 0.2), labels = scales::percent,
        expand = expansion(mult = c(0, 0))
      ) +
      scale_fill_manual(
        values = group_style$colours, drop = FALSE, name = "Variable category"
      ) +
      coord_cartesian(xlim = c(0, 1)) +
      labs(
        title = paste(
          "Zone-level",
          if (niche_value == "climate") "climatic" else "topsoil",
          "importance by variable group"
        ),
        subtitle = paste(
          "Segments partition normalized within-niche importance after variables",
          "are summed within each group and zone."
        ),
        x = "Normalized MeanDecreaseAccuracy (consensus share)", y = NULL
      ) +
      theme_bw(base_size = 9) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.spacing.y = grid::unit(1.2, "mm"),
        axis.text.y = element_text(size = 6),
        strip.text.y.left = element_text(angle = 0, size = 7, face = "bold"),
        legend.position = "right",
        plot.title = element_text(face = "bold")
      )
    
    save_gg(
      p,
      file.path(
        figure_dir, "importance",
        paste0("zone_variable_groups_", niche_value, ".png")
      ),
      12, 14
    )
  }
})

# Level 4: average the within-zone group totals across zones in each vegetation
# category. Climate and soil remain separate compositions.
run_figure("importance_category_groups", {
  grouped <- read_required(
    file.path(importance_table_dir, "FI_reader_06_category_grouped_importance.csv")
  )
  group_plots <- list()
  for (niche_value in c("climate", "soil")) {
    plot_data <- copy(grouped[niche == niche_value])
    plot_data[, category_label := importance_category_factor(category2, reverse = TRUE)]
    group_style <- importance_group_style(plot_data$variable_group)
    plot_data[, variable_group := factor(
      variable_group, levels = group_style$levels
    )]
    
    group_plot <- ggplot(
      plot_data,
      aes(x = consensus_mean, y = category_label, fill = variable_group)
    ) +
      geom_col(
        position = position_stack(reverse = TRUE),
        width = 0.72, colour = "white", linewidth = 0.15
      ) +
      scale_x_continuous(
        breaks = seq(0, 1, by = 0.2),
        labels = scales::percent,
        expand = expansion(mult = c(0, 0))
      ) +
      scale_fill_manual(values = group_style$colours, drop = FALSE) +
      labs(
        x = "Normalized MeanDecreaseAccuracy (consensus share)", y = NULL,
        fill = "Variable category",
        title = paste(
          if (niche_value == "climate") "Climatic" else "Topsoil",
          "importance by vegetation category"
        ),
        subtitle = paste(
          "Each horizontal bar is partitioned by the summed importance",
          "of its variable categories."
        )
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "right",
        plot.title = element_text(face = "bold")
      ) +
      coord_cartesian(xlim = c(0, 1))
    
    group_plots[[niche_value]] <- group_plot
    save_gg(
      group_plot,
      file.path(
        figure_dir, "importance",
        paste0("category_variable_groups_", niche_value, ".png")
      ),
      11, 6.5
    )
  }
  combined_groups <- wrap_plots(group_plots, ncol = 1)
  save_gg(
    combined_groups,
    file.path(figure_dir, "importance", "category_predictor_groups.png"),
    11, 12
  )
  
  p_compare <- ggplot(grouped, aes(binary_mean, multiclass_mean, colour = category2)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45") +
    geom_point(aes(shape = niche), size = 2, alpha = 0.8) +
    facet_wrap(~ niche, scales = "free") +
    labs(
      x = "Binary Multi-Forest grouped importance",
      y = "Multiclass grouped importance", colour = "Category", shape = "Niche"
    ) + theme_manuscript + theme(legend.position = "bottom")
  save_gg(p_compare, file.path(figure_dir, "importance", "grouped_consensus.png"), 10, 6.5)
})

# Level 2: average individual-variable importance across zones in each
# vegetation category. The figure retains the six leading variables per panel;
# the complete category x variable values remain in the source CSV.
run_figure("importance_category_variables", {
  variables <- read_required(
    file.path(importance_table_dir, "FI_reader_07_category_variable_importance.csv")
  )
  top <- copy(variables[consensus_rank <= 6L])
  top[, `:=`(
    category_label = importance_category_factor(category2),
    lower_share = pmax(0, consensus_mean - consensus_se),
    upper_share = consensus_mean + consensus_se
  )]
  
  variable_plots <- list()
  for (niche_value in c("climate", "soil")) {
    plot_data <- copy(top[niche == niche_value])
    setorder(plot_data, category_label, consensus_mean, variable)
    plot_data[, variable_panel := factor(
      paste(category2, variable, sep = "|||"),
      levels = unique(paste(category2, variable, sep = "|||"))
    )]
    group_style <- importance_group_style(plot_data$variable_group)
    plot_data[, variable_group := factor(
      variable_group, levels = group_style$levels
    )]
    
    variable_plot <- ggplot(
      plot_data,
      aes(x = consensus_mean, y = variable_panel, fill = variable_group)
    ) +
      geom_col(width = 0.66, colour = "white", linewidth = 0.15) +
      geom_errorbar(
        data = plot_data[is.finite(lower_share) & is.finite(upper_share)],
        aes(xmin = lower_share, xmax = upper_share),
        orientation = "y", width = 0.16, linewidth = 0.42, colour = "grey25"
      ) +
      facet_wrap(vars(category_label), scales = "free_y", ncol = 2) +
      scale_y_discrete(labels = function(x) sub("^.*\\|\\|\\|", "", x)) +
      scale_x_continuous(
        labels = scales::percent,
        expand = expansion(mult = c(0, 0.08))
      ) +
      scale_fill_manual(
        values = group_style$colours, drop = FALSE, name = "Variable category"
      ) +
      labs(
        x = "Mean normalized MeanDecreaseAccuracy (consensus)", y = NULL,
        title = paste0(
          if (niche_value == "climate") "Climatic" else "Topsoil",
          " variables contributing to ecotype-niche models"
        ),
        subtitle = paste(
          "Top six variables per vegetation category; bars are category means",
          "and error bars show +/-1 SE."
        )
      ) +
      theme_bw(base_size = 9.5) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "right",
        panel.spacing = grid::unit(0.75, "lines"),
        plot.title = element_text(face = "bold")
      )
    
    variable_plots[[niche_value]] <- variable_plot
    save_gg(
      variable_plot,
      file.path(
        figure_dir, "importance",
        paste0("category_variables_", niche_value, ".png")
      ),
      10.5, 10
    )
  }
  combined_variables <- wrap_plots(variable_plots, ncol = 2)
  save_gg(
    combined_variables,
    file.path(figure_dir, "importance", "category_variables.png"),
    21, 10
  )
})

run_figure("importance_climate_period", {
  period <- read_required(
    file.path(importance_table_dir, "FI_reader_09_category_climate_period.csv")
  )
  period[, climate_period := factor(
    climate_period, levels = c("Annual", "Winter", "Spring", "Summer", "Autumn")
  )]
  p <- ggplot(period, aes(category2, consensus_mean, fill = climate_period)) +
    geom_col(position = "stack") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_brewer(palette = "Spectral", direction = -1, na.translate = FALSE) +
    labs(
      x = "Vegetation category", y = "Climate consensus importance",
      fill = "Climate period"
    ) + theme_manuscript +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "top")
  save_gg(p, file.path(figure_dir, "importance", "climate_period.png"), 10.5, 6)
})

run_figure("importance_category_similarity", {
  variables <- read_required(
    file.path(importance_table_dir, "FI_reader_07_category_variable_importance.csv")
  )
  similarity <- rbindlist(lapply(c("climate", "soil"), function(niche_value) {
    wide <- dcast(
      variables[niche == niche_value],
      category2 ~ variable,
      value.var = "consensus_mean",
      fill = 0
    )
    categories <- wide$category2
    matrix_values <- as.matrix(
      wide[, setdiff(names(wide), "category2"), with = FALSE]
    )
    rownames(matrix_values) <- categories
    denominator <- sqrt(rowSums(matrix_values ^ 2))
    cosine <- tcrossprod(matrix_values) / outer(denominator, denominator)
    dimnames(cosine) <- list(categories, categories)
    answer <- as.data.table(as.table(cosine))
    setnames(answer, c("category_x", "category_y", "cosine"))
    answer[, niche := niche_value]
    answer
  }))
  p <- ggplot(similarity, aes(category_x, category_y, fill = cosine)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    facet_wrap(~ niche, nrow = 1) +
    scale_fill_viridis_c(option = "C", limits = c(0, 1)) +
    coord_equal() +
    labs(
      x = NULL, y = NULL,
      fill = "Cosine\nsimilarity",
      title = "Similarity of category-level importance profiles"
    ) + theme_manuscript +
    theme(
      axis.text.x = element_text(angle = 55, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7)
    )
  save_gg(p, file.path(figure_dir, "importance", "category_similarity.png"), 13, 6.5)
})

run_figure("importance_area_concentration", {
  aligned <- read_required(
    file.path(importance_table_dir, "FI_reader_03_aligned_zone_profiles.csv")
  )
  concentration <- aligned[, .(
    hhi = sum(consensus_share ^ 2),
    effective_predictors = if (sum(consensus_share ^ 2) > 0) {
      1 / sum(consensus_share ^ 2)
    } else NA_real_,
    top5_share = sum(sort(consensus_share, decreasing = TRUE)[seq_len(min(5L, .N))])
  ), by = .(zoneID, category2, niche)]
  zone_size <- read_zone_size()
  if (!all(c("zoneID", "n_cells") %in% names(zone_size))) {
    stop("data/zone_lookup.csv needs zoneID and n_cells.")
  }
  concentration <- merge(
    concentration,
    zone_size[, .(zoneID, reference_cells = n_cells)],
    by = "zoneID",
    all.x = TRUE
  )
  p <- ggplot(
    concentration[is.finite(reference_cells) & reference_cells > 0],
    aes(reference_cells, effective_predictors, colour = category2)
  ) +
    geom_point(size = 2, alpha = 0.8) +
    geom_text(aes(label = zoneID), size = 2.4, nudge_y = 0.15, check_overlap = TRUE) +
    facet_wrap(~ niche, scales = "free_y") +
    scale_x_log10(labels = scales::label_number()) +
    labs(
      x = "Observed reference area (raster cells, log scale)",
      y = "Effective number of consensus predictors (1 / HHI)",
      colour = "Vegetation category"
    ) + theme_manuscript + theme(legend.position = "bottom")
  save_gg(p, file.path(figure_dir, "importance", "area_concentration.png"), 10.5, 7)
})

run_figure("importance_representative_zones", {
  aligned <- read_required(
    file.path(importance_table_dir, "FI_reader_03_aligned_zone_profiles.csv")
  )
  profile <- dcast(
    aligned,
    zoneID + category2 ~ variable_key,
    value.var = "consensus_share",
    fill = 0
  )
  profile_columns <- setdiff(names(profile), c("zoneID", "category2"))
  category_values <- sort(setdiff(unique(profile$category2), "no_vegetation"))
  zone_size <- read_zone_size()
  
  for (category_value in category_values) {
    category_profile <- profile[category2 == category_value]
    profile_matrix <- as.matrix(category_profile[, ..profile_columns])
    centroid <- colMeans(profile_matrix)
    distance <- rowSums((profile_matrix - matrix(
      centroid,
      nrow = nrow(profile_matrix), ncol = length(centroid), byrow = TRUE
    )) ^ 2)
    representative <- category_profile$zoneID[which.min(distance)]
    size_table <- zone_size[
      zoneID %in% category_profile$zoneID & is.finite(n_cells),
      .(zoneID, count = n_cells)
    ]
    selected_zones <- unique(c(
      representative,
      size_table[which.min(count), zoneID],
      size_table[which.max(count), zoneID]
    ))
    detail <- aligned[
      zoneID %in% selected_zones & consensus_share > 0
    ][order(zoneID, niche, -consensus_share), head(.SD, 6L), by = .(zoneID, niche)]
    detail[, variable := factor(variable, levels = rev(unique(variable)))]
    p <- ggplot(detail, aes(consensus_share, variable, colour = variable_group)) +
      geom_segment(aes(x = 0, xend = consensus_share, yend = variable), colour = "grey75") +
      geom_point(size = 2) +
      facet_grid(niche ~ zoneID, scales = "free_y", space = "free_y") +
      scale_x_continuous(labels = scales::percent) +
      labs(
        title = tools::toTitleCase(gsub("_", " ", category_value)),
        subtitle = "Category medoid plus smallest and largest observed-area zones",
        x = "Consensus importance", y = NULL, colour = "Predictor group"
      ) + theme_manuscript + theme(legend.position = "bottom")
    save_gg(
      p,
      file.path(
        figure_dir, "importance",
        paste0("category_zones_", safe_name(category_value), ".png")
      ),
      12, 7
    )
  }
})

# 4. Reference and future maps ==================================================

run_figure("reference_maps", {
  files <- c(
    reference_file,
    binary_map_file("normal"),
    multiclass_map_file("normal")
  )
  titles <- c("Observed reference ecotypes", "Binary Multi-Forest overlay", "Multiclass RF")
  panels <- Map(
    function(file, title) zone_map_plot(file, title, legend = FALSE),
    files, titles
  )
  combined <- wrap_plots(panels, nrow = 1, guides = "collect") &
    theme(legend.position = "right")
  save_gg(combined, file.path(figure_dir, "maps", "reference_models.png"), 14, 5.2)
})

run_figure("future_binary_maps", {
  panels <- lapply(future_scenarios, function(scenario) {
    zone_map_plot(
      binary_map_file(scenario),
      scenario_label(scenario)
    )
  })
  combined <- wrap_plots(panels, ncol = 3, guides = "collect") &
    theme(legend.position = "right")
  save_gg(combined, file.path(figure_dir, "maps", "future_binary.png"), 14, 8.5)
})

run_figure("zone_bivariate_niches", {
  if (!length(bivariate_zones) || !length(bivariate_scenarios)) {
    skip_figure(paste(
      "Set ECOCHINA2_BIVARIATE_ZONES and ECOCHINA2_BIVARIATE_SCENARIOS",
      "to draw zone-level climate x soil maps."
    ))
  } else {
    scenario_tag <- if (setequal(bivariate_scenarios, scenarios)) {
      "all_scenarios"
    } else {
      safe_name(paste(bivariate_scenarios, collapse = "__"))
    }
    for (zone in bivariate_zones) {
      component_cache <- new.env(parent = emptyenv())
      panels <- lapply(bivariate_scenarios, function(scenario) {
        bivariate_map_plot(
          bivariate_pair_data(
            scenario, zone, max_cells = 90000L,
            cache = component_cache
          ),
          scenario_label(scenario)
        )
      })
      map_grid <- wrap_plots(
        panels, ncol = min(4L, max(1L, length(panels)))
      )
      combined <- wrap_plots(
        list(map_grid, bivariate_key_plot()),
        ncol = 2, widths = c(6, 1.25)
      ) +
        plot_annotation(
          title = paste0(
            "Zone ", zone,
            ": climate \u00d7 soil suitability underlying the gated dual niche"
          ),
          caption = paste0(
            "Blue identifies climatic opportunity with weak soil support; ",
            "ochre identifies soil support with weak climatic suitability; ",
            "green indicates joint support. Soil is static across scenarios. ",
            "Raster aggregation is for display only."
          )
        )
      save_gg(
        combined,
        file.path(
          figure_dir, "maps",
          paste0(
            "zone", zone, "_climate_soil_bivariate_",
            scenario_tag, ".png"
          )
        ),
        18,
        max(6.5, 4.1 * ceiling(length(panels) / 4) + 1.5)
      )
      rm(component_cache, panels, map_grid, combined)
      gc()
    }
  }
})

run_figure("future_multiclass_maps", {
  panels <- lapply(future_scenarios, function(scenario) {
    zone_map_plot(
      multiclass_map_file(scenario),
      scenario_label(scenario)
    )
  })
  combined <- wrap_plots(panels, ncol = 3, guides = "collect") &
    theme(legend.position = "right")
  save_gg(combined, file.path(figure_dir, "maps", "future_multiclass.png"), 14, 8.5)
})

run_figure("binary_multiclass_projection_comparison", {
  comparison <- read_required(
    multiclass_table_file("binary_mf_vs_multiclass_common_mask.csv")
  )
  comparison[, scenario_label := factor(scenario_label(scenario), levels = scenario_label(scenarios))]
  long <- melt(
    comparison,
    id.vars = "scenario_label",
    measure.vars = c("exact_agreement_all", "exact_agreement_current_zones"),
    variable.name = "scope", value.name = "agreement"
  )
  long[, scope := factor(
    scope,
    levels = c("exact_agreement_all", "exact_agreement_current_zones"),
    labels = c("All common cells", "Current-zone cells")
  )]
  p <- ggplot(long, aes(scenario_label, agreement, fill = scope)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.68) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = NULL, y = "Binary–multiclass agreement", fill = "Comparison") +
    theme_manuscript +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top")
  save_gg(p, file.path(figure_dir, "maps", "binary_multiclass_agreement.png"), 10, 5.8)
})

run_figure("future_area_and_transitions", {
  ecosystem <- keep_multiforest(read_required(niche_table_file("ecosystem")))
  ecosystem[, scenario_label := factor(scenario_label(scenario), levels = scenario_label(scenarios))]
  novel <- ecosystem[zoneID == 99]
  if (!nrow(novel)) stop("No Zone 99 rows in ecosystem_area.csv.")
  p_novel <- ggplot(novel, aes(scenario_label, area_km2 / 1e6, fill = scenario_label)) +
    geom_col(width = 0.68) +
    scale_fill_viridis_d(guide = "none") +
    labs(
      x = NULL, y = expression("Novel area (10"^6*" km"^2*")"),
      caption = "Novel means all 53 current ecotypes have dual suitability < 0.4."
    ) + theme_manuscript + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_gg(p_novel, file.path(figure_dir, "maps", "novel_area.png"), 8.5, 5.2)
  
  transition <- keep_multiforest(read_required(niche_table_file("transition")))
  if (!"change_class" %in% names(transition) && "transition_type" %in% names(transition)) {
    setnames(transition, "transition_type", "change_class")
  }
  transition_summary <- transition[, .(area_km2 = sum(area_km2)), by = .(scenario, change_class)]
  transition_summary[, share := area_km2 / sum(area_km2), by = scenario]
  transition_summary[, scenario_label := factor(
    scenario_label(scenario), levels = scenario_label(future_scenarios)
  )]
  p_transition <- ggplot(
    transition_summary,
    aes(scenario_label, share, fill = factor(change_class, c("stable", "changed", "novel")))
  ) +
    geom_col(width = 0.75) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c(stable = "#3B8C6E", changed = "#E3A43B", novel = "#252525")) +
    labs(x = NULL, y = "Reference-to-future area share", fill = "Class") +
    theme_manuscript + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_gg(p_transition, file.path(figure_dir, "maps", "transition_share.png"), 9, 5.5)
})

run_figure("multiclass_area_change", {
  area <- read_required(
    multiclass_table_file("multiclass_area_change.csv")
  )
  area <- area[scenario != "normal" & is.finite(change_percent)]
  area[, scenario_label := factor(
    scenario_label(scenario), levels = scenario_label(future_scenarios)
  )]
  p <- ggplot(area, aes(scenario_label, factor(zoneID), fill = change_percent)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0) +
    labs(x = NULL, y = "Ecotype zone", fill = "Area change\n(%)") +
    theme_manuscript + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_gg(p, file.path(figure_dir, "maps", "multiclass_area_change.png"), 10, 8)
})


# 5. Population and species niche figures ======================================

run_figure("species_population_reference", {
  populations <- read_required(first_existing(
    c(
      file.path(
        base_dir, "future tree niche var", "tables",
        "population_projection_lookup_var.csv"
      ),
      file.path(base_dir, "species_zone_population_long.csv"),
      file.path(base_dir, "data", "species_zone_population_long.csv"),
      file.path(base_dir, "data", "processed", "species_zone_population_long.csv")
    ),
    "population lookup"
  ))
  if (!"zoneID" %in% names(populations)) {
    if ("source_zone" %in% names(populations)) {
      populations[, zoneID := as.integer(source_zone)]
    } else if ("Zone" %in% names(populations)) {
      populations[, zoneID := as.integer(gsub("[^0-9]", "", Zone))]
    } else {
      stop("Population lookup needs zoneID, source_zone or Zone.")
    }
  }
  if (!"reference_abundance" %in% names(populations) &&
      "Population" %in% names(populations)) {
    populations[, reference_abundance := as.numeric(Population)]
  }
  if (!"reference_abundance" %in% names(populations)) {
    stop("Population lookup needs Population or reference_abundance.")
  }
  p <- ggplot(
    populations,
    aes(factor(zoneID), Species, fill = log10(reference_abundance + 1))
  ) +
    geom_tile(colour = "white", linewidth = 0.15) +
    scale_fill_viridis_c(option = "C") +
    labs(x = "Source ecotype zone", y = "Species", fill = "log10(cells + 1)") +
    theme_manuscript + theme(axis.text.x = element_text(angle = 90, size = 6))
  save_gg(p, file.path(figure_dir, "niche", "reference_populations.png"), 11, 5.5)
})

run_figure("species_niche_area", {
  species <- keep_multiforest(read_required(niche_table_file("species")))
  if (!"n_populations" %in% names(species) && "populations_projected" %in% names(species)) {
    setnames(species, "populations_projected", "n_populations")
  }
  species[, scenario_label := factor(scenario_label(scenario), levels = scenario_label(scenarios))]
  long <- melt(
    species,
    id.vars = c("Species", "scenario_label"),
    measure.vars = c("suitable_area_km2", "suitability_weighted_area_km2"),
    variable.name = "area_type", value.name = "area_km2"
  )
  long[, area_type := factor(
    area_type,
    levels = c("suitable_area_km2", "suitability_weighted_area_km2"),
    labels = c("Thresholded area", "Suitability-weighted area")
  )]
  p <- ggplot(long, aes(scenario_label, area_km2 / 1e6, group = Species, colour = Species)) +
    geom_line(linewidth = 0.55) + geom_point(size = 1) +
    facet_wrap(~ area_type, scales = "free_y", ncol = 1) +
    labs(x = NULL, y = expression("Niche area (10"^6*" km"^2*")"), colour = "Species") +
    theme_manuscript +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
  save_gg(p, file.path(figure_dir, "niche", "species_area.png"), 10.5, 8)
})

run_figure("population_niche_summary", {
  population <- keep_multiforest(read_required(niche_table_file("population")))
  population[, scenario_label := factor(scenario_label(scenario), levels = scenario_label(scenarios))]
  summary <- population[, .(
    mean_area_km2 = mean(suitable_area_km2, na.rm = TRUE),
    se_area_km2 = sd(suitable_area_km2, na.rm = TRUE) / sqrt(.N),
    mean_suitability = mean(mean_dual_suitability, na.rm = TRUE),
    se_suitability = sd(mean_dual_suitability, na.rm = TRUE) / sqrt(.N)
  ), by = .(Species, scenario_label)]
  p_area <- ggplot(
    summary,
    aes(scenario_label, mean_area_km2 / 1e6, group = Species, colour = Species)
  ) +
    geom_line() + geom_point(size = 1) +
    geom_errorbar(aes(
      ymin = pmax(0, mean_area_km2 - se_area_km2) / 1e6,
      ymax = (mean_area_km2 + se_area_km2) / 1e6
    ), width = 0.08, linewidth = 0.3) +
    labs(x = NULL, y = expression("Mean population area (10"^6*" km"^2*")")) +
    theme_manuscript + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  p_suit <- ggplot(
    summary,
    aes(scenario_label, mean_suitability, group = Species, colour = Species)
  ) +
    geom_line() + geom_point(size = 1) +
    geom_errorbar(aes(
      ymin = pmax(0, mean_suitability - se_suitability),
      ymax = pmin(1, mean_suitability + se_suitability)
    ), width = 0.08, linewidth = 0.3) +
    labs(x = NULL, y = "Mean dual suitability", colour = "Species") +
    theme_manuscript + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  combined <- p_area / p_suit + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  save_gg(combined, file.path(figure_dir, "niche", "population_summary.png"), 11, 9)
  
  for (species_name in sort(unique(population$Species))) {
    detail <- population[Species == species_name]
    p_detail <- ggplot(
      detail,
      aes(scenario_label, suitable_area_km2 / 1e6, group = factor(source_zone), colour = factor(source_zone))
    ) +
      geom_line() + geom_point(size = 1) +
      scale_colour_manual(values = zone_colors, name = "Source zone") +
      labs(
        title = species_name, x = NULL,
        y = expression("Population niche area (10"^6*" km"^2*")")
      ) + theme_manuscript +
      theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
    save_gg(
      p_detail,
      file.path(figure_dir, "niche", paste0("population_", safe_name(species_name), ".png")),
      9, 5.5
    )
  }
})

run_figure("species_niche_maps", {
  species_area <- keep_multiforest(read_required(niche_table_file("species")))
  species_names <- sort(unique(species_area$Species))
  available_scenarios <- intersect(scenarios, unique(species_area$scenario))
  for (scenario in available_scenarios) {
    scenario_value <- scenario
    panels <- lapply(species_names, function(species_name) {
      row <- species_area[
        Species == species_name & scenario == scenario_value
      ][1]
      indexed <- if ("species_raster" %in% names(row)) {
        as.character(row[["species_raster"]])
      } else {
        character()
      }
      current_safe <- script6_safe_name(species_name)
      legacy_safe <- script6_legacy_safe_name(species_name)
      file <- first_existing(
        c(
          indexed,
          file.path(
            base_dir, "future tree niche dual suitability var", "mf_var", scenario,
            "species niche", paste0(current_safe, "_species_dual_suitability.tif")
          ),
          file.path(
            base_dir, "future tree niche dual suitability var", "mf_var", scenario,
            "species niche", paste0(legacy_safe, "_species_dual_suitability.tif")
          ),
          file.path(
            base_dir, "population_species", "rasters", scenario,
            paste0(safe_name(species_name), "_species_suitability.tif")
          )
        ),
        paste("species niche", species_name, scenario)
      )
      continuous_map_plot(
        file,
        species_name,
        limits = c(0, 1), max_cells = 55000L
      )
    })
    combined <- wrap_plots(panels, ncol = 5, guides = "collect") +
      plot_annotation(
        title = paste0(
          scenario_label(scenario), ": scalar gated dual suitability"
        ),
        caption = paste(
          "This scalar retains climate suitability only where soil suitability",
          "exceeds 0.2; the bivariate companion separates the two components."
        )
      ) &
      theme(legend.position = "right")
    save_gg(
      combined,
      file.path(
        figure_dir, "niche",
        paste0("species_maps_", scenario, ".png")
      ),
      16, 7.5
    )
  }
})

run_figure("species_bivariate_maps", {
  if (!length(bivariate_scenarios)) {
    skip_figure("Set ECOCHINA2_BIVARIATE_SCENARIOS to draw bivariate species maps.")
  } else {
    species_area <- keep_multiforest(read_required(niche_table_file("species")))
    available_scenarios <- intersect(
      bivariate_scenarios, unique(as.character(species_area$scenario))
    )
    if (!length(available_scenarios)) {
      skip_figure("No requested bivariate scenario is present in the species table.")
    } else for (scenario in available_scenarios) {
      scenario_value <- scenario
      scenario_species <- sort(unique(
        species_area[scenario == scenario_value]$Species
      ))
      component_cache <- new.env(parent = emptyenv())
      panels <- lapply(scenario_species, function(species_name) {
        row <- species_area[
          scenario == scenario_value & Species == species_name
        ][1]
        indexed_rank <- if ("ranked_zone_raster" %in% names(row)) {
          as.character(row[["ranked_zone_raster"]])
        } else {
          character()
        }
        indexed_binary <- if ("species_binary_raster" %in% names(row)) {
          as.character(row[["species_binary_raster"]])
        } else {
          character()
        }
        safe_species <- script6_safe_name(species_name)
        legacy_safe_species <- script6_legacy_safe_name(species_name)
        rank_file <- first_existing(
          c(
            indexed_rank,
            file.path(
              base_dir, "future tree niche dual suitability var", "mf_var",
              scenario, "population ranking",
              paste0(safe_species, "_ranked_population_zones.tif")
            ),
            file.path(
              base_dir, "future tree niche dual suitability var", "mf_var",
              scenario, "population ranking",
              paste0(legacy_safe_species, "_ranked_population_zones.tif")
            )
          ),
          paste("dual-winning population raster", species_name, scenario)
        )
        binary_file <- first_existing(
          c(
            indexed_binary,
            file.path(
              base_dir, "future tree niche dual suitability var", "mf_var",
              scenario, "species niche binary",
              paste0(
                safe_species, "_species_dual_binary_threshold",
                dual_threshold, ".tif"
              )
            ),
            file.path(
              base_dir, "future tree niche dual suitability var", "mf_var",
              scenario, "species niche binary",
              paste0(
                legacy_safe_species, "_species_dual_binary_threshold",
                dual_threshold, ".tif"
              )
            )
          ),
          paste("species dual-niche mask", species_name, scenario)
        )
        map_data <- bivariate_selector_data(
          rank_file,
          scenario,
          mask_file = binary_file,
          max_cells = 55000L,
          cache = component_cache
        )
        bivariate_map_plot(map_data, species_name)
      })
      map_grid <- wrap_plots(panels, ncol = 5)
      combined <- wrap_plots(
        list(map_grid, bivariate_key_plot()),
        ncol = 2, widths = c(7, 1.25)
      ) +
        plot_annotation(
          title = paste0(
            scenario_label(scenario),
            ": species dual niches decomposed into climate \u00d7 soil support"
          ),
          caption = paste0(
            "Only cells with species dual suitability >= ", dual_threshold,
            " are coloured. At each cell, climate and soil belong to the same ",
            "dual-winning source population; maxima are never combined across ",
            "different populations. Raster aggregation is for display only."
          )
        )
      save_gg(
        combined,
        file.path(
          figure_dir, "niche",
          paste0("species_climate_soil_bivariate_", scenario, ".png")
        ),
        18, 8.5
      )
      rm(component_cache, panels, map_grid, combined)
      gc()
    }
  }
})

run_figure("population_bivariate_maps", {
  if (!length(bivariate_scenarios)) {
    skip_figure("Set ECOCHINA2_BIVARIATE_SCENARIOS to draw bivariate population maps.")
  } else {
    population <- keep_multiforest(read_required(niche_table_file("population")))
    required_population_columns <- c(
      "scenario", "Species", "PopulationID", "source_zone"
    )
    missing_population_columns <- setdiff(
      required_population_columns, names(population)
    )
    if (length(missing_population_columns)) {
      stop(
        "Population niche table is missing: ",
        paste(missing_population_columns, collapse = ", ")
      )
    }
    available_scenarios <- intersect(
      bivariate_scenarios, unique(as.character(population$scenario))
    )
    if (!length(available_scenarios)) {
      skip_figure("No requested bivariate scenario is present in the population table.")
    } else for (scenario in available_scenarios) {
      scenario_value <- scenario
      scenario_population <- unique(
        population[scenario == scenario_value, .(
          Species, PopulationID, source_zone = as.integer(source_zone)
        )],
        by = c("Species", "PopulationID", "source_zone")
      )
      component_cache <- new.env(parent = emptyenv())
      for (species_name in sort(unique(scenario_population$Species))) {
        detail <- scenario_population[Species == species_name][
          order(source_zone, PopulationID)
        ]
        panels <- lapply(seq_len(nrow(detail)), function(index) {
          source_zone <- detail$source_zone[[index]]
          map_data <- bivariate_pair_data(
            scenario,
            source_zone,
            max_cells = 55000L,
            cache = component_cache
          )
          bivariate_map_plot(
            map_data,
            paste0(
              detail$PopulationID[[index]], "\nSource zone ", source_zone
            )
          )
        })
        map_grid <- wrap_plots(panels, ncol = min(4L, nrow(detail)))
        combined <- wrap_plots(
          list(map_grid, bivariate_key_plot()),
          ncol = 2, widths = c(7, 1.25)
        ) +
          plot_annotation(
            title = paste0(
              species_name, " | ", scenario_label(scenario),
              ": population climate \u00d7 soil support"
            ),
            caption = paste0(
              "Each population inherits the climate and soil suitability ",
              "surfaces of its source ecotype. Blue shows climatic opportunity ",
              "with weak soil support; green shows joint support. These component ",
              "maps are not thresholded binary niches. Raster aggregation is for ",
              "display only."
            )
          )
        save_gg(
          combined,
          file.path(
            figure_dir, "niche",
            paste0(
              "population_climate_soil_bivariate_", safe_name(species_name),
              "_", scenario, ".png"
            )
          ),
          18,
          max(7, 3.4 * ceiling(nrow(detail) / 4) + 1.5)
        )
        rm(panels, map_grid, combined)
        gc()
      }
      rm(component_cache)
      gc()
    }
  }
})

run_figure("population_overlay_maps", {
  species_area <- keep_multiforest(read_required(niche_table_file("species")))
  species_names <- sort(unique(species_area$Species))
  available_scenarios <- intersect(scenarios, unique(species_area$scenario))
  for (scenario in available_scenarios) {
    scenario_value <- scenario
    panels <- lapply(species_names, function(species_name) {
      row <- species_area[Species == species_name & scenario == scenario_value][1]
      indexed_zone <- if ("ranked_zone_raster" %in% names(row)) {
        as.character(row[["ranked_zone_raster"]])
      } else {
        character()
      }
      current_safe <- script6_safe_name(species_name)
      legacy_safe <- script6_legacy_safe_name(species_name)
      combined_file <- file.path(
        base_dir, "population_species", "rasters", scenario,
        paste0(safe_name(species_name), "_top_population.tif")
      )
      display <- first_existing(
        c(
          indexed_zone,
          file.path(
            base_dir, "future tree niche dual suitability var", "mf_var",
            scenario, "population ranking",
            paste0(current_safe, "_ranked_population_zones.tif")
          ),
          file.path(
            base_dir, "future tree niche dual suitability var", "mf_var",
            scenario, "population ranking",
            paste0(legacy_safe, "_ranked_population_zones.tif")
          ),
          combined_file
        ),
        paste("dual-winning population map", species_name, scenario)
      )
      zone_map_plot(
        display,
        species_name,
        max_cells = 50000L
      )
    })
    combined <- wrap_plots(panels, ncol = 5, guides = "collect") +
      plot_annotation(
        title = paste0(
          scenario_label(scenario), ": dual-winning source population"
        ),
        caption = paste(
          "Colour identifies the source ecotype of the dual-winning population;",
          "the companion bivariate figure shows its climate x soil support."
        )
      ) &
      theme(legend.position = "right")
    save_gg(
      combined,
      file.path(
        figure_dir, "niche",
        paste0("population_overlay_", scenario, ".png")
      ),
      17, 8
    )
  }
})

# 6. Top-k diagnostics (k = 1, 2, 3, 4, 5) ====================================

topk_table_file <- function(name, required = TRUE) {
  first_existing(
    c(
      file.path(base_dir, "assessment", "topk", "tables", name),
      file.path(base_dir, "assessment", "topk", name),
      file.path(base_dir, "topk", "tables", name)
    ),
    paste("Top-k table", name),
    required = required
  )
}
topk_raster_dir <- first_existing(
  c(
    file.path(base_dir, "dual suit ranking var", "mf_var"),
    file.path(base_dir, "topk", "rasters")
  ),
  "Top-k raster directory",
  required = FALSE
)
if (is.na(topk_raster_dir)) {
  topk_raster_dir <- file.path(base_dir, "dual suit ranking var", "mf_var")
}

reference_topk_map_file <- function(k) {
  if (k == 1L) return(binary_map_file("normal"))
  normal_assigned_file <- binary_map_file("normal")
  normal_rank_file <- required_file(
    file.path(topk_raster_dir, "normal", "ranked_zone.tif"),
    "reference ranked-zone raster"
  )
  reference_inputs <- c(reference_file, normal_assigned_file, normal_rank_file)
  first_current_raster(
    c(
      file.path(base_dir, "assessment_var", "future_topk_matched", "cache", "reference_topk_maps", paste0("mf_var_normal_top", k, "_common_domain_v2.tif")),
      file.path(base_dir, "assessment", "topk", "cache", "reference_topk_maps", paste0("mf_var_normal_top", k, "_common_domain_v2.tif"))
    ),
    reference_file,
    reference_inputs,
    paste0("reference-conditioned Top-", k, " diagnostic map")
  )
}

future_topk_map_file <- function(scenario, k) {
  if (k == 1L) {
    reference <- rast(reference_file)
    normal <- rast(binary_map_file("normal"))
    future <- rast(binary_map_file(scenario))
    if (!compareGeom(reference, normal, stopOnError = FALSE) ||
        !compareGeom(reference, future, stopOnError = FALSE)) {
      stop("Reference, normal and future Top-1 maps have incompatible geometry: ", scenario)
    }
    return(ifel(
      reference %in% modeled_zones & !is.na(normal) & !is.na(future),
      future,
      NA
    ))
  }
  normal_assigned_file <- binary_map_file("normal")
  future_assigned_file <- binary_map_file(scenario)
  future_rank_file <- required_file(
    file.path(topk_raster_dir, scenario, "ranked_zone.tif"),
    paste("future ranked-zone raster", scenario)
  )
  future_inputs <- c(
    reference_file, normal_assigned_file, future_assigned_file, future_rank_file
  )
  first_current_raster(
    c(
      file.path(base_dir, "assessment", "topk", "cache", "future_topk_maps", scenario, paste0("future_top", k, ".tif")),
      file.path(topk_raster_dir, scenario, paste0("future_top", k, ".tif"))
    ),
    reference_file,
    future_inputs,
    paste0(
      "future reference-conditioned Top-", k,
      " diagnostic map | ", scenario
    )
  )
}

run_figure("topk_reference_agreement", {
  agreement <- read_required(topk_table_file("reference_agreement.csv"))
  long <- melt(
    agreement, id.vars = "k",
    measure.vars = c("pixel_agreement", "area_agreement"),
    variable.name = "basis", value.name = "agreement"
  )
  long[, basis := factor(
    basis, levels = c("pixel_agreement", "area_agreement"),
    labels = c("Pixel", "Area")
  )]
  p <- ggplot(long, aes(k, agreement, colour = basis)) +
    geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
    scale_x_continuous(breaks = k_values) +
    scale_y_continuous(labels = scales::percent) +
    scale_colour_brewer(palette = "Dark2") +
    labs(x = "Top-k", y = "Observed ecotype recovered", colour = "Basis") +
    theme_manuscript + theme(legend.position = "top")
  save_gg(p, file.path(figure_dir, "topk", "reference_agreement.png"), 7, 4.8)
})

run_figure("topk_reference_confusion", {
  confusion <- read_required(topk_table_file("reference_confusion.csv"))
  confusion[, proportion := pixels / sum(pixels), by = .(k, observed_zone)]
  confusion[, topk_label := factor(
    paste0("Top-", k), levels = paste0("Top-", k_values)
  )]
  p <- ggplot(
    confusion,
    aes(factor(assigned_zone), factor(observed_zone), fill = proportion)
  ) +
    geom_tile() + facet_wrap(~ topk_label, ncol = 5) +
    scale_fill_viridis_c(option = "C", trans = "sqrt", limits = c(0, 1)) +
    coord_equal() +
    labs(x = "Diagnostic assignment", y = "Observed zone", fill = "Row share") +
    theme_manuscript + theme(axis.text = element_text(size = 3.5))
  save_gg(p, file.path(figure_dir, "topk", "reference_confusion.png"), 18, 5.2)
})

run_figure("topk_future_retention", {
  retention <- read_required(topk_table_file("future_retention.csv"))
  retention <- retention[retention_class == "cumulative_topk" & !is.na(k)]
  retention[, scenario_label := factor(
    scenario_label(scenario), levels = scenario_label(future_scenarios)
  )]
  p <- ggplot(retention, aes(k, retained_area_share, colour = scenario_label)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.8) +
    scale_x_continuous(breaks = k_values) +
    scale_y_continuous(labels = scales::percent) +
    scale_colour_viridis_d(option = "D") +
    labs(x = "Top-k", y = "Reference Top-1 ecotype retained", colour = "Scenario") +
    theme_manuscript
  save_gg(p, file.path(figure_dir, "topk", "future_retention.png"), 9, 5.8)
})

run_figure("topk_zone_diagnostics", {
  by_zone <- read_required(topk_table_file("future_retention_by_zone.csv"))
  totals <- by_zone[, .(total_area_km2 = sum(area_km2)), by = .(scenario, k, normal_zone)]
  retained <- merge(
    by_zone[status == "retained", .(retained_area_km2 = sum(area_km2)),
            by = .(scenario, k, normal_zone)],
    totals,
    by = c("scenario", "k", "normal_zone"),
    all = TRUE
  )
  retained[is.na(retained_area_km2), retained_area_km2 := 0]
  retained[, retained_share := retained_area_km2 / total_area_km2]
  retained[, scenario_label := factor(
    scenario_label(scenario), levels = scenario_label(future_scenarios)
  )]
  retained[, topk_label := factor(
    paste0("Top-", k), levels = paste0("Top-", k_values)
  )]
  p <- ggplot(
    retained,
    aes(scenario_label, factor(normal_zone), fill = retained_share)
  ) +
    geom_tile() + facet_wrap(~ topk_label, nrow = 1) +
    scale_fill_viridis_c(option = "C", limits = c(0, 1), labels = scales::percent) +
    labs(x = NULL, y = "Reference Top-1 zone", fill = "Retained\narea share") +
    theme_manuscript +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 6))
  save_gg(p, file.path(figure_dir, "topk", "retention_by_zone.png"), 18, 8)
  
  source_area <- read_required(topk_table_file("source_zone_area.csv"))
  source_area[, scenario_label := factor(
    scenario_label(scenario), levels = scenario_label(scenarios)
  )]
  source_area[, topk_label := factor(
    paste0("Top-", k), levels = paste0("Top-", k_values)
  )]
  p_source <- ggplot(
    source_area,
    aes(scenario_label, factor(source_zone), fill = area_km2 / 1e6)
  ) +
    geom_tile() + facet_wrap(~ topk_label, nrow = 1) +
    scale_fill_viridis_c(option = "C") +
    labs(
      x = NULL, y = "Source ecotype zone",
      fill = expression("Represented area (10"^6*" km"^2*")")
    ) + theme_manuscript +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 6))
  save_gg(p_source, file.path(figure_dir, "topk", "source_zone_area.png"), 18, 8)
})

run_figure("topk_matched_change", {
  change <- read_required(topk_table_file("matched_change.csv"))
  change[, share := area_km2 / sum(area_km2), by = .(scenario, k)]
  change[, scenario_label := factor(
    scenario_label(scenario), levels = scenario_label(future_scenarios)
  )]
  change[, change_class := factor(change_class, c("stable", "changed", "novel"))]
  p <- ggplot(change, aes(factor(k), share, fill = change_class)) +
    geom_col(width = 0.72) +
    facet_wrap(~ scenario_label, ncol = 3) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_discrete(labels = paste0("Top-", k_values)) +
    scale_fill_manual(values = c(stable = "#3B8C6E", changed = "#E3A43B", novel = "#252525")) +
    labs(
      x = "Same k in reference and future", y = "Area share", fill = "Class",
      caption = "Novel area is defined independently of k (max dual suitability < 0.4)."
    ) + theme_manuscript + theme(legend.position = "top")
  save_gg(p, file.path(figure_dir, "topk", "matched_change.png"), 12, 7)
})

run_figure("topk_analogue_availability", {
  analogue <- read_required(topk_table_file("analogue_availability.csv"))
  analogue[, scenario_label := factor(scenario_label(scenario), levels = scenario_label(scenarios))]
  analogue[, area_million_km2 := area_below_km2 / 1e6]
  p <- ggplot(
    analogue,
    aes(minimum_analogues, area_million_km2, colour = scenario_label)
  ) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.8) +
    scale_x_continuous(breaks = k_values) +
    scale_colour_viridis_d(option = "D") +
    labs(
      x = "Minimum number of suitable current ecotypes",
      y = expression("Area below minimum (10"^6*" km"^2*")"),
      colour = "Scenario"
    ) + theme_manuscript
  save_gg(p, file.path(figure_dir, "topk", "analogue_availability.png"), 9, 5.8)
})

run_figure("topk_species_population", {
  species <- read_required(topk_table_file("species_area.csv"))
  species[, scenario_label := factor(scenario_label(scenario), levels = scenario_label(scenarios))]
  p_species <- ggplot(
    species,
    aes(k, area_km2 / 1e6, colour = Species, group = Species)
  ) +
    geom_line() + geom_point(size = 0.8) +
    facet_wrap(~ scenario_label, ncol = 3, scales = "free_y") +
    scale_x_continuous(breaks = k_values) +
    labs(
      x = "Top-k", y = expression("Species represented area (10"^6*" km"^2*")"),
      colour = "Species"
    ) + theme_manuscript + theme(legend.position = "bottom")
  save_gg(p_species, file.path(figure_dir, "topk", "species_area.png"), 12, 8.5)
  
  population_file <- topk_table_file("population_area.csv", required = FALSE)
  if (!is.na(population_file)) {
    population <- fread(population_file)
    population[, scenario_label := factor(scenario_label(scenario), levels = scenario_label(scenarios))]
    for (species_name in sort(unique(population$Species))) {
      detail <- population[Species == species_name]
      p <- ggplot(
        detail,
        aes(k, area_km2 / 1e6, colour = factor(source_zone), group = factor(source_zone))
      ) +
        geom_line() + geom_point(size = 0.9) +
        facet_wrap(~ scenario_label, ncol = 3, scales = "free_y") +
        scale_x_continuous(breaks = k_values) +
        scale_colour_manual(values = zone_colors) +
        labs(
          title = species_name, x = "Top-k",
          y = expression("Population represented area (10"^6*" km"^2*")"),
          colour = "Source zone"
        ) + theme_manuscript + theme(legend.position = "bottom")
      save_gg(
        p,
        file.path(figure_dir, "topk", paste0("population_", safe_name(species_name), ".png")),
        11, 8
      )
    }
  }
})

run_figure("topk_f1_multiclass", {
  comparison_file <- topk_table_file("f1_vs_multiclass.csv", required = FALSE)
  if (length(comparison_file) != 1L || is.na(comparison_file)) {
    skip_figure(
      "f1_vs_multiclass.csv is optional and is not available until script 7 completes"
    )
  } else {
    comparison <- fread(comparison_file)
    comparison[, topk_label := factor(
      paste0("Top-", k), levels = paste0("Top-", k_values)
    )]
    p <- ggplot(comparison, aes(multiclass_f1, f1, colour = factor(k))) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45") +
      geom_point(aes(size = k), alpha = 0.75) +
      scale_colour_viridis_d(option = "D") +
      facet_wrap(~ topk_label, nrow = 1) +
      labs(
        x = "Multiclass reference-map F1", y = "Binary diagnostic F1",
        colour = "Top-k", size = "Top-k"
      ) + theme_manuscript + theme(legend.position = "none")
    save_gg(p, file.path(figure_dir, "topk", "f1_vs_multiclass.png"), 14, 4.2)
  }
})

run_figure("topk_reference_maps", {
  panels <- lapply(k_values, function(k) {
    zone_map_plot(
      reference_topk_map_file(k),
      if (k == 1L) {
        "Reference-period Top-1 assignment"
      } else {
        paste0("Reference-conditioned Top-", k, " diagnostic")
      }
    )
  })
  combined <- wrap_plots(panels, ncol = 5, guides = "collect") &
    theme(legend.position = "right")
  save_gg(
    combined,
    file.path(figure_dir, "topk", "reference_topk_diagnostics.png"),
    18, 4.5
  )
})

run_figure("topk_future_maps", {
  display_zone_ids <- c(modeled_zones, 99L)
  missing_palette_ids <- setdiff(display_zone_ids, zone_palette$zoneID)
  if (length(missing_palette_ids)) {
    stop(
      "Zone palette is missing Top-k display IDs: ",
      paste(missing_palette_ids, collapse = ", ")
    )
  }
  display_zone_colors <- zone_colors[as.character(display_zone_ids)]
  display_zone_labels <- zone_labels[as.character(display_zone_ids)]
  display_zone_labels[["99"]] <- "Novel ecosystem (99)"
  
  # Use fixed factor levels and a fixed legend in every panel. Each map has its
  # own full plotting area; no pixels are hidden by an exploded-stack offset.
  topk_projection_panel <- function(file, title) {
    raster <- if (inherits(file, "SpatRaster")) file else rast(required_file(file))
    x <- raster_points(raster, TRUE, maximum_map_cells)
    x[, zoneID := as.integer(value)]
    x <- x[zoneID %in% display_zone_ids]
    x[, zone := factor(zoneID, levels = display_zone_ids)]
    ggplot(x, aes(x, y, fill = zone)) +
      geom_raster() +
      scale_fill_manual(
        values = display_zone_colors,
        labels = display_zone_labels,
        limits = as.character(display_zone_ids),
        drop = FALSE,
        na.value = "transparent",
        name = "Ecotype",
        guide = guide_legend(
          ncol = 3, byrow = FALSE,
          override.aes = list(alpha = 1)
        )
      ) +
      coord_equal(expand = FALSE) +
      labs(title = title, x = NULL, y = NULL) +
      theme_void(base_size = 9) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 9),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8)
      )
  }
  
  period_title <- function(scenario) {
    period <- gsub("-", "\u2013", sub("SSP.*$", "", scenario))
    ssp <- sub("^.*(SSP[0-9]+)$", "\\1", scenario)
    paste0(period, "\n", ssp)
  }
  
  ssp245_scenarios <- future_scenarios[grepl("SSP245$", future_scenarios)]
  ssp585_scenarios <- future_scenarios[grepl("SSP585$", future_scenarios)]
  if (length(ssp245_scenarios) != 3L || length(ssp585_scenarios) != 3L) {
    stop("Expected three time periods for each of SSP245 and SSP585.")
  }
  
  reference <- rast(reference_file)
  normal_assigned <- rast(binary_map_file("normal"))
  if (!compareGeom(reference, normal_assigned, stopOnError = FALSE)) {
    stop("Observed reference and reference-period assigned map have incompatible geometry.")
  }
  observed_common <- ifel(
    reference %in% modeled_zones & !is.na(normal_assigned),
    reference,
    NA
  )
  observed_panel <- topk_projection_panel(
    observed_common,
    "Observed ecotypes\n1961\u20131990"
  )
  
  projection_outputs <- file.path(
    figure_dir, "topk", paste0("future_topk_k", k_values, ".png")
  )
  
  for (k in k_values) {
    future_panels <- setNames(lapply(future_scenarios, function(scenario) {
      topk_projection_panel(
        future_topk_map_file(scenario, k),
        period_title(scenario)
      )
    }), future_scenarios)
    reference_panel <- topk_projection_panel(
      reference_topk_map_file(k),
      if (k == 1L) {
        "Reference-period\nTop-1 assignment"
      } else {
        paste0("Reference diagnostic\nTop-", k)
      }
    )
    panels <- c(
      list(observed_panel),
      future_panels[ssp245_scenarios],
      list(reference_panel),
      future_panels[ssp585_scenarios]
    )
    
    if (k == 1L) {
      figure_title <- "Top-1 assigned ecosystem projection"
      figure_subtitle <-
        "Each future pixel shows its highest dual-suitability ecotype."
      figure_caption <- paste0(
        "Zone 99 denotes a novel ecosystem (maximum dual suitability < ",
        dual_threshold, ")."
      )
    } else {
      figure_title <- paste0(
        "Reference-conditioned Top-", k, " ecosystem diagnostic"
      )
      figure_subtitle <- paste0(
        "The observed ecotype is retained when it occurs among the first ",
        k, " future dual-suitability ranks; otherwise the future Top-1 ",
        "assignment is shown."
      )
      figure_caption <- paste0(
        "Top-", k, " is a reference-conditioned sensitivity diagnostic, ",
        "not an independent categorical forecast. Zone 99 denotes a novel ",
        "ecosystem (maximum dual suitability < ", dual_threshold,
        ") and is independent of k."
      )
    }
    
    combined <- wrap_plots(panels, ncol = 4, guides = "collect") +
      plot_annotation(
        title = figure_title,
        subtitle = figure_subtitle,
        caption = figure_caption
      ) &
      theme(legend.position = "right")
    save_gg(
      combined,
      projection_outputs[[k]],
      22, 10.5
    )
    rm(future_panels, reference_panel, panels, combined)
    gc()
  }
  projection_info <- file.info(projection_outputs)
  if (any(!file.exists(projection_outputs)) ||
      any(is.na(projection_info$size) | projection_info$size <= 0)) {
    stop(
      "One or more Top-k future figures were not written:\n",
      paste(projection_outputs, collapse = "\n")
    )
  }
  projection_manifest <- data.table(
    k = k_values,
    figure_type = c(
      "unconditional_top1_projection",
      rep("reference_conditioned_topk_diagnostic", 4L)
    ),
    file = normalizePath(
      projection_outputs, winslash = "/", mustWork = TRUE
    ),
    size_bytes = as.numeric(projection_info$size),
    modified = as.character(projection_info$mtime),
    script_revision = script_revision
  )
  fwrite(
    projection_manifest,
    file.path(figure_dir, "topk", "future_topk_manifest.csv")
  )
  cat(
    "[TOP-K FUTURE FIGURES COMPLETE]\n",
    paste(projection_manifest$file, collapse = "\n"), "\n",
    sep = ""
  )
  rm(reference, normal_assigned, observed_common, observed_panel)
  gc()
})


# 7. Top-5 candidate-layer diagnostics =========================================

stack_setting <- Sys.getenv(
  "ECOCHINA2_STACK_SCENARIOS",
  unset = ""
)
stack_scenarios <- trimws(strsplit(stack_setting, ",", fixed = TRUE)[[1]])
stack_scenarios <- stack_scenarios[nzchar(stack_scenarios)]
unknown_stack <- setdiff(stack_scenarios, scenarios)
if (length(unknown_stack)) {
  stop("Unknown ECOCHINA2_STACK_SCENARIOS: ", paste(unknown_stack, collapse = ", "))
}

# Figure IDs are retained for console compatibility. Visible labels describe
# each ordered candidate as part of the Top-5 set, while Top-k remains cumulative.
run_figure("topk_rank_panels", {
  if (!length(stack_scenarios)) {
    skip_figure(
      "Set ECOCHINA2_STACK_SCENARIOS to draw Top-5 candidate details."
    )
  } else {
    for (scenario in stack_scenarios) {
      zone_stack <- rast(required_file(
        file.path(topk_raster_dir, scenario, "ranked_zone.tif")
      ))
      component_cache <- new.env(parent = emptyenv())
      panels <- list()
      for (rank_value in k_values) {
        panels[[length(panels) + 1L]] <- zone_map_plot(
          zone_stack[[rank_value]],
          paste0("Top-5 candidate ", rank_value, ": ecotype")
        )
        component_data <- bivariate_selector_data(
          zone_stack[[rank_value]],
          scenario,
          max_cells = 70000L,
          cache = component_cache
        )
        panels[[length(panels) + 1L]] <- bivariate_map_plot(
          component_data,
          paste0("Top-5 candidate ", rank_value, ": climate \u00d7 soil")
        )
      }
      map_grid <- wrap_plots(panels, ncol = 2, guides = "collect") &
        theme(legend.position = "right")
      combined <- wrap_plots(
        list(map_grid, bivariate_key_plot()),
        ncol = 2, widths = c(7, 1.2)
      ) +
        plot_annotation(
          title = paste0(
            scenario_label(scenario),
            ": Top-5 candidate ecotypes and their climate \u00d7 soil support"
          ),
          caption = paste0(
            "Each row shows one ordered candidate within the Top-5 set. A Top-k ",
            "result is cumulative over ordered candidates 1 through k. The right panel reads ",
            "climate and soil suitability for the same candidate ecotype; it does ",
            "not recolour the one-dimensional gated dual value. Raster aggregation ",
            "is for display only."
          )
        )
      save_gg(
        combined,
        file.path(
          figure_dir, "topk",
          paste0("rank_panel_", scenario, ".png")
        ),
        15, 20
      )
      rm(zone_stack, component_cache, panels, map_grid, combined)
      gc()
    }
  }
})

run_figure("topk_exploded_stack", {
  if (!length(stack_scenarios)) {
    skip_figure(
      "Set ECOCHINA2_STACK_SCENARIOS to draw Top-5 candidate maps."
    )
  } else {
    for (scenario in stack_scenarios) {
      zone_stack <- rast(required_file(
        file.path(topk_raster_dir, scenario, "ranked_zone.tif")
      ))
      candidates <- rbindlist(lapply(k_values, function(rank_value) {
        x <- raster_points(zone_stack[[rank_value]], TRUE, 70000L)
        x[, `:=`(
          candidate_rank = rank_value,
          zone = factor(as.integer(value), levels = zone_palette$zoneID)
        )]
        x
      }))
      candidates[, candidate_label := factor(
        candidate_rank,
        levels = k_values,
        labels = paste0("Candidate ", k_values, " of Top-5")
      )]
      
      p <- ggplot(
        candidates,
        aes(x = x, y = y, fill = zone)
      ) +
        geom_raster() +
        facet_wrap(~ candidate_label, ncol = 2, drop = FALSE) +
        scale_fill_manual(
          values = zone_colors, labels = zone_labels, drop = TRUE,
          name = "Ecotype zone",
          guide = guide_legend(ncol = 3)
        ) +
        coord_equal(expand = FALSE) +
        labs(
          title = paste0(
            scenario_label(scenario), ": Top-5 candidate ecotype layers"
          ),
          subtitle = paste(
            "Each ordered candidate is a complete categorical map;",
            "no raster cells are hidden by another layer."
          ),
          caption = paste(
            "Each panel shows one ordered candidate within the Top-5 set.",
            "Top-k results are cumulative over ordered candidates 1 through k."
          ),
          x = NULL, y = NULL
        ) + theme_void(base_size = 10) +
        theme(
          legend.position = "right",
          strip.text = element_text(size = 10, face = "bold"),
          strip.background = element_rect(fill = "grey95", colour = NA),
          panel.spacing = grid::unit(5, "mm"),
          plot.title = element_text(size = 13, face = "bold"),
          plot.subtitle = element_text(size = 10)
        )
      save_gg(
        p,
        file.path(
          figure_dir, "topk",
          paste0("rank_exploded_", scenario, ".png")
        ),
        15, 13
      )
      rm(zone_stack, candidates, p)
      gc()
    }
  }
})


# 8. Audit =====================================================================

figure_audit <- rbindlist(figure_log, fill = TRUE)
fwrite(figure_audit, file.path(figure_dir, "figure_log.csv"))

if (!identical(selected_figure_ids, "all")) {
  unknown_figure_ids <- setdiff(selected_figure_ids, unique(figure_audit$figure))
  if (length(unknown_figure_ids)) {
    stop(
      "Unknown ECOCHINA2_FIGURE_IDS: ",
      paste(unknown_figure_ids, collapse = ", ")
    )
  }
}

failed <- figure_audit[status == "failed"]
missing_inputs <- figure_audit[status == "missing_input"]
fwrite(
  figure_audit[status %in% c("failed", "missing_input")],
  file.path(figure_dir, "figure_pending.csv")
)
cat(
  "\nFigure families completed: ", sum(figure_audit$status == "complete"),
  "; skipped: ", sum(figure_audit$status == "skipped"),
  "; missing input: ", nrow(missing_inputs),
  "; failed: ", nrow(failed), "\n",
  sep = ""
)
if (nrow(missing_inputs)) {
  cat(
    "Pending input (see figures/figure_pending.csv): ",
    paste(missing_inputs$figure, collapse = ", "), "\n", sep = ""
  )
}
if (nrow(failed)) {
  stop(
    nrow(failed), " figure family/families failed. See figures/figure_log.csv: ",
    paste(failed$figure, collapse = ", ")
  )
}

cat(
  if (nrow(missing_inputs)) "\nVISUALIZATION PARTIAL - INPUTS PENDING\n" else "\nVISUALIZATION COMPLETE\n",
  "Figure families: ", nrow(figure_audit), "\n",
  "Output: ", figure_dir, "\n",
  sep = ""
)
