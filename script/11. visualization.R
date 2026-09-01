# Centralized manuscript visualization for the clean Multi-Forest workflow
# ==============================================================================
# Run after scripts 5, 5.3, 6, 7 and 8. Analysis tables are never recalculated
# here; this script turns their canonical outputs into manuscript-ready figures.
# Set ECOCHINA2_STACK_SCENARIOS to a comma-separated subset of scenarios when
# choosing the detailed Top-k raster-stack figures.

library(data.table)
library(terra)
library(ggplot2)
library(patchwork)

rm(list = ls())
gc()


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
  ggsave(
    filename = filename, plot = plot, width = width, height = height,
    units = "in", dpi = dpi, bg = "white"
  )
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
  started <- Sys.time()
  result <- tryCatch({
    force(expression)
    data.table(
      figure = id, section = section, status = "complete", message = NA_character_,
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
    facet_wrap(~ category_label, scales = "free_x", ncol = 2) +
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

run_figure("assessment_model_performance", {
  metrics <- read_required(file.path(assessment_dir, "rf_test_zone_metrics.csv"))
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
  confusion <- read_required(file.path(assessment_dir, "normal_map_confusion_long.csv"))
  confusion <- row_proportion(confusion, "original_zone", "pixels")
  p <- ggplot(confusion, aes(factor(assigned_zone), factor(original_zone), fill = proportion)) +
    geom_tile() +
    scale_fill_viridis_c(option = "C", trans = "sqrt", labels = scales::percent) +
    coord_equal() +
    labs(x = "Assigned zone", y = "Observed zone", fill = "Row share") +
    theme_manuscript + theme(axis.text = element_text(size = 5))
  save_gg(p, file.path(figure_dir, "assessment", "reference_confusion.png"), 9, 8)

  category <- read_required(
    file.path(assessment_dir, "normal_map_category_confusion_long.csv")
  )
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
  binary <- read_required(file.path(assessment_dir, "normal_map_zone_metrics.csv"))
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

  overall <- read_required(file.path(assessment_dir, "normal_map_overall_metrics.csv"))
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
  zone_flow <- read_required(file.path(assessment_dir, "normal_map_confusion_long.csv"))
  make_chord(
    zone_flow, "original_zone", "assigned_zone", "pixels",
    file.path(figure_dir, "assessment", "reference_zone_chord.pdf"),
    zone_colors
  )
  category_flow <- read_required(
    file.path(assessment_dir, "normal_map_category_confusion_long.csv")
  )
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
  gate <- read_required(file.path(assessment_dir, "soil_gate_sensitivity.csv"))
  gate[, scenario_label := factor(scenario_label(scenario), levels = scenario_label(scenarios))]
  p <- ggplot(gate, aes(soil_gate, below_0.4_share, colour = scenario_label)) +
    geom_vline(xintercept = 0.2, linetype = 2, colour = "grey35") +
    geom_line(linewidth = 0.75) + geom_point(size = 1.5) +
    scale_x_continuous(breaks = unique(gate$soil_gate)) +
    scale_y_continuous(labels = scales::percent) +
    scale_colour_viridis_d(option = "D", end = 0.9) +
    labs(
      x = "Soil gate", y = "Area below dual-suitability 0.4",
      colour = "Scenario",
      caption = "Gate 0 is climate-only; the dashed line is the primary gate (0.2)."
    ) + theme_manuscript + theme(legend.position = "right")
  save_gg(p, file.path(figure_dir, "assessment", "soil_gate_sensitivity.png"), 9.5, 6)
})


# 3. Feature importance =========================================================

importance_table_dir <- file.path(assessment_dir, "feature_importance", "tables")

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

run_figure("importance_category_groups", {
  grouped <- read_required(
    file.path(importance_table_dir, "FI_reader_06_category_grouped_importance.csv")
  )
  p <- ggplot(grouped, aes(variable_group, category2, fill = consensus_mean)) +
    geom_tile(colour = "white", linewidth = 0.15) +
    facet_wrap(~ niche, scales = "free_x", ncol = 1) +
    scale_fill_viridis_c(option = "C", labels = scales::percent) +
    labs(
      x = "Predictor group", y = "Vegetation category",
      fill = "Consensus\nimportance"
    ) + theme_manuscript +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
  save_gg(p, file.path(figure_dir, "importance", "category_predictor_groups.png"), 12, 9)

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

run_figure("importance_category_variables", {
  variables <- read_required(
    file.path(importance_table_dir, "FI_reader_07_category_variable_importance.csv")
  )
  top <- variables[order(category2, niche, consensus_rank)][consensus_rank <= 8]
  top[, variable_label := factor(variable, levels = rev(unique(variable)))]
  p <- ggplot(top, aes(variable_label, category2, size = consensus_mean, colour = variable_group)) +
    geom_point(alpha = 0.85) + facet_wrap(~ niche, scales = "free_x", ncol = 1) +
    scale_size_area(max_size = 7, labels = scales::percent) +
    coord_cartesian(clip = "off") +
    labs(
      x = "Top category-level variables", y = "Vegetation category",
      size = "Consensus\nimportance", colour = "Predictor group"
    ) + theme_manuscript +
    theme(axis.text.x = element_text(angle = 55, hjust = 1), legend.position = "bottom")
  save_gg(p, file.path(figure_dir, "importance", "category_variables.png"), 13, 10)
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
      indexed <- species_area[Species == species_name & scenario == scenario_value, species_raster]
      file <- first_existing(
        c(
          indexed,
          file.path(
            base_dir, "future tree niche dual suitability var", "mf_var", scenario,
            "species niche", paste0(safe_name(species_name), "_species_dual_suitability.tif")
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
      plot_annotation(title = scenario_label(scenario)) &
      theme(legend.position = "right")
    save_gg(
      combined,
      file.path(figure_dir, "niche", paste0("species_maps_", scenario, ".png")),
      16, 7.5
    )
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
      legacy_zone <- if ("ranked_zone_raster" %in% names(row)) row$ranked_zone_raster else NA_character_
      legacy_suit <- if ("ranked_suitability_raster" %in% names(row)) row$ranked_suitability_raster else NA_character_
      combined_file <- file.path(
        base_dir, "population_species", "rasters", scenario,
        paste0(safe_name(species_name), "_top_population.tif")
      )
      display <- if (!is.na(legacy_zone) && !is.na(legacy_suit) &&
        file.exists(legacy_zone) && file.exists(legacy_suit)) {
        c(rast(legacy_zone)[[1]], rast(legacy_suit)[[1]])
      } else {
        combined_file
      }
      zone_suitability_plot(
        display,
        species_name,
        max_cells = 50000L
      )
    })
    combined <- wrap_plots(panels, ncol = 5, guides = "collect") +
      plot_annotation(title = paste0(scenario_label(scenario), ": dominant source population")) &
      theme(legend.position = "right")
    save_gg(
      combined,
      file.path(figure_dir, "niche", paste0("population_overlay_", scenario, ".png")),
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
    paste("reference Top-k map", k)
  )
}

future_topk_map_file <- function(scenario, k) {
  if (k == 1L) {
    normal <- rast(binary_map_file("normal"))
    future <- rast(binary_map_file(scenario))
    if (!compareGeom(normal, future, stopOnError = FALSE)) {
      stop("Normal and future Top-1 maps have incompatible geometry: ", scenario)
    }
    return(ifel(!is.na(normal) & !is.na(future), future, NA))
  }
  first_existing(
    c(
      file.path(base_dir, "assessment", "topk", "cache", "future_topk_maps", scenario, paste0("future_top", k, ".tif")),
      file.path(topk_raster_dir, scenario, paste0("future_top", k, ".tif"))
    ),
    paste("future Top-k map", scenario, k)
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
  p <- ggplot(
    confusion,
    aes(factor(assigned_zone), factor(observed_zone), fill = proportion)
  ) +
    geom_tile() + facet_wrap(~ paste0("k = ", k), ncol = 5) +
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
  p <- ggplot(
    retained,
    aes(scenario_label, factor(normal_zone), fill = retained_share)
  ) +
    geom_tile() + facet_wrap(~ paste0("k = ", k), nrow = 1) +
    scale_fill_viridis_c(option = "C", limits = c(0, 1), labels = scales::percent) +
    labs(x = NULL, y = "Reference Top-1 zone", fill = "Retained\narea share") +
    theme_manuscript +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 6))
  save_gg(p, file.path(figure_dir, "topk", "retention_by_zone.png"), 18, 8)

  source_area <- read_required(topk_table_file("source_zone_area.csv"))
  source_area[, scenario_label := factor(
    scenario_label(scenario), levels = scenario_label(scenarios)
  )]
  p_source <- ggplot(
    source_area,
    aes(scenario_label, factor(source_zone), fill = area_km2 / 1e6)
  ) +
    geom_tile() + facet_wrap(~ paste0("k = ", k), nrow = 1) +
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
  comparison <- read_required(topk_table_file("f1_vs_multiclass.csv"))
  p <- ggplot(comparison, aes(multiclass_f1, f1, colour = factor(k))) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45") +
    geom_point(aes(size = k), alpha = 0.75) +
    scale_colour_viridis_d(option = "D") +
    facet_wrap(~ k, nrow = 1) +
    labs(
      x = "Multiclass reference-map F1", y = "Binary diagnostic F1",
      colour = "Top-k", size = "Top-k"
    ) + theme_manuscript + theme(legend.position = "none")
  save_gg(p, file.path(figure_dir, "topk", "f1_vs_multiclass.png"), 14, 4.2)
})

run_figure("topk_reference_maps", {
  panels <- lapply(k_values, function(k) {
    zone_map_plot(
      reference_topk_map_file(k),
      paste0("Reference diagnostic k = ", k)
    )
  })
  combined <- wrap_plots(panels, ncol = 5, guides = "collect") &
    theme(legend.position = "right")
  save_gg(combined, file.path(figure_dir, "topk", "reference_maps.png"), 18, 4.5)
})

run_figure("topk_future_maps", {
  for (scenario in future_scenarios) {
    panels <- lapply(k_values, function(k) {
      zone_map_plot(
        future_topk_map_file(scenario, k),
        paste0("k = ", k)
      )
    })
    combined <- wrap_plots(panels, ncol = 5, guides = "collect") +
      plot_annotation(title = scenario_label(scenario)) &
      theme(legend.position = "right")
    save_gg(
      combined,
      file.path(figure_dir, "topk", paste0("future_maps_", scenario, ".png")),
      18, 4.8
    )
  }
})


# 7. Five-rank raster-stack visualizations =====================================

stack_setting <- Sys.getenv(
  "ECOCHINA2_STACK_SCENARIOS",
  unset = "normal,2071-2100SSP585"
)
stack_scenarios <- trimws(strsplit(stack_setting, ",", fixed = TRUE)[[1]])
unknown_stack <- setdiff(stack_scenarios, scenarios)
if (length(unknown_stack)) {
  stop("Unknown ECOCHINA2_STACK_SCENARIOS: ", paste(unknown_stack, collapse = ", "))
}

rank_pair_data <- function(scenario, rank_value, max_cells = 70000L) {
  zone_stack <- rast(required_file(
    file.path(topk_raster_dir, scenario, "ranked_zone.tif")
  ))
  suit_stack <- rast(required_file(
    file.path(topk_raster_dir, scenario, "ranked_suitability.tif")
  ))
  zone <- raster_points(zone_stack[[rank_value]], TRUE, max_cells)
  suitability <- raster_points(suit_stack[[rank_value]], FALSE, max_cells)
  setnames(zone, "value", "zoneID")
  setnames(suitability, "value", "dual_suitability")
  answer <- merge(zone, suitability, by = c("x", "y"), all = FALSE)
  answer[, `:=`(
    rank = rank_value,
    zone = factor(as.integer(zoneID), levels = zone_palette$zoneID)
  )]
  answer
}

run_figure("topk_rank_panels", {
  for (scenario in stack_scenarios) {
    zone_stack <- rast(required_file(
      file.path(topk_raster_dir, scenario, "ranked_zone.tif")
    ))
    suit_stack <- rast(required_file(
      file.path(topk_raster_dir, scenario, "ranked_suitability.tif")
    ))
    panels <- list()
    for (rank_value in k_values) {
      panels[[length(panels) + 1L]] <- zone_map_plot(
        zone_stack[[rank_value]], paste0("Rank ", rank_value, " zone")
      )
      panels[[length(panels) + 1L]] <- continuous_map_plot(
        suit_stack[[rank_value]], paste0("Rank ", rank_value, " dual suitability"),
        limits = c(0, 1)
      )
    }
    combined <- wrap_plots(panels, ncol = 2, guides = "collect") +
      plot_annotation(
        title = paste0(scenario_label(scenario), ": ranked zones and dual suitability"),
        caption = "Rows are ranks 1–5. Each zone map is paired with its continuous suitability values."
      ) & theme(legend.position = "right")
    save_gg(
      combined,
      file.path(figure_dir, "topk", paste0("rank_panel_", scenario, ".png")),
      11, 20
    )
  }
})

run_figure("topk_exploded_stack", {
  for (scenario in stack_scenarios) {
    all_ranks <- rbindlist(lapply(5:1, function(rank_value) {
      rank_pair_data(scenario, rank_value)
    }))
    x_span <- diff(range(all_ranks$x, na.rm = TRUE))
    y_span <- diff(range(all_ranks$y, na.rm = TRUE))
    all_ranks[, layer_height := 5L - rank]
    all_ranks[, `:=`(
      x_plot = x + layer_height * x_span * 0.045,
      y_plot = y + layer_height * y_span * 0.075
    )]

    p <- ggplot()
    for (rank_value in 5:1) {
      layer_data <- all_ranks[rank == rank_value]
      p <- p + geom_raster(
        data = layer_data,
        aes(x_plot, y_plot, fill = zone, alpha = dual_suitability)
      )
    }
    labels <- all_ranks[, .(
      x_plot = min(x_plot), y_plot = max(y_plot)
    ), by = rank]
    labels[, label := paste0("Rank ", rank)]
    p <- p +
      geom_label(
        data = labels, aes(x_plot, y_plot, label = label),
        inherit.aes = FALSE, size = 3, label.size = 0.15
      ) +
      scale_fill_manual(values = zone_colors, labels = zone_labels, drop = TRUE) +
      scale_alpha_continuous(
        limits = c(0, 1), range = c(0.12, 1),
        name = "Dual suitability"
      ) +
      coord_equal(expand = FALSE) +
      labs(
        title = paste0(scenario_label(scenario), ": exploded Top-5 dual-suitability stack"),
        subtitle = "Rank 5 is the bottom layer; ranks 4, 3, 2 and 1 are offset upward in that order.",
        caption = "Zone is encoded by hue and its dual-suitability value by opacity; the paired 5×2 panel provides the exact continuous scale.",
        x = NULL, y = NULL, fill = "Ecotype zone"
      ) + theme_void(base_size = 10) +
      theme(
        legend.position = "right",
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10)
      )
    save_gg(
      p,
      file.path(figure_dir, "topk", paste0("rank_exploded_", scenario, ".png")),
      13, 10
    )
  }
})


# 8. Audit =====================================================================

figure_audit <- rbindlist(figure_log, fill = TRUE)
fwrite(figure_audit, file.path(figure_dir, "figure_log.csv"))

failed <- figure_audit[status == "failed"]
if (nrow(failed)) {
  stop(
    nrow(failed), " figure family/families failed. See figures/figure_log.csv: ",
    paste(failed$figure, collapse = ", ")
  )
}

cat(
  "\nVISUALIZATION COMPLETE\n",
  "Figure families: ", nrow(figure_audit), "\n",
  "Output: ", figure_dir, "\n",
  sep = ""
)
