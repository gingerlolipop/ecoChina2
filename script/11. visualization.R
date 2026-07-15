# 11. visualization.R | rebuilt robust version
# ============================================================
# Generate manuscript tables and figures from existing results.
#
# Main analysis line:
#   One-hot binary RF models for climate and soil niches:
#     1. plain_rf
#     2. plain_mf_rf
#     3. optimized_rf
#     4. optimized_mf_rf
#
# Key outputs:
#   visualization/tables/
#   visualization/figures/
#   visualization/visualization_step_log.csv
#
# Notes:
#   - Each table/figure is an independent step, but the run is strict: any
#     incomplete or invalid step makes the script stop after writing the audit log.
#   - Missing files, partial method/scenario grids, ambiguous inputs, duplicated
#     records, and inconsistent metrics are treated as errors.
#   - Zone-related figures use color_palette_China.csv whenever possible.
#   - Figure 1 has two panels:
#       Figure 1a = mean performance across modeled zones.
#       Figure 1b = zone-level performance values.
#   - Figure 2a-d show predicted reference maps only; the original map is saved separately.
#   - Figure 6 now has Figure 6a (area change) and Figure 6b (zone-transition sankey).
#   - Figures 8 and 9 plot all projection-eligible tree species, 10 species per page
#     (2 rows x 5 columns), for all scenarios.
#   - Table 3 and Figure 7 use reference abundance from the code 6 population
#     lookup; abundance is never inferred from unique species-zone rows. The
#     code 6 `projected` flag means eligible for projection because the source
#     zone was modeled; it is not a future survival or persistence result.
#   - All reference-map metrics are recalculated directly from the original
#     raster and the assigned-zone rasters. Existing metrics CSV files are audit
#     inputs only and are never used as the source of manuscript results.
#   - Binary-only comparisons use the common valid mask of the four binary maps;
#     binary-versus-multiclass comparisons use the common valid mask of all five
#     maps. Testing-set and map-reconstruction metrics are never mixed.
#   - Critical inputs use explicit paths and structural checks. Ambiguous or
#     duplicated records stop the affected step instead of being averaged.
# ============================================================

library(terra)
library(data.table)
library(ggplot2)

rm(list = ls())
gc()

# 0. Settings ==================================================================

base_dir <- "H:/Jing/ecoChina2"

vis_dir <- file.path(base_dir, "visualization")
fig_dir <- file.path(vis_dir, "figures")
tab_dir <- file.path(vis_dir, "tables")

dir.create(vis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

reference_file <- file.path(base_dir, "raster/ecosys_ori.tif")
result_map_root <- file.path(base_dir, "result maps")
dual_root <- file.path(base_dir, "dual suit")

assessment_dir <- file.path(base_dir, "assessment")
binary_test_metrics_file <- file.path(assessment_dir, "rf_test_zone_metrics.csv")
binary_map_overall_file <- file.path(assessment_dir, "normal_map_overall_metrics.csv")
binary_map_zone_file <- file.path(assessment_dir, "normal_map_zone_metrics.csv")
binary_map_confusion_file <- file.path(assessment_dir, "normal_map_confusion_long.csv")

multiclass_assessment_dir <- file.path(assessment_dir, "multiclass_rf")
map_comparison_overall_file <- file.path(
  multiclass_assessment_dir,
  "multiclass_vs_overlay_map_overall_metrics.csv"
)
map_comparison_zone_file <- file.path(
  multiclass_assessment_dir,
  "multiclass_vs_overlay_map_zone_metrics.csv"
)
multiclass_reference_map_file <- file.path(
  result_map_root,
  "multiclass_rf",
  "assigned_zone_normal_multiclass_rf.tif"
)

population_lookup_file <- file.path(
  base_dir,
  "future tree niche",
  "tables",
  "population_projection_lookup.csv"
)
dual_population_lookup_file <- file.path(
  base_dir,
  "future tree niche dual suitability",
  "tables",
  "dual_population_projection_lookup_from_code6.csv"
)

preferred_method <- "optimized_rf"
tree_plot_method <- preferred_method
tree_plot_scenarios <- c(
  "normal",
  "2011-2040SSP245",
  "2041-2070SSP245",
  "2071-2100SSP245",
  "2011-2040SSP585",
  "2041-2070SSP585",
  "2071-2100SSP585"
)

# Dual suitability threshold used only for plotting population/species niches.
dual_plot_threshold <- 0.2
# Code 1.2 retains a species-zone population only when at least 10 occupied
# raster cells are present. Code 6 carries those retained abundances forward.
min_reference_population_abundance <- 10L

# Modeled zones. Zones 8 and 51 were not modeled, but remain part of the
# original reference map and reference population-abundance summaries.
model_zoneID <- c(1:7, 9:50, 52:55)
unmodeled_zoneID <- c(8L, 51L)
reference_zoneID <- sort(unique(c(model_zoneID, unmodeled_zoneID)))
novel_value <- 99

# Model folders used by assigned-zone maps.
map_method_order <- c(
  "plain_rf",
  "plain_mf",
  "optimized_rf",
  "optimized_mf"
)

map_method_labels <- c(
  plain_rf = "Plain RF",
  plain_mf = "Plain MF RF",
  optimized_rf = "Optimized RF",
  optimized_mf = "Optimized MF RF"
)

all_map_method_order <- c(map_method_order, "multiclass_rf")
all_map_method_labels <- c(
  map_method_labels,
  multiclass_rf = "Multiclass RF"
)

# Workflow labels used by binary climate/soil RF summaries.
workflow_order <- c(
  "plain_rf",
  "plain_mf_rf",
  "optimized_rf",
  "optimized_mf_rf"
)

workflow_labels <- c(
  plain_rf = "Plain RF",
  plain_mf_rf = "Plain MF RF",
  optimized_rf = "Optimized RF",
  optimized_mf_rf = "Optimized MF RF"
)

future_order <- c(
  "2011-2040SSP245",
  "2041-2070SSP245",
  "2071-2100SSP245",
  "2011-2040SSP585",
  "2041-2070SSP585",
  "2071-2100SSP585"
)

scenario_levels <- c("normal", future_order)

fig_dpi <- 320
display_max_cells <- 300000
species_page_ncol <- 5
species_page_nrow <- 2

# Never report COMPLETE when any manuscript output failed or was incomplete.
strict_mode <- TRUE
# Continue through independent figure/table steps so one run reveals every
# remaining problem. At the end, the script still stops without COMPLETE if
# any step failed. Figure 10a remains the first output.
stop_on_first_error <- FALSE
# Resume mode preserves completed expensive outputs, especially Figures 8 and 9.
# Each step still replaces its own outputs when that step is actually rerun.
resume_existing_outputs <- TRUE

# Do not globally delete completed outputs at startup. Failed or rerun steps clean
# their own managed files. Set TRUE only for a deliberate clean rebuild.
purge_previous_outputs <- FALSE

terraOptions(memfrac = 0.10)

# Non-zone colors. Zone-related figures use the vegetation-zone palette.
niche_cols <- c(
  climate = "#2E5E8C",
  soil = "#9A7B32"
)

# 1. General helpers ============================================================

cat0 <- function(...) {
  cat(..., "\n", sep = "")
}

norm_name <- function(x) {
  tolower(gsub("[^a-z0-9]+", "_", x))
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

pick_col <- function(dt, candidates) {
  if (is.null(dt) || ncol(dt) == 0) return(NA_character_)
  nms <- names(dt)
  nms_norm <- norm_name(nms)
  cand_norm <- norm_name(candidates)
  hit <- match(cand_norm, nms_norm)
  hit <- hit[!is.na(hit)]
  if (length(hit) == 0) return(NA_character_)
  nms[hit[1]]
}

resolve_exact_file <- function(candidates, label, required = TRUE) {
  candidates <- unique(as.character(candidates))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  existing <- candidates[file.exists(candidates)]
  
  if (length(existing) > 0) return(existing[1])
  
  if (required) {
    stop(
      "Missing ", label, ". Expected one of:\n",
      paste0("  - ", candidates, collapse = "\n")
    )
  }
  
  NA_character_
}

remove_existing_outputs <- function(files) {
  files <- unique(as.character(files))
  files <- files[!is.na(files) & file.exists(files)]
  if (length(files) > 0) {
    unlink(files, force = TRUE)
    cat0("[REMOVED OLD OUTPUT] ", paste(basename(files), collapse = ", "))
  }
  invisible(files)
}

require_columns <- function(dt, required, label) {
  missing_cols <- setdiff(required, names(dt))
  if (length(missing_cols) > 0) {
    stop(
      label,
      " is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  invisible(TRUE)
}

assert_unique_keys <- function(dt, keys, label) {
  require_columns(dt, keys, label)
  dup <- dt[, .N, by = keys][N > 1]
  if (nrow(dup) > 0) {
    stop(
      label,
      " contains duplicated key combinations for: ",
      paste(keys, collapse = " + "),
      ". The first duplicated keys are: ",
      paste(capture.output(print(head(dup, 10))), collapse = " ")
    )
  }
  invisible(TRUE)
}

require_methods <- function(dt, required_methods, method_col = "method_key", label) {
  require_columns(dt, method_col, label)
  missing_methods <- setdiff(required_methods, unique(as.character(dt[[method_col]])))
  if (length(missing_methods) > 0) {
    stop(
      label,
      " is missing required methods: ",
      paste(missing_methods, collapse = ", ")
    )
  }
  invisible(TRUE)
}

require_complete_zone_grid <- function(
    dt,
    required_methods,
    required_zones,
    method_col = "method_key",
    zone_col = "zoneID",
    label) {
  require_columns(dt, c(method_col, zone_col), label)
  
  expected <- CJ(
    method_value = as.character(required_methods),
    zone_value = as.integer(required_zones),
    unique = TRUE
  )
  observed <- unique(data.table(
    method_value = as.character(dt[[method_col]]),
    zone_value = as.integer(dt[[zone_col]])
  ))
  
  missing <- fsetdiff(expected, observed)
  extra <- fsetdiff(observed, expected)
  
  if (nrow(missing) > 0 || nrow(extra) > 0) {
    stop(
      label,
      " does not contain the exact method-by-zone grid. Missing: ",
      paste(capture.output(print(head(missing, 12))), collapse = " "),
      "; extra: ",
      paste(capture.output(print(head(extra, 12))), collapse = " ")
    )
  }
  
  invisible(TRUE)
}



recompute_classification_metrics <- function(dt, label) {
  dt <- copy(dt)
  require_columns(dt, c("TP", "TN", "FP", "FN"), label)
  
  count_cols <- c("TP", "TN", "FP", "FN")
  for (cc in count_cols) {
    set(dt, j = cc, value = as.numeric(dt[[cc]]))
  }
  
  # Extract ordinary numeric vectors once. Do not refer to TP/TN/FP/FN as
  # unquoted symbols inside data.table expressions; that caused the previous
  # "object 'TP' not found" error on some data.table versions/environments.
  tp <- dt[["TP"]]
  tn <- dt[["TN"]]
  fp <- dt[["FP"]]
  fn <- dt[["FN"]]
  
  bad_index <- (
    !complete.cases(as.data.frame(dt[, ..count_cols])) |
      !is.finite(tp) |
      !is.finite(tn) |
      !is.finite(fp) |
      !is.finite(fn) |
      tp < 0 |
      tn < 0 |
      fp < 0 |
      fn < 0
  )
  
  if (any(bad_index)) {
    bad <- dt[which(bad_index)]
    stop(
      label,
      " contains missing, non-finite, or negative confusion counts. First rows: ",
      paste(capture.output(print(head(bad, 10))), collapse = " ")
    )
  }
  
  total_n <- tp + tn + fp + fn
  accuracy_value <- ifelse(total_n > 0, (tp + tn) / total_n, NA_real_)
  recall_value <- ifelse(tp + fn > 0, tp / (tp + fn), NA_real_)
  specificity_value <- ifelse(tn + fp > 0, tn / (tn + fp), NA_real_)
  precision_value <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
  f1_value <- ifelse(
    2 * tp + fp + fn > 0,
    2 * tp / (2 * tp + fp + fn),
    NA_real_
  )
  balanced_accuracy_value <- (recall_value + specificity_value) / 2
  tss_value <- recall_value + specificity_value - 1
  
  set(dt, j = "accuracy", value = as.numeric(accuracy_value))
  set(dt, j = "recall", value = as.numeric(recall_value))
  set(dt, j = "specificity", value = as.numeric(specificity_value))
  set(dt, j = "precision", value = as.numeric(precision_value))
  set(dt, j = "f1", value = as.numeric(f1_value))
  set(
    dt,
    j = "balanced_accuracy",
    value = as.numeric(balanced_accuracy_value)
  )
  set(dt, j = "tss", value = as.numeric(tss_value))
  
  dt[]
}


validate_confusion_derived_counts <- function(dt, label, tolerance = 1e-10) {
  require_columns(dt, c("TP", "FP", "FN"), label)
  
  tp <- as.numeric(dt[["TP"]])
  fp <- as.numeric(dt[["FP"]])
  fn <- as.numeric(dt[["FN"]])
  
  if ("original_pixels" %in% names(dt)) {
    original_pixels_value <- as.numeric(dt[["original_pixels"]])
    bad_index <- (
      !is.na(original_pixels_value) &
        abs(original_pixels_value - (tp + fn)) > tolerance
    )
    
    if (any(bad_index)) {
      bad <- dt[which(bad_index)]
      stop(
        label,
        " has original_pixels inconsistent with TP + FN. First rows: ",
        paste(capture.output(print(head(bad, 10))), collapse = " ")
      )
    }
  }
  
  if ("predicted_pixels" %in% names(dt)) {
    predicted_pixels_value <- as.numeric(dt[["predicted_pixels"]])
    bad_index <- (
      !is.na(predicted_pixels_value) &
        abs(predicted_pixels_value - (tp + fp)) > tolerance
    )
    
    if (any(bad_index)) {
      bad <- dt[which(bad_index)]
      stop(
        label,
        " has predicted_pixels inconsistent with TP + FP. First rows: ",
        paste(capture.output(print(head(bad, 10))), collapse = " ")
      )
    }
  }
  
  if ("predicted_to_original" %in% names(dt)) {
    observed_ratio <- as.numeric(dt[["predicted_to_original"]])
    expected_ratio <- ifelse(
      tp + fn > 0,
      (tp + fp) / (tp + fn),
      NA_real_
    )
    
    bad_index <- (
      xor(is.na(observed_ratio), is.na(expected_ratio)) |
        (
          !is.na(observed_ratio) &
            !is.na(expected_ratio) &
            abs(observed_ratio - expected_ratio) > tolerance
        )
    )
    
    if (any(bad_index)) {
      bad <- dt[which(bad_index)]
      stop(
        label,
        " has predicted_to_original inconsistent with confusion counts. ",
        "First rows: ",
        paste(capture.output(print(head(bad, 10))), collapse = " ")
      )
    }
  }
  
  invisible(TRUE)
}

validate_metric_ranges <- function(dt, unit_cols = character(), tss_cols = character(), label) {
  for (cc in intersect(unit_cols, names(dt))) {
    bad <- dt[!is.na(get(cc)) & (!is.finite(get(cc)) | get(cc) < -1e-10 | get(cc) > 1 + 1e-10)]
    if (nrow(bad) > 0) {
      stop(
        label, " has values outside [0, 1] in ", cc, ". First rows: ",
        paste(capture.output(print(head(bad, 10))), collapse = " ")
      )
    }
  }
  for (cc in intersect(tss_cols, names(dt))) {
    bad <- dt[!is.na(get(cc)) & (!is.finite(get(cc)) | get(cc) < -1 - 1e-10 | get(cc) > 1 + 1e-10)]
    if (nrow(bad) > 0) {
      stop(
        label, " has values outside [-1, 1] in ", cc, ". First rows: ",
        paste(capture.output(print(head(bad, 10))), collapse = " ")
      )
    }
  }
  invisible(TRUE)
}


validate_overall_against_zone_counts <- function(
    overall,
    zone,
    label,
    tolerance = 1e-10) {
  
  require_columns(
    overall,
    c("method_key", "coverage", "exact_zone_accuracy"),
    label
  )
  require_columns(
    zone,
    c("method_key", "TP", "TN", "FP", "FN"),
    label
  )
  
  method_values <- unique(as.character(zone[["method_key"]]))
  
  derived_list <- lapply(method_values, function(method_value) {
    current <- zone[
      as.character(zone[["method_key"]]) == method_value
    ]
    
    tp <- as.numeric(current[["TP"]])
    tn <- as.numeric(current[["TN"]])
    fp <- as.numeric(current[["FP"]])
    fn <- as.numeric(current[["FN"]])
    
    totals <- unique(tp + tn + fp + fn)
    
    if (
      length(totals) != 1L ||
      !is.finite(totals[1]) ||
      totals[1] <= 0
    ) {
      stop(
        label,
        " has inconsistent per-zone comparison totals for method ",
        method_value
      )
    }
    
    data.table(
      method_key = method_value,
      compared_from_zone_counts = as.numeric(totals[1]),
      exact_from_zone_counts = as.numeric(sum(tp) / totals[1])
    )
  })
  
  derived <- rbindlist(
    derived_list,
    use.names = TRUE,
    fill = FALSE
  )
  
  check <- merge(
    overall,
    derived,
    by = "method_key",
    all.x = TRUE,
    sort = FALSE
  )
  
  check[, exact_difference := abs(
    as.numeric(exact_zone_accuracy) -
      as.numeric(exact_from_zone_counts)
  )]
  
  if ("compared_pixels" %in% names(check)) {
    check[, compared_difference := abs(
      as.numeric(compared_pixels) -
        as.numeric(compared_from_zone_counts)
    )]
  } else {
    check[, compared_difference := 0]
  }
  
  if (
    all(
      c(
        "valid_original_pixels",
        "compared_pixels"
      ) %in% names(check)
    )
  ) {
    check[, coverage_from_counts := (
      as.numeric(compared_pixels) /
        as.numeric(valid_original_pixels)
    )]
    check[, coverage_difference := abs(
      as.numeric(coverage) -
        as.numeric(coverage_from_counts)
    )]
  } else {
    check[, coverage_difference := 0]
  }
  
  bad <- check[
    !is.finite(exact_difference) |
      exact_difference > tolerance |
      !is.finite(compared_difference) |
      compared_difference > tolerance |
      !is.finite(coverage_difference) |
      coverage_difference > tolerance
  ]
  
  if (nrow(bad) > 0) {
    stop(
      label,
      " overall metrics are inconsistent with zone confusion counts. ",
      "First rows: ",
      paste(capture.output(print(head(bad, 10))), collapse = " ")
    )
  }
  
  invisible(TRUE)
}

relative_path <- function(x) {
  base_norm <- normalizePath(base_dir, winslash = "/", mustWork = FALSE)
  x_norm <- normalizePath(x, winslash = "/", mustWork = FALSE)
  gsub(paste0("^", base_norm, "/?"), "", x_norm)
}

theme_ms <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = "grey70"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      axis.text = element_text(colour = "grey20"),
      legend.title = element_text(face = "bold")
    )
}

theme_fig1 <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.35),
      strip.background = element_rect(fill = "grey95", colour = "grey78", linewidth = 0.5),
      strip.text = element_text(face = "bold", colour = "grey15"),
      axis.text = element_text(colour = "grey20"),
      axis.title = element_text(colour = "grey10", face = "bold"),
      plot.title = element_text(face = "bold", size = rel(1.18), colour = "grey10"),
      plot.subtitle = element_text(colour = "grey20"),
      plot.caption = element_text(colour = "grey35", size = rel(0.88)),
      legend.title = element_text(face = "bold"),
      panel.spacing.x = grid::unit(1.05, "lines"),
      panel.spacing.y = grid::unit(0.42, "lines"),
      plot.margin = margin(8, 10, 8, 8)
    )
}

save_gg <- function(p, file, width = 7.5, height = 5.5) {
  ggsave(
    filename = file,
    plot = p,
    width = width,
    height = height,
    dpi = fig_dpi,
    bg = "white"
  )
  cat0("[SAVED] ", file)
}

step_log <- list()

run_step <- function(step_name, code) {
  cat0("\n============================================================")
  cat0(step_name)
  cat0("============================================================")
  t0 <- Sys.time()
  tryCatch(
    {
      force(code)
      elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
      step_log[[length(step_log) + 1L]] <<- data.table(
        step = step_name,
        status = "done",
        message = NA_character_,
        elapsed_sec = elapsed
      )
      cat0("[DONE] ", step_name, " | ", elapsed, " sec")
    },
    error = function(e) {
      elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
      step_log[[length(step_log) + 1L]] <<- data.table(
        step = step_name,
        status = "skipped_or_error",
        message = conditionMessage(e),
        elapsed_sec = elapsed
      )
      cat0("[ERROR] ", step_name)
      cat0("  ", conditionMessage(e))
      
      # Persist the failure immediately. In strict mode, stop at the first
      # invalid manuscript output rather than continuing with dependent steps.
      partial_log <- rbindlist(step_log, fill = TRUE)
      partial_log_file <- file.path(vis_dir, "visualization_step_log.csv")
      fwrite(partial_log, partial_log_file)
      if (strict_mode && stop_on_first_error) {
        stop(
          "Visualization stopped at: ", step_name,
          ". See ", partial_log_file,
          " | ", conditionMessage(e),
          call. = FALSE
        )
      } else {
        cat0("[CONTINUE] Later independent steps will still be attempted.")
      }
    }
  )
  invisible(NULL)
}

extract_zone_id <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(as.integer(x))
  miss <- is.na(out)
  if (any(miss)) {
    xx <- x[miss]
    xx <- sub(".*?(\\d+).*", "\\1", xx)
    out[miss] <- suppressWarnings(as.integer(xx))
  }
  out
}

# 2. Zone palette ==============================================================

default_zone_colors <- function(vals) {
  vals <- sort(unique(as.integer(vals)))
  vals <- vals[!is.na(vals)]
  cols <- grDevices::hcl.colors(n = length(vals), palette = "Dark 3")
  names(cols) <- as.character(vals)
  if ("8" %in% names(cols)) cols["8"] <- "#BDBDBD"
  if ("99" %in% names(cols)) cols["99"] <- "#333333"
  cols
}

load_zone_lookup <- function() {
  vals <- sort(unique(c(reference_zoneID, novel_value)))
  cols <- default_zone_colors(vals)
  lookup <- data.table(
    zoneID = vals,
    zone_label = ifelse(vals == novel_value, "Novel ecotype", paste0("Zone ", vals)),
    color = unname(cols[as.character(vals)])
  )
  
  # Use the project palette only. Do not select a similarly named CSV elsewhere.
  pal_file <- file.path(base_dir, "color_palette_China.csv")
  
  if (file.exists(pal_file)) {
    pal <- fread(pal_file)
    z_col <- pick_col(pal, c("zoneID", "zone_id", "zoneid", "id", "value"))
    c_col <- pick_col(pal, c("COLOR", "color", "colour", "hex", "hex_color"))
    n_col <- pick_col(
      pal,
      c("zone_name", "zone_label", "zone", "name", "vegetation", "ecosystem", "type")
    )
    
    if (is.na(z_col) || is.na(c_col)) {
      stop("Cannot identify zoneID and COLOR columns in: ", pal_file)
    }
    
    tmp <- data.table(zoneID = as.integer(pal[[z_col]]))
    
    if (!is.na(c_col)) {
      tmp[, color_new := as.character(pal[[c_col]])]
      tmp[!grepl("^#", color_new) & !is.na(color_new), color_new := paste0("#", color_new)]
    }
    
    if (!is.na(n_col)) {
      tmp[, label_new := as.character(pal[[n_col]])]
    }
    
    assert_unique_keys(tmp[!is.na(zoneID)], "zoneID", "color_palette_China.csv")
    missing_palette_zones <- setdiff(reference_zoneID, tmp$zoneID)
    if (length(missing_palette_zones) > 0) {
      stop(
        "color_palette_China.csv is missing reference zones: ",
        paste(missing_palette_zones, collapse = ", ")
      )
    }
    bad_colors <- tmp[
      zoneID %in% reference_zoneID &
        (is.na(color_new) | !grepl("^#[0-9A-Fa-f]{6}$", color_new))
    ]
    if (nrow(bad_colors) > 0) {
      stop(
        "Invalid hexadecimal colors in color_palette_China.csv: ",
        paste(capture.output(print(head(bad_colors, 12))), collapse = " ")
      )
    }
    
    lookup <- merge(lookup, tmp, by = "zoneID", all.x = TRUE, sort = FALSE)
    
    if ("color_new" %in% names(lookup)) {
      lookup[!is.na(color_new) & nzchar(color_new), color := color_new]
      lookup[, color_new := NULL]
    }
    
    if ("label_new" %in% names(lookup)) {
      lookup[!is.na(label_new) & nzchar(label_new), zone_label := label_new]
      lookup[, label_new := NULL]
    }
  } else {
    stop("Missing required project vegetation-zone palette: ", pal_file)
  }
  
  lookup[zoneID == novel_value, `:=`(
    zone_label = "Novel ecotype",
    color = "#333333"
  )]
  
  lookup[]
}

zone_lookup <- load_zone_lookup()

zone_color_vector <- function(vals = c(model_zoneID, novel_value)) {
  vals <- sort(unique(as.integer(vals)))
  vals <- vals[!is.na(vals)]
  
  lk <- merge(
    data.table(zoneID = vals),
    zone_lookup[, .(zoneID, color)],
    by = "zoneID",
    all.x = TRUE
  )
  
  miss <- is.na(lk$color)
  if (any(miss)) {
    extra_cols <- default_zone_colors(lk$zoneID[miss])
    lk[miss, color := unname(extra_cols[as.character(zoneID)])]
  }
  
  cols <- lk$color
  names(cols) <- as.character(lk$zoneID)
  cols
}

zone_label_vector <- function(vals = c(model_zoneID, novel_value)) {
  vals <- sort(unique(as.integer(vals)))
  vals <- vals[!is.na(vals)]
  
  lk <- merge(
    data.table(zoneID = vals),
    zone_lookup[, .(zoneID, zone_label)],
    by = "zoneID",
    all.x = TRUE
  )
  
  lk[is.na(zone_label), zone_label := paste0("Zone ", zoneID)]
  labs <- lk$zone_label
  names(labs) <- as.character(lk$zoneID)
  labs
}

# 3. Binary RF helpers ==========================================================

infer_workflow <- function(x) {
  x <- tolower(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("optimized.*mf|optimized_mf", x)] <- "optimized_mf_rf"
  out[grepl("optimized.*rf|optimized_rf", x) & is.na(out)] <- "optimized_rf"
  out[grepl("plain.*mf|plain_mf", x) & is.na(out)] <- "plain_mf_rf"
  out[grepl("plain.*rf|plain_rf", x) & is.na(out)] <- "plain_rf"
  out
}

normalise_map_method <- function(x) {
  x <- infer_workflow(x)
  x[x == "plain_mf_rf"] <- "plain_mf"
  x[x == "optimized_mf_rf"] <- "optimized_mf"
  x
}

canonical_map_method <- function(x) {
  raw <- as.character(x)
  out <- normalise_map_method(raw)
  out[grepl("multiclass", raw, ignore.case = TRUE)] <- "multiclass_rf"
  out
}

metric_from_col <- function(x) {
  x0 <- norm_name(x)
  
  if (grepl("oob.*acc|acc.*oob|oob.*accuracy|accuracy.*oob", x0)) return("OOB accuracy")
  if (grepl("oob.*auc|auc.*oob", x0)) return("OOB AUC")
  if (grepl("oob.*f1|f1.*oob", x0)) return("OOB F1")
  if (grepl("oob.*precision|precision.*oob", x0)) return("OOB precision")
  if (grepl("oob.*recall|recall.*oob|oob.*sensitivity|sensitivity.*oob", x0)) return("OOB recall")
  if (grepl("oob.*specificity|specificity.*oob", x0)) return("OOB specificity")
  
  if (grepl("train.*acc|acc.*train|train.*accuracy|accuracy.*train", x0)) return("Train accuracy")
  if (grepl("train.*auc|auc.*train", x0)) return("Train AUC")
  if (grepl("train.*f1|f1.*train", x0)) return("Train F1")
  if (grepl("train.*precision|precision.*train", x0)) return("Train precision")
  if (grepl("train.*recall|recall.*train|train.*sensitivity|sensitivity.*train", x0)) return("Train recall")
  if (grepl("train.*specificity|specificity.*train", x0)) return("Train specificity")
  
  NA_character_
}

find_binary_summary_file <- function(niche_type) {
  if (niche_type == "climate") {
    candidates <- file.path(
      base_dir, "accuracy_climate", "climate_rf_accuracy_summary.csv"
    )
  } else if (niche_type == "soil") {
    candidates <- file.path(
      base_dir, "accuracy_soil", "soil_rf_accuracy_summary.csv"
    )
  } else {
    stop("Unsupported niche type: ", niche_type)
  }
  
  resolve_exact_file(
    candidates,
    paste0(niche_type, " binary RF accuracy summary"),
    required = TRUE
  )
}

read_train_oob_metrics <- function(niche_type, file) {
  if (is.na(file) || !file.exists(file)) {
    stop("Missing binary RF summary file for ", niche_type)
  }
  
  dt <- fread(file)
  
  z_col <- pick_col(dt, c("zoneID", "zone_id", "zone", "ecotype", "class", "veg_zone"))
  if (is.na(z_col)) stop("Cannot find zone column in: ", file)
  
  workflow_col <- pick_col(dt, c("workflow", "method", "model", "model_version", "version", "model_type", "rf_type"))
  
  dt[, zoneID_tmp := extract_zone_id(get(z_col))]
  dt <- dt[zoneID_tmp %in% model_zoneID]
  if (nrow(dt) == 0) stop("No modeled zones found in: ", file)
  
  if (is.na(workflow_col)) {
    stop(
      "Cannot identify the workflow/model column in: ",
      file,
      ". Workflow identity is never inferred from row order."
    )
  }
  
  dt[, workflow := infer_workflow(get(workflow_col))]
  dt[, workflow_source := workflow_col]
  dt <- dt[workflow %in% workflow_order]
  
  if (nrow(dt) == 0) {
    stop("No recognized binary workflows found in: ", file)
  }
  
  assert_unique_keys(
    dt,
    c("zoneID_tmp", "workflow"),
    paste0("Binary RF summary: ", relative_path(file))
  )
  require_complete_zone_grid(
    dt,
    workflow_order,
    model_zoneID,
    method_col = "workflow",
    zone_col = "zoneID_tmp",
    label = paste0("Binary RF summary: ", relative_path(file))
  )
  
  numeric_cols <- names(dt)[sapply(dt, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, "zoneID_tmp")
  
  out <- list()
  
  for (mc in numeric_cols) {
    metric <- metric_from_col(mc)
    if (is.na(metric)) next
    
    tmp <- data.table(
      niche_type = niche_type,
      zoneID = dt$zoneID_tmp,
      workflow = dt$workflow,
      metric = metric,
      value = as.numeric(dt[[mc]]),
      metric_raw = mc,
      workflow_source = dt$workflow_source,
      source_file = relative_path(file)
    )
    
    out[[length(out) + 1L]] <- tmp
  }
  
  if (length(out) == 0) stop("No usable OOB/Train metric columns found in: ", file)
  
  out <- rbindlist(out, fill = TRUE)
  out <- out[workflow %in% workflow_order & !is.na(value) & is.finite(value)]
  
  if (nrow(out) == 0) stop("No valid OOB/Train workflow metrics extracted from: ", file)
  validate_metric_ranges(
    out,
    unit_cols = "value",
    label = paste0("Binary RF OOB/Train metrics: ", relative_path(file))
  )
  
  out[]
}

find_test_metrics_file <- function() {
  resolve_exact_file(
    binary_test_metrics_file,
    "binary independent-test zone metrics",
    required = TRUE
  )
}

read_test_metrics <- function(file) {
  if (is.na(file) || !file.exists(file)) {
    stop("Missing binary independent-test metrics: assessment/rf_test_zone_metrics.csv")
  }
  
  dt <- fread(file)
  
  niche_col <- pick_col(dt, c("niche", "niche_type", "type"))
  workflow_col <- pick_col(dt, c("method", "workflow", "model", "model_type"))
  zone_col <- pick_col(dt, c("zoneID", "zone_id", "zone", "class", "veg_zone"))
  
  if (is.na(niche_col) || is.na(workflow_col) || is.na(zone_col)) {
    stop("Cannot identify niche/method/zone columns in: ", file)
  }
  
  dt[, niche_type := tolower(as.character(get(niche_col)))]
  dt[, workflow := infer_workflow(get(workflow_col))]
  dt[, zoneID := extract_zone_id(get(zone_col))]
  
  dt <- dt[
    niche_type %in% c("climate", "soil") &
      workflow %in% workflow_order &
      zoneID %in% model_zoneID
  ]
  
  if (nrow(dt) == 0) stop("No valid rows found in: ", file)
  
  assert_unique_keys(
    dt,
    c("niche_type", "workflow", "zoneID"),
    paste0("Binary independent-test metrics: ", relative_path(file))
  )
  
  for (nt in c("climate", "soil")) {
    require_complete_zone_grid(
      dt[niche_type == nt],
      workflow_order,
      model_zoneID,
      method_col = "workflow",
      zone_col = "zoneID",
      label = paste0("Binary independent-test metrics (", nt, ")")
    )
  }
  
  dt <- recompute_classification_metrics(
    dt,
    paste0("Binary independent-test metrics: ", relative_path(file))
  )
  validate_metric_ranges(
    dt,
    unit_cols = c(
      "accuracy", "balanced_accuracy", "recall", "specificity",
      "precision", "f1", "auc"
    ),
    tss_cols = "tss",
    label = "Binary independent-test metrics"
  )
  
  metric_cols <- c(
    accuracy = "Test accuracy",
    balanced_accuracy = "Test balanced accuracy",
    f1 = "Test F1",
    precision = "Test precision",
    recall = "Test recall",
    specificity = "Test specificity",
    tss = "Test TSS"
  )
  
  # AUC cannot be reconstructed from the confusion matrix; retain the exact
  # value produced by code 5 when present.
  if ("auc" %in% names(dt)) metric_cols <- c(metric_cols, auc = "Test AUC")
  
  out <- rbindlist(
    lapply(names(metric_cols), function(mc) {
      data.table(
        niche_type = dt$niche_type,
        zoneID = dt$zoneID,
        workflow = dt$workflow,
        metric = unname(metric_cols[mc]),
        value = as.numeric(dt[[mc]]),
        metric_raw = mc,
        workflow_source = "assessment/rf_test_zone_metrics.csv",
        source_file = relative_path(file)
      )
    }),
    fill = TRUE
  )
  
  if (any(!is.finite(out$value))) {
    bad <- out[!is.finite(value)]
    stop(
      "Non-finite binary test metrics remain after count-based recalculation. First rows: ",
      paste(capture.output(print(head(bad, 10))), collapse = " ")
    )
  }
  
  out[]
}

make_binary_rf_table1 <- function() {
  climate_file <- find_binary_summary_file("climate")
  soil_file <- find_binary_summary_file("soil")
  test_file <- find_test_metrics_file()
  
  cat0("[TABLE 1 SOURCE] climate: ", climate_file)
  cat0("[TABLE 1 SOURCE] soil: ", soil_file)
  cat0("[TABLE 1 SOURCE] test: ", test_file)
  source_table <- data.table(
    niche_type = c("climate", "soil", "test"),
    source_file = c(climate_file, soil_file, test_file)
  )
  fwrite(source_table, file.path(tab_dir, "Table1_source_files_used.csv"))
  
  long <- rbindlist(
    list(
      read_train_oob_metrics("climate", climate_file),
      read_train_oob_metrics("soil", soil_file),
      read_test_metrics(test_file)
    ),
    fill = TRUE
  )
  
  assert_unique_keys(
    long,
    c("niche_type", "workflow", "zoneID", "metric"),
    "Combined binary RF performance table"
  )
  metric_grid_audit <- long[
    ,
    .(
      n_zones = uniqueN(zoneID),
      missing_zones = paste(
        setdiff(model_zoneID, sort(unique(as.integer(zoneID)))),
        collapse = ","
      )
    ),
    by = .(niche_type, workflow, metric)
  ]
  metric_grid_audit[, complete_grid := n_zones == length(model_zoneID)]
  fwrite(
    metric_grid_audit,
    file.path(tab_dir, "Table1_metric_grid_audit.csv")
  )
  
  incomplete_metric_groups <- metric_grid_audit[complete_grid == FALSE]
  if (nrow(incomplete_metric_groups) > 0) {
    cat0(
      "[TABLE 1 EXCLUDED INCOMPLETE METRICS] ",
      paste(
        paste0(
          incomplete_metric_groups$niche_type, " | ",
          incomplete_metric_groups$workflow, " | ",
          incomplete_metric_groups$metric, " (",
          incomplete_metric_groups$n_zones, "/",
          length(model_zoneID), " zones)"
        ),
        collapse = "; "
      )
    )
    
    # Anti-join removes only metric groups that cannot be compared across the
    # full modeled-zone grid. This is appropriate for unavailable multi-Forest
    # OOB statistics and prevents a single repaired Zone from creating a false
    # manuscript-wide OOB comparison.
    long <- long[
      !incomplete_metric_groups,
      on = .(niche_type, workflow, metric)
    ]
  }
  
  if (nrow(long) == 0) {
    stop("No complete binary RF metric groups remain after grid auditing.")
  }
  
  # No duplicate rows are averaged. Every retained value has the complete
  # modeled-zone grid and preserves its original source.
  long[, workflow := factor(workflow, levels = workflow_order)]
  
  summary_long <- long[
    ,
    {
      n_z <- uniqueN(zoneID)
      sd_v <- if (n_z > 1) sd(value, na.rm = TRUE) else 0
      .(
        n_zones = n_z,
        mean = mean(value, na.rm = TRUE),
        sd = sd_v,
        se = sd_v / sqrt(n_z),
        min = min(value, na.rm = TRUE),
        max = max(value, na.rm = TRUE)
      )
    },
    by = .(niche_type, workflow, metric)
  ]
  
  summary_long[, `:=`(
    ci95_low = mean - 1.96 * se,
    ci95_high = mean + 1.96 * se
  )]
  
  summary_long[!grepl("TSS", metric), `:=`(
    ci95_low = pmax(ci95_low, 0),
    ci95_high = pmin(ci95_high, 1)
  )]
  
  summary_long[grepl("TSS", metric), `:=`(
    ci95_low = pmax(ci95_low, -1),
    ci95_high = pmin(ci95_high, 1)
  )]
  summary_long[, interval_type :=
                 "Descriptive across-zone mean +/- 1.96 SE; not a model-based confidence interval"]
  
  summary_long[, workflow_label := workflow_labels[as.character(workflow)]]
  
  setorder(summary_long, niche_type, workflow, metric)
  setorder(long, niche_type, workflow, zoneID, metric)
  
  list(long = long, summary_long = summary_long)
}

# 4. Raster plotting helpers ===================================================

thin_raster_for_plot <- function(x, categorical = TRUE, max_cells = display_max_cells) {
  x <- x[[1]]
  if (ncell(x) <= max_cells) return(x)
  
  fact <- ceiling(sqrt(ncell(x) / max_cells))
  fun <- if (categorical) "modal" else "mean"
  
  aggregate(x, fact = fact, fun = fun, na.rm = TRUE)
}

raster_to_plot_dt <- function(
    file,
    categorical = TRUE,
    population_mode = FALSE,
    allowed_values = NULL) {
  if (is.na(file) || !file.exists(file)) stop("Missing raster file: ", file)
  
  x <- rast(file)[[1]]
  x <- thin_raster_for_plot(x, categorical = categorical)
  
  dt <- as.data.table(as.data.frame(x, xy = TRUE, na.rm = TRUE))
  if (nrow(dt) == 0) return(data.table())
  
  v_col <- setdiff(names(dt), c("x", "y"))[1]
  setnames(dt, v_col, "value")
  
  if (categorical) {
    dt[, value := as.integer(round(value))]
    
    if (!population_mode) {
      if (is.null(allowed_values)) {
        allowed_values <- c(model_zoneID, novel_value)
      }
      allowed_values <- unique(as.integer(allowed_values))
      dt <- dt[value %in% allowed_values]
    }
  }
  
  dt[!is.na(value) & is.finite(value)]
}

plot_zone_panel_gg <- function(
    files,
    titles,
    outfile,
    categorical = TRUE,
    population_mode = FALSE,
    show_legend = FALSE,
    ncol = 2,
    width = 10,
    height = 6,
    draw_as_points = FALSE,
    point_alpha = 1,
    point_size = 0.08,
    allowed_values = NULL,
    plot_title = NULL) {
  
  plot_list <- list()
  
  if (length(files) != length(titles)) {
    stop("files and titles must have the same length.")
  }
  if (length(files) == 0) stop("No raster panels were supplied.")
  
  for (i in seq_along(files)) {
    cat0("  Preparing raster panel: ", titles[i])
    
    dt <- raster_to_plot_dt(
      files[[i]],
      categorical = categorical,
      population_mode = population_mode,
      allowed_values = allowed_values
    )
    
    if (nrow(dt) == 0) {
      stop("Raster panel has no valid cells after filtering: ", titles[i])
    }
    dt[, panel := titles[i]]
    plot_list[[length(plot_list) + 1L]] <- dt
  }
  
  plot_dt <- rbindlist(plot_list, fill = TRUE)
  plot_dt[, panel := factor(panel, levels = titles[titles %in% unique(panel)])]
  
  if (categorical) {
    vals <- sort(unique(plot_dt$value))
    vals <- vals[!is.na(vals)]
    cols <- zone_color_vector(vals)
    if (novel_value %in% vals) cols[as.character(novel_value)] <- "#333333"
    plot_dt[, value_chr := factor(as.character(value), levels = as.character(vals))]
    
    if (draw_as_points) {
      p <- ggplot(plot_dt, aes(x = x, y = y, colour = value_chr)) +
        geom_point(size = point_size, alpha = point_alpha) +
        scale_colour_manual(values = cols, drop = FALSE) +
        labs(colour = "Zone")
    } else {
      p <- ggplot(plot_dt, aes(x = x, y = y, fill = value_chr)) +
        geom_raster(alpha = point_alpha) +
        scale_fill_manual(values = cols, drop = FALSE) +
        labs(fill = "Zone")
    }
    
  } else {
    if (draw_as_points) {
      p <- ggplot(plot_dt, aes(x = x, y = y, colour = value)) +
        geom_point(size = point_size, alpha = point_alpha) +
        scale_colour_gradientn(
          colours = grDevices::hcl.colors(60, palette = "YlGnBu"),
          na.value = NA
        ) +
        labs(colour = "Suitability")
    } else {
      p <- ggplot(plot_dt, aes(x = x, y = y, fill = value)) +
        geom_raster(alpha = point_alpha) +
        scale_fill_gradientn(
          colours = grDevices::hcl.colors(60, palette = "YlGnBu"),
          na.value = NA
        ) +
        labs(fill = "Suitability")
    }
  }
  
  p <- p +
    facet_wrap(~ panel, ncol = ncol) +
    coord_equal(expand = FALSE) +
    labs(title = plot_title) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 13),
      strip.text = element_text(face = "bold", size = 9.5),
      legend.position = ifelse(show_legend, "right", "none"),
      legend.text = element_text(size = 7),
      legend.key.size = grid::unit(0.35, "cm")
    )
  
  ggsave(outfile, p, width = width, height = height, dpi = fig_dpi, bg = "white")
  cat0("[SAVED] ", outfile)
}

# 5. Scenario, assigned-map, and dual-niche helpers ============================

scenario_label <- function(s) {
  ifelse(
    s == "normal",
    "Predicted reference",
    paste0(sub("SSP[0-9]+$", "", s), "\n", sub("^.*(SSP[0-9]+)$", "\\1", s))
  )
}

scenario_period <- function(s) {
  ifelse(s == "normal", "Predicted reference", sub("SSP[0-9]+$", "", s))
}

scenario_ssp <- function(s) {
  ifelse(s == "normal", "Predicted reference", sub("^.*(SSP[0-9]+)$", "\\1", s))
}

assigned_map_cache <- new.env(parent = emptyenv())

validate_assigned_map_file <- function(file, method, scenario) {
  if (is.na(file) || !file.exists(file)) return(NA_character_)
  if (!file.exists(reference_file)) stop("Missing reference raster: ", reference_file)
  
  x <- rast(file)[[1]]
  ref <- rast(reference_file)[[1]]
  if (!isTRUE(compareGeom(x, ref, stopOnError = FALSE))) {
    stop("Assigned map geometry does not match the reference raster: ", file)
  }
  
  fr <- as.data.table(freq(x))
  if (nrow(fr) == 0) stop("Assigned map has no non-NA cells: ", file)
  value_col <- pick_col(fr, c("value", "zone", "zoneID"))
  if (is.na(value_col)) stop("Cannot identify raster values in frequency table: ", file)
  vals <- suppressWarnings(
    as.integer(round(as.numeric(as.character(fr[[value_col]]))))
  )
  vals <- sort(unique(vals[!is.na(vals)]))
  
  allowed <- if (scenario == "normal") model_zoneID else c(model_zoneID, novel_value)
  unexpected <- setdiff(vals, allowed)
  if (length(unexpected) > 0) {
    stop(
      "Assigned map contains unexpected zone values for ", method, " | ", scenario,
      ": ", paste(unexpected, collapse = ", "), " | ", file
    )
  }
  if (scenario == "normal" && novel_value %in% vals) {
    stop("Predicted reference map must not contain Novel 99: ", file)
  }
  
  file
}

find_assigned_map <- function(method, scenario) {
  cache_key <- paste(method, scenario, sep = "|")
  if (exists(cache_key, envir = assigned_map_cache, inherits = FALSE)) {
    return(get(cache_key, envir = assigned_map_cache, inherits = FALSE))
  }
  
  method_dir <- file.path(result_map_root, method)
  if (!dir.exists(method_dir)) return(NA_character_)
  
  expected_basenames <- c(
    paste0(
      "assigned_zone_", scenario,
      "_threshold0.2_tol1e-04_novel99_maskNA8_noNovelNormal.tif"
    ),
    paste0(
      "assigned_zone_", scenario,
      "_threshold0.2_tol0.0001_novel99_maskNA8_noNovelNormal.tif"
    )
  )
  expected_files <- file.path(method_dir, expected_basenames)
  existing <- expected_files[file.exists(expected_files)]
  
  if (length(existing) == 1) {
    selected <- validate_assigned_map_file(existing, method, scenario)
    assign(cache_key, selected, envir = assigned_map_cache)
    return(selected)
  }
  if (length(existing) > 1) {
    stop(
      "Multiple exact assigned maps found for ",
      method,
      " | ",
      scenario,
      ": ",
      paste(existing, collapse = "; ")
    )
  }
  
  # Strict fallback for equivalent formatting only. Never select threshold 0.1.
  files <- list.files(
    method_dir,
    pattern = paste0("^assigned_zone_", scenario, ".*\\.tif$"),
    recursive = FALSE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  files <- files[
    grepl("threshold0\\.2", basename(files), ignore.case = TRUE) &
      grepl("tol(1e-04|0\\.0001)", basename(files), ignore.case = TRUE) &
      grepl("novel99", basename(files), ignore.case = TRUE) &
      grepl("maskNA8", basename(files), ignore.case = TRUE) &
      grepl("noNovelNormal", basename(files), ignore.case = TRUE) &
      !grepl("(color|rgb|plot|legend)", basename(files), ignore.case = TRUE)
  ]
  
  if (length(files) == 1) {
    selected <- validate_assigned_map_file(files, method, scenario)
    assign(cache_key, selected, envir = assigned_map_cache)
    return(selected)
  }
  if (length(files) > 1) {
    stop(
      "Ambiguous threshold-0.2 assigned maps for ",
      method,
      " | ",
      scenario,
      ": ",
      paste(files, collapse = "; ")
    )
  }
  
  NA_character_
}

find_dual_file <- function(method, scenario, zone) {
  f <- file.path(
    dual_root,
    method,
    scenario,
    paste0("dual_suitability_zone", zone, ".tif")
  )
  
  if (file.exists(f)) f else NA_character_
}

area_by_zone <- function(file) {
  x <- rast(file)[[1]]
  a <- cellSize(x, unit = "km")
  z <- zonal(a, x, fun = "sum", na.rm = TRUE)
  dt <- as.data.table(z)
  
  if (ncol(dt) < 2) stop("Cannot calculate area by zone: ", file)
  
  setnames(dt, names(dt)[1], "zoneID")
  setnames(dt, names(dt)[2], "area_km2")
  dt[, zoneID := as.integer(round(zoneID))]
  dt <- dt[!is.na(zoneID)]
  dt <- dt[zoneID %in% c(model_zoneID, novel_value)]
  dt[]
}

zone_area_from_reference <- function() {
  if (!file.exists(reference_file)) stop("Missing reference raster: ", reference_file)
  out <- area_by_zone(reference_file)[zoneID %in% model_zoneID]
  assert_unique_keys(out, "zoneID", "Reference-map zone areas")
  missing_zones <- setdiff(model_zoneID, out$zoneID)
  if (length(missing_zones) > 0) {
    stop(
      "Reference raster is missing modeled zones required for area scaling: ",
      paste(missing_zones, collapse = ", ")
    )
  }
  if (any(!is.finite(out$area_km2) | out$area_km2 <= 0)) {
    stop("Reference-map zone areas must all be finite and positive.")
  }
  out[]
}

validate_common_prediction_mask <- function(files, label) {
  files <- as.character(files)
  if (length(files) < 2) return(invisible(TRUE))
  if (any(is.na(files)) || any(!file.exists(files))) {
    stop(label, " includes missing map files.")
  }
  
  ref_valid <- !is.na(rast(files[1])[[1]])
  for (i in 2:length(files)) {
    current <- rast(files[i])[[1]]
    if (!isTRUE(compareGeom(ref_valid, current, stopOnError = FALSE))) {
      stop(label, " contains geometrically inconsistent maps.")
    }
    mismatch <- global(ref_valid != !is.na(current), "sum", na.rm = TRUE)[1, 1]
    mismatch <- as.numeric(mismatch)
    if (!is.finite(mismatch) || mismatch > 0) {
      stop(
        label, " does not use one common prediction mask. ",
        "Mask mismatch cells between ", basename(files[1]), " and ",
        basename(files[i]), ": ", mismatch
      )
    }
  }
  invisible(TRUE)
}

reference_map_comparison_cache <- new.env(parent = emptyenv())

load_reference_category_lut <- function() {
  pal_file <- resolve_exact_file(
    file.path(base_dir, "color_palette_China.csv"),
    "vegetation-zone palette",
    required = TRUE
  )
  pal <- fread(pal_file)
  zone_col <- pick_col(pal, c("zoneID", "zone_id", "zoneid", "id"))
  category_col <- pick_col(
    pal,
    c("category2", "broad_category", "vegetation_category", "category")
  )
  if (is.na(zone_col) || is.na(category_col)) {
    stop(
      "color_palette_China.csv must contain zoneID and category2 (or an ",
      "equivalent broad-category column) for map reconstruction assessment."
    )
  }
  
  lut_dt <- data.table(
    zoneID = extract_zone_id(pal[[zone_col]]),
    category = trimws(as.character(pal[[category_col]]))
  )
  lut_dt <- lut_dt[zoneID %in% model_zoneID]
  assert_unique_keys(lut_dt, "zoneID", "Reference broad-category lookup")
  
  missing_zones <- setdiff(model_zoneID, lut_dt$zoneID)
  bad_categories <- lut_dt[is.na(category) | !nzchar(category)]
  if (length(missing_zones) > 0 || nrow(bad_categories) > 0) {
    stop(
      "Broad-category lookup is incomplete. Missing zones: ",
      paste(missing_zones, collapse = ", "),
      "; blank categories: ",
      paste(bad_categories$zoneID, collapse = ", ")
    )
  }
  
  setNames(lut_dt$category, as.character(lut_dt$zoneID))
}

reference_map_files <- function(methods) {
  methods <- as.character(methods)
  invalid <- setdiff(methods, all_map_method_order)
  if (length(invalid) > 0) {
    stop("Unsupported reference-map methods: ", paste(invalid, collapse = ", "))
  }
  
  files <- setNames(rep(NA_character_, length(methods)), methods)
  binary_methods <- intersect(methods, map_method_order)
  if (length(binary_methods) > 0) {
    files[binary_methods] <- vapply(
      binary_methods,
      function(m) find_assigned_map(m, "normal"),
      character(1)
    )
  }
  if ("multiclass_rf" %in% methods) {
    files["multiclass_rf"] <- validate_assigned_map_file(
      multiclass_reference_map_file,
      "multiclass_rf",
      "normal"
    )
  }
  
  missing <- names(files)[is.na(files) | !file.exists(files)]
  if (length(missing) > 0) {
    stop(
      "Reference-map comparison is missing assigned maps: ",
      paste(missing, collapse = ", ")
    )
  }
  files
}

build_reference_map_comparison <- function(
    methods,
    cache_key,
    output_prefix) {
  methods <- as.character(methods)
  if (exists(cache_key, envir = reference_map_comparison_cache, inherits = FALSE)) {
    return(get(cache_key, envir = reference_map_comparison_cache, inherits = FALSE))
  }
  
  if (!file.exists(reference_file)) stop("Missing reference raster: ", reference_file)
  files <- reference_map_files(methods)
  
  original_raw <- rast(reference_file)[[1]]
  original <- subst(
    original_raw,
    from = model_zoneID,
    to = model_zoneID,
    others = NA
  )
  names(original) <- "ori"
  
  predictions <- lapply(methods, function(m) {
    x <- rast(files[[m]])[[1]]
    if (!isTRUE(compareGeom(original, x, stopOnError = FALSE))) {
      stop("Reference-map geometry mismatch for method: ", m)
    }
    x <- subst(x, from = model_zoneID, to = model_zoneID, others = NA)
    names(x) <- m
    x
  })
  names(predictions) <- methods
  
  valid_original_pixels <- as.numeric(
    global(!is.na(original), "sum", na.rm = TRUE)[1, 1]
  )
  if (!is.finite(valid_original_pixels) || valid_original_pixels <= 0) {
    stop("The modeled-zone reference raster contains no valid pixels.")
  }
  
  individual_mask_audit <- rbindlist(lapply(methods, function(m) {
    compared <- as.numeric(
      global(!is.na(original) & !is.na(predictions[[m]]), "sum", na.rm = TRUE)[1, 1]
    )
    data.table(
      method_key = m,
      valid_original_pixels = valid_original_pixels,
      individual_compared_pixels = compared,
      individual_missing_predictions = valid_original_pixels - compared,
      individual_coverage = compared / valid_original_pixels,
      source_map = relative_path(files[[m]])
    )
  }))
  
  common_valid <- !is.na(original)
  for (m in methods) common_valid <- common_valid & !is.na(predictions[[m]])
  common_compared_pixels <- as.numeric(
    global(common_valid, "sum", na.rm = TRUE)[1, 1]
  )
  if (!is.finite(common_compared_pixels) || common_compared_pixels <= 0) {
    stop("No common valid pixels remain for reference-map comparison: ", cache_key)
  }
  
  original_common <- ifel(common_valid, original, NA)
  names(original_common) <- "ori"
  category_lut <- load_reference_category_lut()
  
  confusion_list <- list()
  zone_list <- list()
  overall_list <- list()
  
  for (m in methods) {
    pred_common <- ifel(common_valid, predictions[[m]], NA)
    names(pred_common) <- "pred"
    
    ct <- as.data.table(
      crosstab(c(original_common, pred_common), long = TRUE, useNA = FALSE)
    )
    if (ncol(ct) != 3) {
      stop("Unexpected crosstab structure for method: ", m)
    }
    setnames(ct, names(ct), c("ori", "pred", "n"))
    ct[, `:=`(
      method_key = m,
      ori = as.integer(ori),
      pred = as.integer(pred),
      n = as.numeric(n)
    )]
    ct <- ct[
      ori %in% model_zoneID & pred %in% model_zoneID &
        is.finite(n) & n > 0
    ]
    if (sum(ct$n) != common_compared_pixels) {
      stop(
        "Crosstab total does not equal the common comparison mask for ", m,
        ": ", sum(ct$n), " versus ", common_compared_pixels
      )
    }
    
    zone <- rbindlist(lapply(model_zoneID, function(z) {
      # Use lower-case local scalars and assign them explicitly to the output
      # columns. This avoids data.table non-standard-evaluation collisions with
      # the uppercase TP/TN/FP/FN column names.
      tp <- as.numeric(ct[ori == z & pred == z, sum(n, na.rm = TRUE)])
      fn <- as.numeric(ct[ori == z & pred != z, sum(n, na.rm = TRUE)])
      fp <- as.numeric(ct[ori != z & pred == z, sum(n, na.rm = TRUE)])
      tn <- as.numeric(common_compared_pixels - tp - fn - fp)
      
      if (any(!is.finite(c(tp, tn, fp, fn))) || any(c(tp, tn, fp, fn) < 0)) {
        stop(
          "Invalid confusion counts for method ", m,
          " and zone ", z,
          ": TP=", tp, ", TN=", tn, ", FP=", fp, ", FN=", fn
        )
      }
      
      original_pixels_z <- tp + fn
      predicted_pixels_z <- tp + fp
      
      data.table(
        method_key = as.character(m),
        zoneID = as.integer(z),
        original_pixels = original_pixels_z,
        predicted_pixels = predicted_pixels_z,
        predicted_to_original = if (
          is.finite(original_pixels_z) && original_pixels_z > 0
        ) {
          predicted_pixels_z / original_pixels_z
        } else {
          NA_real_
        },
        TP = tp,
        TN = tn,
        FP = fp,
        FN = fn
      )
    }))
    if (any(zone$original_pixels <= 0)) {
      missing_common_zones <- zone[original_pixels <= 0, zoneID]
      stop(
        "The common comparison mask removes all original pixels for zones: ",
        paste(missing_common_zones, collapse = ", "),
        " | method set: ", cache_key
      )
    }
    zone <- recompute_classification_metrics(
      zone,
      paste0("Direct reference-map zone metrics: ", m)
    )
    validate_confusion_derived_counts(
      zone,
      paste0("Direct reference-map zone metrics: ", m)
    )
    if (any(!is.finite(zone$f1))) {
      stop("Non-finite direct reference-map F1 for method: ", m)
    }
    
    ct[, `:=`(
      ori_category = unname(category_lut[as.character(ori)]),
      pred_category = unname(category_lut[as.character(pred)])
    )]
    if (anyNA(ct$ori_category) || anyNA(ct$pred_category)) {
      stop("Missing broad-category labels in direct map assessment for: ", m)
    }
    
    mask_row <- individual_mask_audit[method_key == m]
    overall <- data.table(
      method_key = m,
      valid_original_pixels = valid_original_pixels,
      individual_compared_pixels = mask_row$individual_compared_pixels,
      individual_missing_predictions = mask_row$individual_missing_predictions,
      individual_coverage = mask_row$individual_coverage,
      compared_pixels = common_compared_pixels,
      missing_predictions = valid_original_pixels - common_compared_pixels,
      coverage = common_compared_pixels / valid_original_pixels,
      exact_zone_accuracy = ct[ori == pred, sum(n)] / common_compared_pixels,
      broad_category_accuracy = ct[ori_category == pred_category, sum(n)] /
        common_compared_pixels,
      macro_balanced_accuracy = mean(zone$balanced_accuracy),
      macro_recall = mean(zone$recall),
      macro_specificity = mean(zone$specificity),
      macro_precision = mean(zone$precision),
      macro_f1 = mean(zone$f1),
      macro_tss = mean(zone$tss),
      source_map = relative_path(files[[m]]),
      comparison_scope = cache_key
    )
    
    confusion_list[[m]] <- ct
    zone_list[[m]] <- zone
    overall_list[[m]] <- overall
  }
  
  confusion <- rbindlist(confusion_list, fill = TRUE)
  zone <- rbindlist(zone_list, fill = TRUE)
  overall <- rbindlist(overall_list, fill = TRUE)
  
  require_complete_zone_grid(
    zone,
    methods,
    model_zoneID,
    label = paste0("Direct reference-map metrics: ", cache_key)
  )
  validate_metric_ranges(
    zone,
    unit_cols = c(
      "accuracy", "balanced_accuracy", "recall", "specificity",
      "precision", "f1"
    ),
    tss_cols = "tss",
    label = paste0("Direct reference-map zone metrics: ", cache_key)
  )
  validate_metric_ranges(
    overall,
    unit_cols = c(
      "individual_coverage", "coverage", "exact_zone_accuracy",
      "broad_category_accuracy", "macro_balanced_accuracy", "macro_recall",
      "macro_specificity", "macro_precision", "macro_f1"
    ),
    tss_cols = "macro_tss",
    label = paste0("Direct reference-map overall metrics: ", cache_key)
  )
  validate_overall_against_zone_counts(
    overall,
    zone,
    paste0("Direct reference-map metrics: ", cache_key)
  )
  
  overall[, method_label := all_map_method_labels[method_key]]
  zone[, method_label := all_map_method_labels[method_key]]
  confusion[, method_label := all_map_method_labels[method_key]]
  
  output_files <- c(
    overall = file.path(tab_dir, paste0(output_prefix, "_overall_metrics.csv")),
    zone = file.path(tab_dir, paste0(output_prefix, "_zone_metrics.csv")),
    confusion = file.path(tab_dir, paste0(output_prefix, "_confusion_long.csv")),
    mask = file.path(tab_dir, paste0(output_prefix, "_mask_audit.csv"))
  )
  remove_existing_outputs(output_files)
  fwrite(overall, output_files[["overall"]])
  fwrite(zone, output_files[["zone"]])
  fwrite(confusion, output_files[["confusion"]])
  fwrite(
    overall[, .(
      method_key, method_label, valid_original_pixels,
      individual_compared_pixels, individual_missing_predictions,
      individual_coverage, compared_pixels, missing_predictions, coverage,
      source_map, comparison_scope
    )],
    output_files[["mask"]]
  )
  
  result <- list(
    overall = overall,
    zone = zone,
    confusion = confusion,
    files = files,
    output_files = output_files
  )
  assign(cache_key, result, envir = reference_map_comparison_cache)
  result
}

get_binary_reference_comparison <- function() {
  build_reference_map_comparison(
    methods = map_method_order,
    cache_key = "binary_common_mask",
    output_prefix = "Binary_reference_map_recalculated"
  )
}

get_unified_reference_comparison <- function() {
  build_reference_map_comparison(
    methods = all_map_method_order,
    cache_key = "binary_multiclass_common_mask",
    output_prefix = "Unified_reference_map_recalculated"
  )
}

validate_unified_reference_map_masks <- function() {
  invisible(get_unified_reference_comparison())
}

read_binary_map_overall <- function() {
  out <- copy(get_binary_reference_comparison()$overall)
  out <- out[method_key %in% map_method_order]
  out[, source_file := source_map]
  out[, method_order_tmp := match(method_key, map_method_order)]
  setorder(out, method_order_tmp)
  out[, method_order_tmp := NULL]
  out[]
}

read_binary_map_zone <- function() {
  out <- copy(get_binary_reference_comparison()$zone)
  out <- out[method_key %in% map_method_order]
  source_lut <- setNames(
    get_binary_reference_comparison()$overall$source_map,
    get_binary_reference_comparison()$overall$method_key
  )
  out[, source_file := unname(source_lut[method_key])]
  out[, method_order_tmp := match(method_key, map_method_order)]
  setorder(out, method_order_tmp, zoneID)
  out[, method_order_tmp := NULL]
  out[]
}

read_population_lookup <- function() {
  f <- resolve_exact_file(
    population_lookup_file,
    "full reference population lookup generated by code 6",
    required = TRUE
  )
  dt <- fread(f)
  
  species_col <- pick_col(dt, c("Species", "species", "species_name", "taxon"))
  zone_col <- pick_col(dt, c("source_zone", "population_zone", "zoneID", "zone_id", "Zone"))
  abundance_col <- pick_col(
    dt,
    c("reference_abundance", "Population", "freq_pop", "frequency", "abundance")
  )
  projected_col <- pick_col(dt, c("projected", "is_projected"))
  
  if (is.na(species_col) || is.na(zone_col)) {
    stop("Cannot identify Species and source_zone columns in: ", f)
  }
  
  out <- data.table(
    species = trimws(as.character(dt[[species_col]])),
    source_zone = extract_zone_id(dt[[zone_col]])
  )
  
  if (!is.na(abundance_col)) {
    out[, reference_abundance := as.numeric(dt[[abundance_col]])]
  } else {
    out[, reference_abundance := NA_real_]
  }
  
  if (!is.na(projected_col)) {
    out[, projected := as.logical(dt[[projected_col]])]
  } else {
    out[, projected := source_zone %in% model_zoneID]
  }
  
  out[, source_file := relative_path(f)]
  out <- out[
    !is.na(species) &
      nzchar(species) &
      source_zone %in% reference_zoneID
  ]
  
  if (nrow(out) == 0) {
    stop("No valid species-zone populations found in: ", f)
  }
  
  if (any(out$reference_abundance < 0, na.rm = TRUE)) {
    stop("Negative reference abundance values were found in: ", f)
  }
  low_retained <- out[
    !is.na(reference_abundance) &
      is.finite(reference_abundance) &
      reference_abundance < min_reference_population_abundance
  ]
  if (nrow(low_retained) > 0) {
    stop(
      "Population lookup contains abundances below the retained-population threshold (",
      min_reference_population_abundance, " occupied cells). Rerun code 1.2 and code 6. ",
      "First rows: ",
      paste(capture.output(print(head(low_retained, 12))), collapse = " ")
    )
  }
  if (anyNA(out$projected)) {
    stop("Missing projection-eligibility flags were found in: ", f)
  }
  projected_expected <- out$source_zone %in% model_zoneID
  if (any(out$projected != projected_expected)) {
    bad <- out[projected != (source_zone %in% model_zoneID)]
    stop(
      "Population projection-eligibility flags do not match the modeled-zone definition. ",
      "First rows: ",
      paste(capture.output(print(head(bad, 12))), collapse = " ")
    )
  }
  
  assert_unique_keys(
    out,
    c("species", "source_zone"),
    paste0("Population lookup: ", relative_path(f))
  )
  
  setorder(out, species, source_zone)
  
  # When code 6.2 has already written its projection-eligible lookup, verify that
  # it is exactly the modeled-source-zone subset of the full code 6 table. This catches stale
  # population definitions before Figures 8-9 are generated.
  if (file.exists(dual_population_lookup_file)) {
    dual_dt <- fread(dual_population_lookup_file)
    sp_col <- pick_col(dual_dt, c("Species", "species"))
    z_col <- pick_col(dual_dt, c("source_zone", "zoneID", "zone_id", "Zone"))
    a_col <- pick_col(
      dual_dt,
      c("reference_abundance", "Population", "freq_pop", "frequency", "abundance")
    )
    if (is.na(sp_col) || is.na(z_col) || is.na(a_col)) {
      stop(
        "Cannot audit code 6.2 population lookup; required columns are missing: ",
        dual_population_lookup_file
      )
    }
    
    dual_norm <- data.table(
      species = trimws(as.character(dual_dt[[sp_col]])),
      source_zone = extract_zone_id(dual_dt[[z_col]]),
      reference_abundance = as.numeric(dual_dt[[a_col]])
    )
    dual_norm <- dual_norm[source_zone %in% model_zoneID]
    assert_unique_keys(
      dual_norm,
      c("species", "source_zone"),
      "Code 6.2 projection-eligible population lookup"
    )
    
    code6_projected <- out[
      projected == TRUE & source_zone %in% model_zoneID,
      .(species, source_zone, reference_abundance)
    ]
    cmp <- merge(
      code6_projected,
      dual_norm,
      by = c("species", "source_zone"),
      all = TRUE,
      suffixes = c("_code6", "_code6_2")
    )
    bad <- cmp[
      is.na(reference_abundance_code6) |
        is.na(reference_abundance_code6_2) |
        abs(reference_abundance_code6 - reference_abundance_code6_2) > 1e-10
    ]
    if (nrow(bad) > 0) {
      stop(
        "Code 6 and code 6.2 population lookups are inconsistent. Rerun code 6.2. ",
        "First mismatches: ",
        paste(capture.output(print(head(bad, 12))), collapse = " ")
      )
    }
  }
  
  out[]
}

get_species_population_table <- function() {
  pop <- read_population_lookup()
  pop <- pop[
    projected == TRUE &
      source_zone %in% model_zoneID,
    .(species, source_zone, reference_abundance, source_file)
  ]
  
  if (nrow(pop) == 0) {
    stop("No projection-eligible species-zone populations remain after excluding Zones 8 and 51.")
  }
  
  pop[]
}

dual_to_points <- function(file, source_zone, threshold = dual_plot_threshold) {
  if (is.na(file) || !file.exists(file)) return(data.table())
  
  x <- rast(file)[[1]]
  ref <- rast(reference_file)[[1]]
  
  rng <- global(x, c("min", "max"), na.rm = TRUE)
  min_v <- as.numeric(rng[1, "min"])
  max_v <- as.numeric(rng[1, "max"])
  if (!is.finite(min_v) || !is.finite(max_v)) {
    stop("Dual-suitability raster has no finite values: ", file)
  }
  if (min_v < -1e-6 || max_v > 1 + 1e-6) {
    stop(
      "Dual-suitability values fall outside [0, 1]: ", file,
      " | range = [", min_v, ", ", max_v, "]"
    )
  }
  
  if (!isTRUE(compareGeom(x, ref, stopOnError = FALSE))) {
    x <- resample(x, ref, method = "bilinear")
  }
  x <- mask(x, ref)
  
  # Apply the scientific suitability threshold before display aggregation.
  # Averaging the full raster first can erase small suitable patches and shrink
  # the displayed niche. Max aggregation preserves every suitable source cell
  # while retaining the maximum dual suitability for the display cell.
  x <- ifel(!is.na(x) & x > threshold, x, NA)
  if (ncell(x) > display_max_cells) {
    fact <- ceiling(sqrt(ncell(x) / display_max_cells))
    x <- aggregate(x, fact = fact, fun = "max", na.rm = TRUE)
    x <- ifel(!is.na(x) & x > threshold, x, NA)
  }
  
  dt <- as.data.table(as.data.frame(x, xy = TRUE, na.rm = TRUE))
  
  if (nrow(dt) == 0) return(data.table())
  
  val_col <- setdiff(names(dt), c("x", "y"))[1]
  setnames(dt, val_col, "dual_suitability")
  
  dt <- dt[
    !is.na(dual_suitability) &
      is.finite(dual_suitability) &
      dual_suitability > threshold
  ]
  
  if (nrow(dt) == 0) return(data.table())
  
  dt[, source_zone := source_zone]
  dt[]
}

reference_plot_background <- function(max_cells = 120000L) {
  if (!file.exists(reference_file)) stop("Missing reference raster: ", reference_file)
  r <- rast(reference_file)[[1]]
  r_plot <- thin_raster_for_plot(r, categorical = TRUE, max_cells = max_cells)
  bg <- as.data.table(as.data.frame(r_plot, xy = TRUE, na.rm = TRUE))
  if (nrow(bg) == 0) stop("Reference raster has no valid cells for map backgrounds.")
  bg <- bg[, .(x, y)]
  list(
    data = bg,
    xlim = c(xmin(r), xmax(r)),
    ylim = c(ymin(r), ymax(r))
  )
}

build_dual_cache <- function(method, scenario, zones) {
  zones <- sort(unique(as.integer(zones)))
  files <- setNames(
    vapply(zones, function(z) find_dual_file(method, scenario, z), character(1)),
    as.character(zones)
  )
  
  missing_zones <- names(files)[is.na(files) | !file.exists(files)]
  if (length(missing_zones) > 0) {
    stop(
      "Incomplete dual-suitability inputs for ", method, " | ", scenario,
      ". Missing source zones: ", paste(missing_zones, collapse = ", ")
    )
  }
  
  cache <- lapply(zones, function(z) {
    dual_to_points(files[[as.character(z)]], z, threshold = dual_plot_threshold)
  })
  names(cache) <- as.character(zones)
  cache
}

validate_dual_scenario_grid <- function(method, scenarios, zones) {
  jobs <- CJ(
    scenario = as.character(scenarios),
    source_zone = as.integer(zones),
    unique = TRUE
  )
  jobs[, file := mapply(
    function(s, z) find_dual_file(method, s, z),
    scenario,
    source_zone,
    USE.NAMES = FALSE
  )]
  missing <- jobs[is.na(file) | !file.exists(file)]
  if (nrow(missing) > 0) {
    stop(
      "Incomplete dual-suitability method-scenario-zone grid for ", method, ". ",
      "First missing inputs: ",
      paste(capture.output(print(head(missing, 20))), collapse = " ")
    )
  }
  invisible(jobs)
}

plot_species_population_dual_pages <- function(pop_tbl, method, scenario, outfile_prefix) {
  species_order <- sort(unique(pop_tbl$species))
  zones_needed <- sort(unique(pop_tbl$source_zone))
  cache <- build_dual_cache(method, scenario, zones_needed)
  map_bg <- reference_plot_background()
  
  page_size <- species_page_ncol * species_page_nrow
  page_id <- ceiling(seq_along(species_order) / page_size)
  cols <- zone_color_vector(zones_needed)
  
  for (pg in sort(unique(page_id))) {
    sp_sel <- species_order[page_id == pg]
    plot_list <- list()
    
    for (sp in sp_sel) {
      z_sp <- pop_tbl[species == sp, sort(unique(source_zone))]
      dts <- cache[as.character(z_sp)]
      dts <- dts[lengths(dts) > 0]
      
      if (length(dts) > 0) {
        dt <- rbindlist(dts, fill = TRUE)
        if (nrow(dt) > 0) {
          dt[, species := sp]
          plot_list[[length(plot_list) + 1L]] <- dt
        }
      }
    }
    
    if (length(plot_list) > 0) {
      plot_dt <- rbindlist(plot_list, fill = TRUE)
    } else {
      plot_dt <- data.table(
        x = numeric(), y = numeric(), dual_suitability = numeric(),
        source_zone = integer(), species = character()
      )
    }
    plot_dt[, species := factor(species, levels = sp_sel)]
    plot_dt[, source_zone_chr := factor(
      as.character(source_zone),
      levels = as.character(zones_needed)
    )]
    
    panel_frame <- data.table(
      species = factor(sp_sel, levels = sp_sel),
      x = mean(map_bg$xlim),
      y = mean(map_bg$ylim)
    )
    
    p <- ggplot() +
      geom_raster(
        data = map_bg$data,
        aes(x = x, y = y),
        inherit.aes = FALSE,
        fill = "grey95"
      ) +
      geom_blank(data = panel_frame, aes(x = x, y = y)) +
      geom_point(
        data = plot_dt,
        aes(x = x, y = y, colour = source_zone_chr),
        size = 0.075,
        alpha = 0.30
      ) +
      facet_wrap(~ species, ncol = species_page_ncol, drop = FALSE) +
      coord_equal(xlim = map_bg$xlim, ylim = map_bg$ylim, expand = FALSE) +
      scale_colour_manual(values = cols, drop = FALSE) +
      labs(
        title = paste0("Tree population dual-niche projections: ", map_method_labels[method], " | ", scenario),
        subtitle = paste0(
          "Full reference-map extent. Each population uses its source-zone dual niche; ",
          "cells with dual suitability > ", dual_plot_threshold,
          " are colored by source zone with 30% opacity."
        ),
        colour = "Source zone"
      ) +
      theme_void(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0, size = 13),
        plot.subtitle = element_text(hjust = 0, size = 9),
        strip.text = element_text(face = "bold", size = 9.2),
        legend.position = "right",
        legend.text = element_text(size = 6.6),
        legend.key.size = grid::unit(0.30, "cm")
      ) +
      guides(colour = guide_legend(
        ncol = 2,
        byrow = TRUE,
        override.aes = list(size = 2, alpha = 1)
      ))
    
    outfile <- paste0(
      outfile_prefix, "_", safe_name(method), "_", safe_name(scenario),
      "_page", sprintf("%02d", pg), ".png"
    )
    ggsave(outfile, p, width = 15.5, height = 7.0, dpi = fig_dpi, bg = "white")
    cat0("[SAVED] ", outfile)
  }
}

plot_species_level_dual_pages <- function(pop_tbl, method, scenario, outfile_prefix) {
  species_order <- sort(unique(pop_tbl$species))
  zones_needed <- sort(unique(pop_tbl$source_zone))
  cache <- build_dual_cache(method, scenario, zones_needed)
  map_bg <- reference_plot_background()
  
  page_size <- species_page_ncol * species_page_nrow
  page_id <- ceiling(seq_along(species_order) / page_size)
  
  for (pg in sort(unique(page_id))) {
    sp_sel <- species_order[page_id == pg]
    plot_list <- list()
    
    for (sp in sp_sel) {
      z_sp <- pop_tbl[species == sp, sort(unique(source_zone))]
      dts <- cache[as.character(z_sp)]
      dts <- dts[lengths(dts) > 0]
      
      if (length(dts) > 0) {
        dt <- rbindlist(dts, fill = TRUE)
        if (nrow(dt) > 0) {
          dt <- dt[
            ,
            .(dual_suitability = max(dual_suitability, na.rm = TRUE)),
            by = .(x, y)
          ]
          dt[, species := sp]
          plot_list[[length(plot_list) + 1L]] <- dt
        }
      }
    }
    
    if (length(plot_list) > 0) {
      plot_dt <- rbindlist(plot_list, fill = TRUE)
    } else {
      plot_dt <- data.table(
        x = numeric(), y = numeric(), dual_suitability = numeric(),
        species = character()
      )
    }
    plot_dt[, species := factor(species, levels = sp_sel)]
    
    panel_frame <- data.table(
      species = factor(sp_sel, levels = sp_sel),
      x = mean(map_bg$xlim),
      y = mean(map_bg$ylim)
    )
    
    p <- ggplot() +
      geom_raster(
        data = map_bg$data,
        aes(x = x, y = y),
        inherit.aes = FALSE,
        fill = "grey95"
      ) +
      geom_blank(data = panel_frame, aes(x = x, y = y)) +
      geom_point(
        data = plot_dt,
        aes(x = x, y = y, colour = dual_suitability),
        size = 0.075,
        alpha = 0.30
      ) +
      facet_wrap(~ species, ncol = species_page_ncol, drop = FALSE) +
      coord_equal(xlim = map_bg$xlim, ylim = map_bg$ylim, expand = FALSE) +
      scale_colour_gradientn(
        colours = grDevices::hcl.colors(60, palette = "YlGnBu"),
        limits = c(dual_plot_threshold, 1),
        oob = scales::squish,
        name = "Max dual\nsuitability"
      ) +
      labs(
        title = paste0("Species-level dual-niche projections: ", map_method_labels[method], " | ", scenario),
        subtitle = paste0(
          "Full reference-map extent. Each species is the pixel-wise maximum across ",
          "its population dual niches; cells with suitability > ",
          dual_plot_threshold, " are shown with 30% opacity."
        )
      ) +
      theme_void(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0, size = 13),
        plot.subtitle = element_text(hjust = 0, size = 9),
        strip.text = element_text(face = "bold", size = 9.2),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.key.size = grid::unit(0.35, "cm")
      )
    
    outfile <- paste0(
      outfile_prefix, "_", safe_name(method), "_", safe_name(scenario),
      "_page", sprintf("%02d", pg), ".png"
    )
    ggsave(outfile, p, width = 15.5, height = 7.0, dpi = fig_dpi, bg = "white")
    cat0("[SAVED] ", outfile)
  }
}

# 6. Managed-output cleanup and priority execution =========================================================

# Remove every output owned by this script before the strict run begins. This
# prevents an old, invalid PNG/CSV from surviving when a later step stops. The
# patterns are deliberately restricted to this script's Figure/Table products.
if (purge_previous_outputs) {
  previous_outputs <- c(
    list.files(
      fig_dir,
      pattern = "^Figure([1-9]|10).*\\.png$",
      full.names = TRUE,
      recursive = FALSE,
      ignore.case = TRUE
    ),
    list.files(
      tab_dir,
      pattern = paste0(
        "^(Table[1-4]|Figure([1-9]|10)|",
        "Binary_reference_map_recalculated|",
        "Unified_reference_map_recalculated).*\\.csv$"
      ),
      full.names = TRUE,
      recursive = FALSE,
      ignore.case = TRUE
    ),
    file.path(vis_dir, "visualization_step_log.csv")
  )
  previous_outputs <- unique(previous_outputs[
    file.exists(previous_outputs) & !dir.exists(previous_outputs)
  ])
  if (length(previous_outputs) > 0) {
    unlink(previous_outputs, force = TRUE)
    cat0(
      "[PURGED MANAGED VISUALIZATION OUTPUTS] ",
      length(previous_outputs),
      " files"
    )
  }
}

# Priority output helpers | Unified reference-map comparison ===================

read_map_comparison_overall <- function() {
  out <- copy(get_unified_reference_comparison()$overall)
  out[, source_file := source_map]
  out[, method_order_tmp := match(method_key, all_map_method_order)]
  setorder(out, method_order_tmp)
  out[, method_order_tmp := NULL]
  out[]
}

read_map_comparison_zone <- function() {
  out <- copy(get_unified_reference_comparison()$zone)
  source_lut <- setNames(
    get_unified_reference_comparison()$overall$source_map,
    get_unified_reference_comparison()$overall$method_key
  )
  out[, source_file := unname(source_lut[method_key])]
  out[, method_order_tmp := match(method_key, all_map_method_order)]
  setorder(out, method_order_tmp, zoneID)
  out[, method_order_tmp := NULL]
  out[]
}

audit_map_metric_consistency <- function() {
  binary <- read_binary_map_overall()[
    , .(
      method_key,
      binary_scope_compared_pixels = compared_pixels,
      binary_scope_coverage = coverage,
      binary_scope_exact_zone_accuracy = exact_zone_accuracy,
      binary_scope_macro_f1 = macro_f1
    )
  ]
  unified <- read_map_comparison_overall()[
    method_key %in% map_method_order,
    .(
      method_key,
      unified_scope_compared_pixels = compared_pixels,
      unified_scope_coverage = coverage,
      unified_scope_exact_zone_accuracy = exact_zone_accuracy,
      unified_scope_macro_f1 = macro_f1,
      individual_coverage,
      source_map
    )
  ]
  audit <- merge(binary, unified, by = "method_key", all = TRUE, sort = FALSE)
  audit[, `:=`(
    comparison_note = paste0(
      "Binary-scope values use the common valid mask of four binary maps; ",
      "unified-scope values use the common valid mask of four binary maps plus multiclass."
    ),
    compared_pixel_difference =
      unified_scope_compared_pixels - binary_scope_compared_pixels,
    exact_accuracy_difference =
      unified_scope_exact_zone_accuracy - binary_scope_exact_zone_accuracy,
    macro_f1_difference = unified_scope_macro_f1 - binary_scope_macro_f1,
    saved_code5_overall_exists = file.exists(binary_map_overall_file),
    saved_code5_zone_exists = file.exists(binary_map_zone_file),
    saved_code7_overall_exists = file.exists(map_comparison_overall_file),
    saved_code7_zone_exists = file.exists(map_comparison_zone_file),
    manuscript_values_source = "Recalculated directly from raster maps"
  )]
  audit[, method_order_tmp := match(method_key, map_method_order)]
  setorder(audit, method_order_tmp)
  audit[, method_order_tmp := NULL]
  audit[]
}

# PRIORITY OUTPUT 1. Figure 10a ================================================

run_step("Figure 10a | Binary workflows vs multiclass reference-map F1 comparison", {
  figure10a_png <- file.path(
    fig_dir,
    "Figure10a_all_binary_workflows_vs_multiclass_reference_map_F1_area_scaled.png"
  )
  figure10a_csv <- file.path(
    tab_dir,
    "Figure10a_all_binary_workflows_vs_multiclass_reference_map_F1.csv"
  )
  remove_existing_outputs(c(figure10a_png, figure10a_csv))
  
  map_zone <- read_map_comparison_zone()
  
  binary_f1 <- map_zone[
    method_key %in% map_method_order,
    .(method_key, zoneID, binary_f1 = f1)
  ]
  multiclass_f1 <- map_zone[
    method_key == "multiclass_rf",
    .(zoneID, multiclass_f1 = f1)
  ]
  
  cmp <- merge(
    binary_f1,
    multiclass_f1,
    by = "zoneID",
    all.x = TRUE,
    sort = FALSE
  )
  
  area_dt <- zone_area_from_reference()[, .(zoneID, area_km2)]
  cmp <- merge(cmp, area_dt, by = "zoneID", all.x = TRUE, sort = FALSE)
  
  cmp[, `:=`(
    method_label = map_method_labels[method_key],
    zoneID_chr = as.character(zoneID),
    valid_pair = is.finite(binary_f1) & is.finite(multiclass_f1),
    f1_difference = multiclass_f1 - binary_f1,
    source_file = "Recalculated directly from original and assigned-zone rasters"
  )]
  cmp[, higher_f1 := fifelse(
    !valid_pair,
    NA_character_,
    fifelse(
      f1_difference > 0,
      "Multiclass",
      fifelse(f1_difference < 0, "Binary", "Equal")
    )
  )]
  cmp[, method_label := factor(
    method_label,
    levels = map_method_labels[map_method_order]
  )]
  
  assert_unique_keys(
    cmp,
    c("method_key", "zoneID"),
    "Figure 10a comparison data"
  )
  require_complete_zone_grid(
    cmp,
    map_method_order,
    model_zoneID,
    label = "Figure 10a comparison data"
  )
  
  if (any(!cmp$valid_pair) || any(!is.finite(cmp$area_km2))) {
    bad <- cmp[!valid_pair | !is.finite(area_km2)]
    stop(
      "Figure 10a requires finite F1 and reference area for every method-zone pair. ",
      "First invalid rows: ",
      paste(capture.output(print(head(bad, 12))), collapse = " ")
    )
  }
  plot_cmp <- cmp
  
  cols <- zone_color_vector(plot_cmp$zoneID)
  axis_lim <- c(0, 1)
  
  p <- ggplot(
    plot_cmp,
    aes(
      x = binary_f1,
      y = multiclass_f1,
      colour = zoneID_chr,
      size = area_km2
    )
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      colour = "grey45",
      linewidth = 0.45
    ) +
    geom_point(alpha = 0.80) +
    geom_text(
      aes(label = zoneID),
      check_overlap = TRUE,
      size = 2.2,
      nudge_y = 0.012,
      alpha = 1,
      show.legend = FALSE
    ) +
    facet_wrap(~ method_label, nrow = 2) +
    scale_colour_manual(values = cols, guide = "none") +
    scale_size_continuous(
      range = c(1.2, 6.0),
      name = expression("Original area (km"^2*")")
    ) +
    scale_x_continuous(
      limits = axis_lim,
      breaks = seq(0, 1, by = 0.2),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      limits = axis_lim,
      breaks = seq(0, 1, by = 0.2),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    coord_equal(xlim = axis_lim, ylim = axis_lim, expand = FALSE) +
    labs(
      title = "Binary workflows vs multiclass reference-map F1",
      subtitle = paste0(
        "Both axes are zone-level reconstruction F1 values recalculated from the ",
        "same original raster and the common valid mask of all five maps. Above the line ",
        "favors multiclass; below the line favors binary."
      ),
      x = "Binary workflow reference-map F1",
      y = "Multiclass reference-map F1"
    ) +
    theme_ms(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  
  save_gg(
    p,
    figure10a_png,
    width = 9.5,
    height = 8
  )
  
  fwrite(
    cmp,
    figure10a_csv
  )
})

# 2. Table 1 and Figure 1 ======================================================

run_step("Table 1 and Figure 1 | One-hot climate/soil binary RF performance", {
  remove_existing_outputs(c(
    file.path(tab_dir, "Table1_source_files_used.csv"),
    file.path(tab_dir, "Table1_binary_RF_zone_level_metrics_long.csv"),
    file.path(tab_dir, "Table1_binary_RF_zone_level_metrics.csv"),
    file.path(tab_dir, "Table1_binary_RF_workflow_summary_long.csv"),
    file.path(tab_dir, "Table1_binary_RF_workflow_summary.csv"),
    file.path(tab_dir, "Table1_metric_grid_audit.csv"),
    file.path(fig_dir, "Figure1a_onehot_binary_RF_summary.png"),
    file.path(
      fig_dir,
      "Figure1_onehot_binary_RF_OOB_train_test_performance_dotrange.png"
    ),
    file.path(fig_dir, "Figure1b_onehot_binary_RF_zone_level_metrics.png")
  ))
  
  tbl <- make_binary_rf_table1()
  long <- tbl$long
  summary_long <- tbl$summary_long
  
  long_file <- file.path(tab_dir, "Table1_binary_RF_zone_level_metrics_long.csv")
  fwrite(long, long_file)
  cat0("[SAVED] ", long_file)
  
  zone_wide <- dcast(
    long,
    niche_type + workflow + zoneID ~ metric,
    value.var = "value"
  )
  setorder(zone_wide, niche_type, workflow, zoneID)
  
  zone_file <- file.path(tab_dir, "Table1_binary_RF_zone_level_metrics.csv")
  fwrite(zone_wide, zone_file)
  cat0("[SAVED] ", zone_file)
  
  summary_file <- file.path(tab_dir, "Table1_binary_RF_workflow_summary_long.csv")
  fwrite(summary_long, summary_file)
  cat0("[SAVED] ", summary_file)
  
  table_metrics <- c(
    "OOB accuracy",
    "Train accuracy",
    "Test accuracy",
    "Test balanced accuracy",
    "Test AUC",
    "Test F1",
    "Test precision",
    "Test recall",
    "Test specificity",
    "Test TSS"
  )
  
  table_main <- summary_long[metric %in% table_metrics]
  table_main[, mean_se := sprintf("%.3f +/- %.3f", mean, se)]
  table_main <- dcast(
    table_main,
    niche_type + workflow + workflow_label + n_zones ~ metric,
    value.var = "mean_se"
  )
  setorder(table_main, niche_type, workflow)
  
  main_file <- file.path(tab_dir, "Table1_binary_RF_workflow_summary.csv")
  fwrite(table_main, main_file)
  cat0("[SAVED] ", main_file)
  
  metric_label_map <- c(
    "accuracy" = "Accuracy",
    "balanced accuracy" = "Balanced accuracy",
    "AUC" = "AUC",
    "F1" = "F1",
    "precision" = "Precision",
    "recall" = "Recall",
    "specificity" = "Specificity",
    "TSS" = "TSS"
  )
  
  metric_levels <- c(
    "Accuracy",
    "Balanced accuracy",
    "AUC",
    "F1",
    "Precision",
    "Recall",
    "Specificity",
    "TSS"
  )
  
  plot_summary <- copy(summary_long)
  plot_summary <- plot_summary[
    grepl("^(OOB|Train|Test)", metric) &
      !grepl("error|threshold", metric, ignore.case = TRUE)
  ]
  
  if (nrow(plot_summary) == 0) stop("No OOB/Train/Test metrics found.")
  
  plot_summary[, eval_set := sub(" .*", "", metric)]
  plot_summary[, metric_name := sub("^(OOB|Train|Test) ", "", metric)]
  plot_summary[, metric_label := unname(metric_label_map[metric_name])]
  plot_summary <- plot_summary[!is.na(metric_label)]
  
  if (nrow(plot_summary) == 0) stop("No supported Figure 1 metrics found.")
  
  plot_summary[, metric_label := factor(metric_label, levels = rev(metric_levels))]
  plot_summary[, eval_set := factor(eval_set, levels = c("OOB", "Train", "Test"))]
  plot_summary[, workflow := factor(workflow, levels = workflow_order)]
  plot_summary[, workflow_label := factor(
    workflow_labels[as.character(workflow)],
    levels = workflow_labels[workflow_order]
  )]
  plot_summary[, niche_type := factor(niche_type, levels = c("climate", "soil"))]
  
  plot_zone <- copy(long)
  plot_zone <- plot_zone[
    grepl("^(OOB|Train|Test)", metric) &
      !grepl("error|threshold", metric, ignore.case = TRUE)
  ]
  plot_zone[, eval_set := sub(" .*", "", metric)]
  plot_zone[, metric_name := sub("^(OOB|Train|Test) ", "", metric)]
  plot_zone[, metric_label := unname(metric_label_map[metric_name])]
  plot_zone <- plot_zone[!is.na(metric_label)]
  plot_zone[, metric_label := factor(metric_label, levels = rev(metric_levels))]
  plot_zone[, eval_set := factor(eval_set, levels = c("OOB", "Train", "Test"))]
  plot_zone[, workflow := factor(workflow, levels = workflow_order)]
  plot_zone[, workflow_label := factor(
    workflow_labels[as.character(workflow)],
    levels = workflow_labels[workflow_order]
  )]
  plot_zone[, niche_type := factor(niche_type, levels = c("climate", "soil"))]
  plot_zone[, zoneID_chr := factor(as.character(zoneID), levels = as.character(model_zoneID))]
  
  zone_cols <- zone_color_vector(model_zoneID)
  
  x_all <- c(
    plot_zone$value,
    plot_summary$ci95_low,
    plot_summary$ci95_high,
    plot_summary$mean
  )
  x_all <- x_all[is.finite(x_all)]
  if (length(x_all) == 0) stop("No finite values are available for plotting.")
  
  # Use the requested zoom only when it does not hide a point or confidence interval.
  observed_lower <- floor(min(x_all, na.rm = TRUE) / 0.05) * 0.05
  x_lower_fig1a <- max(-1, min(0.6, observed_lower))
  x_lower_fig1b <- max(-1, min(0.4, observed_lower))
  x_upper <- min(1, max(1, ceiling(max(x_all, na.rm = TRUE) / 0.05) * 0.05))
  if (x_upper <= x_lower_fig1a) x_upper <- 1
  x_breaks_fig1a <- seq(x_lower_fig1a, x_upper, by = 0.1)
  x_breaks_fig1b <- seq(x_lower_fig1b, x_upper, by = 0.1)
  
  # Figure 1a: summary means.
  p1 <- ggplot(plot_summary, aes(x = mean, y = metric_label)) +
    geom_errorbarh(
      aes(xmin = ci95_low, xmax = ci95_high),
      height = 0.22,
      linewidth = 0.75,
      colour = "grey35",
      alpha = 0.95,
      na.rm = TRUE
    ) +
    geom_point(
      aes(fill = niche_type),
      shape = 21,
      size = 2.8,
      stroke = 0.35,
      colour = "grey12",
      na.rm = TRUE
    ) +
    facet_grid(
      niche_type + eval_set ~ workflow_label,
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    scale_fill_manual(values = niche_cols, guide = "none") +
    scale_x_continuous(
      limits = c(x_lower_fig1a, x_upper),
      breaks = x_breaks_fig1a,
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = "Figure 1a. One-hot binary RF performance across vegetation zones",
      subtitle = paste0(
        "Points are means across modeled vegetation zones; horizontal lines show mean +/- 1.96 SE as a descriptive across-zone interval. ",
        "OOB/Train metrics come from RF training summaries, and Test metrics come from independent balanced test-set assessment. ",
        "The x-axis starts at ", sprintf("%.2f", x_lower_fig1a),
        " so no point or interval is hidden."
      ),
      x = "Metric value",
      y = NULL,
      caption = "Modeled zones: 1-7, 9-50, 52-55."
    ) +
    theme_fig1(base_size = 10.2) +
    theme(
      axis.text.x = element_text(size = 8.5),
      axis.text.y = element_text(size = 8.6),
      strip.text.x = element_text(size = 9.6),
      strip.text.y = element_text(size = 9.1),
      axis.ticks.y = element_blank()
    )
  
  save_gg(
    p1,
    file.path(fig_dir, "Figure1a_onehot_binary_RF_summary.png"),
    width = 13.5,
    height = 9.6
  )
  
  save_gg(
    p1,
    file.path(fig_dir, "Figure1_onehot_binary_RF_OOB_train_test_performance_dotrange.png"),
    width = 13.5,
    height = 9.6
  )
  
  # Figure 1b: zone-level values.
  p2 <- ggplot(plot_zone, aes(x = value, y = metric_label)) +
    geom_point(
      aes(colour = zoneID_chr),
      position = position_jitter(width = 0, height = 0.14),
      size = 1.75,
      alpha = 0.75,
      stroke = 0,
      show.legend = FALSE
    ) +
    geom_errorbarh(
      data = plot_summary,
      aes(xmin = ci95_low, xmax = ci95_high, y = metric_label),
      inherit.aes = FALSE,
      height = 0.22,
      linewidth = 0.72,
      colour = "grey15",
      alpha = 0.95
    ) +
    geom_point(
      data = plot_summary,
      aes(x = mean, y = metric_label, fill = niche_type),
      inherit.aes = FALSE,
      shape = 23,
      size = 2.8,
      stroke = 0.35,
      colour = "grey10",
      show.legend = FALSE
    ) +
    facet_grid(
      niche_type + eval_set ~ workflow_label,
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    scale_colour_manual(values = zone_cols, guide = "none") +
    scale_fill_manual(values = niche_cols, guide = "none") +
    scale_x_continuous(
      limits = c(x_lower_fig1b, x_upper),
      breaks = x_breaks_fig1b,
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = "Figure 1b. Zone-level variation in one-hot binary RF performance",
      subtitle = "Small points are individual vegetation zones. Diamonds and horizontal lines show the corresponding across-zone mean and mean +/- 1.96 SE as a descriptive interval. The x-axis is zoomed without hiding any value.",
      x = "Metric value",
      y = NULL,
      caption = "Each point represents one modeled vegetation zone."
    ) +
    theme_fig1(base_size = 10.0) +
    theme(
      axis.text.x = element_text(size = 8.3),
      axis.text.y = element_text(size = 8.4),
      strip.text.x = element_text(size = 9.5),
      strip.text.y = element_text(size = 9.0),
      axis.ticks.y = element_blank()
    )
  
  save_gg(
    p2,
    file.path(fig_dir, "Figure1b_onehot_binary_RF_zone_level_metrics.png"),
    width = 13.5,
    height = 10.4
  )
})

# 7. Figure 2 ==================================================================

run_step("Figure 2 | Predicted reference maps for all workflows and original reference map", {
  if (!file.exists(reference_file)) stop("Missing reference raster: ", reference_file)
  
  # Remove old Figure 2 files so stale maps cannot survive a failed rerun.
  old_figure2 <- c(
    list.files(
      fig_dir,
      pattern = "^Figure2[a-d]_(original_vs_predicted|predicted)_reference_map_.*\\.png$",
      full.names = TRUE
    ),
    file.path(fig_dir, "Figure2_reference_vegetation_map.png")
  )
  remove_existing_outputs(old_figure2)
  
  figure_letters <- letters[seq_along(map_method_order)]
  
  # Figure 2a-d: one predicted reference map per workflow.
  predicted_files <- setNames(
    vapply(map_method_order, function(m) find_assigned_map(m, "normal"), character(1)),
    map_method_order
  )
  missing_methods <- names(predicted_files)[
    is.na(predicted_files) | !file.exists(predicted_files)
  ]
  if (length(missing_methods) > 0) {
    stop(
      "Figure 2 requires all four predicted reference maps. Missing: ",
      paste(missing_methods, collapse = ", ")
    )
  }
  
  for (i in seq_along(map_method_order)) {
    m <- map_method_order[i]
    pred_file <- predicted_files[[m]]
    
    plot_zone_panel_gg(
      files = list(pred_file),
      titles = map_method_labels[m],
      outfile = file.path(
        fig_dir,
        paste0("Figure2", figure_letters[i], "_predicted_reference_map_", m, ".png")
      ),
      categorical = TRUE,
      show_legend = FALSE,
      ncol = 1,
      width = 6.2,
      height = 4.8,
      plot_title = paste0("Predicted reference map: ", map_method_labels[m])
    )
  }
  
  # Original vegetation map: generated once as a separate figure.
  plot_zone_panel_gg(
    files = list(reference_file),
    titles = "Original vegetation map",
    outfile = file.path(fig_dir, "Figure2_reference_vegetation_map.png"),
    categorical = TRUE,
    show_legend = FALSE,
    ncol = 1,
    width = 6.2,
    height = 4.8,
    allowed_values = reference_zoneID,
    plot_title = "Original vegetation map"
  )
})

# 8. Table 2 | direct raster reconstruction metrics ===========================

run_step("Table 2 | Reference-map reconstruction performance for four binary workflows", {
  out_file <- file.path(tab_dir, "Table2_reference_map_accuracy_binary_workflows.csv")
  remove_existing_outputs(out_file)
  
  overall <- read_binary_map_overall()
  zone <- read_binary_map_zone()
  
  zone_summary <- zone[
    ,
    .(n_zones = uniqueN(zoneID)),
    by = method_key
  ]
  
  tab2 <- merge(
    overall,
    zone_summary,
    by = "method_key",
    all.x = TRUE,
    sort = FALSE
  )
  
  tab2[, `:=`(
    method_label = map_method_labels[method_key],
    individual_coverage_percent = 100 * individual_coverage,
    coverage_percent = 100 * coverage,
    exact_zone_agreement_percent = 100 * exact_zone_accuracy,
    broad_category_agreement_percent = 100 * broad_category_accuracy,
    overall_source_file = "Recalculated directly from four binary assigned-zone rasters",
    zone_source_file = "Recalculated directly from four binary assigned-zone rasters"
  )]
  
  drop_cols <- intersect(c("method_raw", "source_file"), names(tab2))
  if (length(drop_cols) > 0) tab2[, (drop_cols) := NULL]
  
  tab2[, method_order_tmp := match(method_key, map_method_order)]
  setorder(tab2, method_order_tmp)
  tab2[, method_order_tmp := NULL]
  setcolorder(
    tab2,
    c(
      "method_key", "method_label",
      "individual_coverage", "individual_coverage_percent",
      "coverage", "coverage_percent",
      "exact_zone_accuracy", "exact_zone_agreement_percent",
      "broad_category_accuracy", "broad_category_agreement_percent",
      "n_zones",
      "macro_balanced_accuracy", "macro_recall", "macro_specificity",
      "macro_precision", "macro_f1", "macro_tss",
      intersect(
        c("valid_original_pixels", "compared_pixels", "missing_predictions"),
        names(tab2)
      ),
      "overall_source_file", "zone_source_file"
    )
  )
  
  fwrite(tab2, out_file)
  cat0("[SAVED] ", out_file)
})

# 9. Figure 3 ==================================================================

run_step("Figure 3 | Climate and soil binary RF zone-level metrics for all workflows", {
  remove_existing_outputs(c(
    file.path(fig_dir, "Figure3a_climate_binary_RF_zone_metrics_all_workflows.png"),
    file.path(fig_dir, "Figure3b_soil_binary_RF_zone_metrics_all_workflows.png")
  ))
  
  long_file <- file.path(tab_dir, "Table1_binary_RF_zone_level_metrics_long.csv")
  
  if (file.exists(long_file)) {
    long <- fread(long_file)
  } else {
    tbl <- make_binary_rf_table1()
    long <- tbl$long
    fwrite(long, long_file)
  }
  
  metric_priority <- c(
    "Test F1", "Test precision", "Test recall",
    "F1", "Precision", "Recall",
    "OOB F1", "OOB precision", "OOB recall"
  )
  
  metrics_keep <- metric_priority[metric_priority %in% unique(long$metric)]
  if (length(metrics_keep) == 0) stop("No F1/precision/recall metrics found for Figure 3.")
  
  plot_dt <- long[metric %in% metrics_keep]
  plot_dt[, workflow := factor(workflow, levels = workflow_order)]
  plot_dt[, workflow_label := factor(
    workflow_labels[as.character(workflow)],
    levels = workflow_labels[workflow_order]
  )]
  plot_dt[, metric := factor(metric, levels = metrics_keep)]
  plot_dt[, zoneID_chr := as.character(zoneID)]
  plot_dt[, zoneID_fac := factor(zoneID_chr, levels = as.character(rev(model_zoneID)))]
  
  cols <- zone_color_vector(model_zoneID)
  
  for (nt in c("climate", "soil")) {
    dt_sub <- plot_dt[niche_type == nt]
    
    if (nrow(dt_sub) == 0) {
      stop("Figure 3 has no records for niche type: ", nt)
    }
    
    title_text <- ifelse(
      nt == "climate",
      "Climate binary RF zone-level performance",
      "Soil binary RF zone-level performance"
    )
    
    p <- ggplot(dt_sub, aes(x = value, y = zoneID_fac, colour = zoneID_chr)) +
      geom_point(size = 1.5, alpha = 0.75) +
      facet_grid(metric ~ workflow_label, scales = "free_x") +
      scale_colour_manual(values = cols, guide = "none") +
      labs(
        title = title_text,
        subtitle = "Each point is one vegetation zone; colors follow the vegetation-zone palette.",
        x = "Metric value",
        y = "Vegetation zone"
      ) +
      theme_ms(base_size = 9.5) +
      theme(
        axis.text.y = element_text(size = 6.2),
        axis.text.x = element_text(size = 7),
        panel.grid.major.y = element_blank()
      )
    
    save_gg(
      p,
      file.path(
        fig_dir,
        paste0("Figure3", ifelse(nt == "climate", "a", "b"), "_", nt, "_binary_RF_zone_metrics_all_workflows.png")
      ),
      width = 12.5,
      height = 8.5
    )
  }
})

# 10. Figure 4 =================================================================

run_step("Figure 4 | Major ecotype confusion flows", {
  out_file <- file.path(fig_dir, "Figure4_major_ecotype_confusion_flows.png")
  out_csv <- file.path(tab_dir, "Figure4_major_ecotype_confusion_flows.csv")
  remove_existing_outputs(c(out_file, out_csv))
  
  flow <- copy(get_binary_reference_comparison()$confusion)[
    method_key == preferred_method &
      ori %in% model_zoneID & pred %in% model_zoneID &
      ori != pred & is.finite(n) & n > 0,
    .(count = sum(n)),
    by = .(from = ori, to = pred)
  ]
  setorder(flow, -count, from, to)
  if (nrow(flow) == 0) stop("No off-diagonal confusion flows found.")
  
  top_n <- min(20L, nrow(flow))
  flow_top <- flow[seq_len(top_n)]
  flow_top[, `:=`(
    method_key = preferred_method,
    comparison_scope = "Common valid mask of the four binary reference maps",
    from_chr = as.character(from),
    to_chr = as.character(to),
    flow_label = paste0("Zone ", from, " -> Zone ", to)
  )]
  flow_top[, flow_label := factor(flow_label, levels = rev(flow_label))]
  fwrite(flow_top, out_csv)
  
  cols <- zone_color_vector(unique(flow_top$from))
  
  if (requireNamespace("ggalluvial", quietly = TRUE)) {
    p <- ggplot(flow_top, aes(y = count, axis1 = from_chr, axis2 = to_chr)) +
      ggalluvial::geom_alluvium(
        aes(fill = from_chr),
        width = 1 / 12,
        alpha = 0.78,
        show.legend = FALSE
      ) +
      ggalluvial::geom_stratum(width = 1 / 8, fill = "grey94", colour = "grey40") +
      ggplot2::geom_text(
        stat = "stratum",
        aes(label = after_stat(stratum)),
        size = 2.6
      ) +
      scale_fill_manual(values = cols) +
      scale_x_discrete(limits = c("Original", "Predicted"), expand = c(0.08, 0.08)) +
      labs(
        title = paste0("Major reference-map confusion flows: ", map_method_labels[preferred_method]),
        subtitle = "Top off-diagonal transitions on the common valid mask of the four binary workflows.",
        x = NULL,
        y = "Pixel count"
      ) +
      theme_ms()
  } else {
    p <- ggplot(flow_top, aes(x = count, y = flow_label, fill = from_chr)) +
      geom_col(width = 0.72, alpha = 0.80) +
      scale_fill_manual(values = cols, guide = "none") +
      labs(
        title = paste0("Major reference-map confusion flows: ", map_method_labels[preferred_method]),
        subtitle = "Top off-diagonal transitions on the common valid mask of the four binary workflows.",
        x = "Pixel count",
        y = NULL
      ) +
      theme_ms(base_size = 10) +
      theme(panel.grid.major.y = element_blank())
  }
  
  save_gg(p, out_file, width = 8.8, height = 6.6)
})

# 11. Figure 5 =================================================================

run_step("Figure 5 | Ecosystem niche maps for available workflows and periods", {
  figure5_outputs <- c(
    list.files(
      fig_dir,
      pattern = "^Figure5[a-d]_ecosystem_niche_maps_.*_2x4\\.png$",
      full.names = TRUE
    ),
    file.path(tab_dir, "Figure5_ecosystem_niche_map_index_all_models_periods.csv"),
    file.path(tab_dir, "Figure5_input_availability_audit.csv"),
    file.path(tab_dir, "Figure5_method_availability_summary.csv")
  )
  remove_existing_outputs(figure5_outputs)
  
  if (!file.exists(reference_file)) {
    stop("Missing reference raster: ", reference_file)
  }
  
  panel_titles <- c(
    "Original",
    "SSP245\n2011-2040",
    "SSP245\n2041-2070",
    "SSP245\n2071-2100",
    "Predicted\nreference",
    "SSP585\n2011-2040",
    "SSP585\n2041-2070",
    "SSP585\n2071-2100"
  )
  
  required_map_scenarios <- c("normal", future_order)
  
  preflight <- CJ(
    method = map_method_order,
    scenario = required_map_scenarios,
    unique = TRUE
  )
  
  preflight[, file := mapply(
    find_assigned_map,
    method,
    scenario,
    USE.NAMES = FALSE
  )]
  
  preflight[, exists := !is.na(file) & file.exists(file)]
  preflight[, figure_letter := letters[match(method, map_method_order)]]
  preflight[, method_label := unname(map_method_labels[method])]
  
  availability <- preflight[
    ,
    .(
      n_required = .N,
      n_available = sum(exists),
      complete_workflow = all(exists),
      missing_scenarios = paste(scenario[!exists], collapse = "; ")
    ),
    by = .(method, method_label, figure_letter)
  ]
  
  fwrite(
    preflight,
    file.path(tab_dir, "Figure5_input_availability_audit.csv")
  )
  fwrite(
    availability,
    file.path(tab_dir, "Figure5_method_availability_summary.csv")
  )
  
  complete_methods <- availability[
    complete_workflow == TRUE,
    method
  ]
  
  incomplete_methods <- availability[
    complete_workflow == FALSE
  ]
  
  if (nrow(incomplete_methods) > 0L) {
    cat0(
      "[FIGURE 5 NONFATAL MISSING INPUTS] ",
      paste(
        paste0(
          incomplete_methods$method,
          ": ",
          incomplete_methods$missing_scenarios
        ),
        collapse = " | "
      )
    )
    cat0(
      "[FIGURE 5 NOTE] Missing future assigned maps cannot be created by ",
      "visualization code. Complete workflows will still be plotted, and the ",
      "missing-input audit is saved."
    )
  }
  
  if (length(complete_methods) == 0L) {
    stop(
      "Figure 5 has no workflow with all six future assigned maps plus the ",
      "predicted reference map. See Figure5_input_availability_audit.csv."
    )
  }
  
  index_list <- list()
  
  for (m in complete_methods) {
    method_position <- match(m, map_method_order)
    figure_letter <- letters[method_position]
    
    panel_files <- c(
      reference_file,
      find_assigned_map(m, "2011-2040SSP245"),
      find_assigned_map(m, "2041-2070SSP245"),
      find_assigned_map(m, "2071-2100SSP245"),
      find_assigned_map(m, "normal"),
      find_assigned_map(m, "2011-2040SSP585"),
      find_assigned_map(m, "2041-2070SSP585"),
      find_assigned_map(m, "2071-2100SSP585")
    )
    
    if (any(is.na(panel_files)) || any(!file.exists(panel_files))) {
      stop(
        "Figure 5 preflight inconsistency for workflow ", m,
        ". A map disappeared after the availability audit."
      )
    }
    
    validate_common_prediction_mask(
      panel_files[-1],
      paste0("Figure 5 assigned maps for ", m)
    )
    
    output_file <- file.path(
      fig_dir,
      paste0(
        "Figure5",
        figure_letter,
        "_ecosystem_niche_maps_",
        m,
        "_2x4.png"
      )
    )
    
    plot_zone_panel_gg(
      files = as.list(panel_files),
      titles = panel_titles,
      outfile = output_file,
      categorical = TRUE,
      population_mode = FALSE,
      show_legend = FALSE,
      ncol = 4,
      width = 13.2,
      height = 6.8,
      draw_as_points = FALSE,
      point_alpha = 1,
      allowed_values = c(reference_zoneID, novel_value),
      plot_title = paste0(
        "Ecosystem niche projection: ",
        map_method_labels[m]
      )
    )
    
    index_list[[length(index_list) + 1L]] <- data.table(
      figure = paste0("Figure5", figure_letter),
      method = m,
      method_label = unname(map_method_labels[m]),
      panel = panel_titles,
      file = panel_files,
      output_file = output_file,
      output_generated = TRUE
    )
  }
  
  idx <- rbindlist(index_list, fill = TRUE)
  
  fwrite(
    idx,
    file.path(
      tab_dir,
      "Figure5_ecosystem_niche_map_index_all_models_periods.csv"
    )
  )
  
  cat0(
    "[FIGURE 5 GENERATED WORKFLOWS] ",
    paste(complete_methods, collapse = ", ")
  )
})



# Figure 6b helpers =============================================================

zone_short_label <- function(z) {
  z <- as.integer(z)
  ifelse(z == novel_value, "Novel", as.character(z))
}

build_transition_paths <- function(files, stage_names, family_label) {
  files <- as.character(files)
  if (length(files) != length(stage_names)) {
    stop("files and stage_names must have the same length.")
  }
  if (any(is.na(files)) || any(!file.exists(files))) {
    stop("One or more transition rasters are missing for ", family_label)
  }
  
  r_list <- lapply(files, function(f) rast(f)[[1]])
  for (i in 2:length(r_list)) {
    ok <- try(compareGeom(r_list[[1]], r_list[[i]], stopOnError = FALSE), silent = TRUE)
    if (inherits(ok, "try-error") || !isTRUE(ok)) {
      stop("Transition rasters are not geometrically aligned for ", family_label)
    }
  }
  
  rs <- r_list[[1]]
  if (length(r_list) > 1) {
    for (i in 2:length(r_list)) rs <- c(rs, r_list[[i]])
  }
  area_r <- cellSize(r_list[[1]], unit = "km")
  rs <- c(rs, area_r)
  names(rs) <- c(paste0("stage", seq_along(stage_names)), "area_km2")
  
  dt <- as.data.table(as.data.frame(rs, na.rm = FALSE))
  if (nrow(dt) == 0) return(data.table())
  
  stage_cols <- paste0("stage", seq_along(stage_names))
  for (sc in stage_cols) dt[, (sc) := as.integer(round(get(sc)))]
  
  allowed_vals <- c(model_zoneID, novel_value)
  keep <- complete.cases(dt[, ..stage_cols]) & !is.na(dt$area_km2) & is.finite(dt$area_km2) & dt$area_km2 > 0
  for (sc in stage_cols) keep <- keep & dt[[sc]] %in% allowed_vals
  dt <- dt[keep]
  if (nrow(dt) == 0) return(data.table())
  
  out <- dt[, .(area_km2 = sum(area_km2, na.rm = TRUE)), by = stage_cols]
  if (nrow(out) == 0) return(data.table())
  
  out[, family := family_label]
  out[, family := factor(family, levels = c("SSP245", "SSP585"))]
  out[, area_million_km2 := area_km2 / 1e6]
  
  for (i in seq_along(stage_names)) {
    sc <- paste0("stage", i)
    ax <- paste0("axis", i)
    out[, (ax) := zone_short_label(get(sc))]
  }
  
  out[]
}

build_pairwise_transitions <- function(path_dt, stage_names) {
  if (nrow(path_dt) == 0) return(data.table())
  
  out_list <- vector(
    "list",
    max(0L, length(stage_names) - 1L)
  )
  
  for (i in seq_len(length(stage_names) - 1L)) {
    from_col <- paste0("stage", i)
    to_col <- paste0("stage", i + 1L)
    interval_label <- paste0(
      stage_names[i],
      " -> ",
      stage_names[i + 1L]
    )
    
    # Create full-length columns before grouping. A scalar interval expression
    # inside data.table's by-list is not recycled and caused the previous
    # lengths [n, 1, n, n] error.
    tmp <- data.table(
      family = path_dt[["family"]],
      interval = rep(interval_label, nrow(path_dt)),
      from_zone = as.integer(path_dt[[from_col]]),
      to_zone = as.integer(path_dt[[to_col]]),
      area_km2 = as.numeric(path_dt[["area_km2"]])
    )
    
    tmp <- tmp[
      is.finite(area_km2) &
        area_km2 > 0 &
        !is.na(from_zone) &
        !is.na(to_zone)
    ]
    
    tmp <- tmp[
      ,
      .(
        area_km2 = sum(area_km2, na.rm = TRUE)
      ),
      by = .(
        family,
        interval,
        from_zone,
        to_zone
      )
    ]
    
    out_list[[i]] <- tmp
  }
  
  out <- rbindlist(out_list, fill = TRUE)
  
  if (nrow(out) > 0L) {
    out[, interval_order := match(
      interval,
      paste0(
        stage_names[-length(stage_names)],
        " -> ",
        stage_names[-1L]
      )
    )]
    setorder(out, family, interval_order, from_zone, to_zone)
    out[, interval_order := NULL]
  }
  
  out[]
}

plot_transition_sankey <- function(path_dt, stage_names, outfile) {
  if (nrow(path_dt) == 0) stop("No transition paths are available for Figure 6b.")
  
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    # Dependency-safe fallback: preserve the transition information as a
    # pairwise stacked-flow summary rather than failing the entire script.
    pair_dt <- build_pairwise_transitions(path_dt, stage_names)
    if (nrow(pair_dt) == 0) {
      stop("No pairwise transition data are available for Figure 6b fallback.")
    }
    pair_dt[, from_zone_chr := factor(
      zone_short_label(from_zone),
      levels = zone_short_label(c(model_zoneID, novel_value))
    )]
    pair_dt[, area_million_km2 := area_km2 / 1e6]
    cols <- zone_color_vector(c(model_zoneID, novel_value))
    names(cols) <- zone_short_label(as.integer(names(cols)))
    
    p_fallback <- ggplot(
      pair_dt,
      aes(x = interval, y = area_million_km2, fill = from_zone_chr)
    ) +
      geom_col(width = 0.72) +
      facet_grid(family ~ ., scales = "free_x") +
      scale_fill_manual(values = cols, drop = FALSE) +
      labs(
        title = paste0(
          "Figure 6b. Pairwise ecotype-transition summaries: ",
          map_method_labels[preferred_method]
        ),
        subtitle = paste0(
          "Fallback used because package 'ggalluvial' is unavailable. ",
          "Bars preserve all pairwise transition areas."
        ),
        x = NULL,
        y = expression("Transition area (" * 10^6 * " km"^2 * ")"),
        fill = "Source zone"
      ) +
      theme_ms(base_size = 9.5) +
      theme(
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.text = element_text(size = 6),
        legend.key.size = grid::unit(0.28, "cm")
      ) +
      guides(fill = guide_legend(ncol = 3, byrow = TRUE))
    
    save_gg(p_fallback, outfile, width = 12.8, height = 8.8)
    return(invisible(NULL))
  }
  
  plot_dt <- copy(path_dt)
  plot_dt[, alluvium_id := .I]
  plot_dt[, reference_zone_chr := as.character(stage1)]
  
  cols <- zone_color_vector(unique(plot_dt$stage1))
  
  p <- ggplot(
    plot_dt,
    aes(
      y = area_million_km2,
      axis1 = axis1,
      axis2 = axis2,
      axis3 = axis3,
      axis4 = axis4
    )
  ) +
    ggalluvial::geom_alluvium(
      aes(fill = reference_zone_chr),
      width = 0.14,
      alpha = 0.62,
      knot.pos = 0.35,
      decreasing = FALSE,
      show.legend = FALSE
    ) +
    ggalluvial::geom_stratum(
      width = 0.16,
      fill = "grey95",
      colour = "grey45",
      linewidth = 0.28
    ) +
    ggplot2::geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum)),
      size = 2.0,
      colour = "grey10",
      check_overlap = TRUE
    ) +
    facet_grid(family ~ ., switch = "y") +
    scale_fill_manual(values = cols) +
    scale_x_discrete(
      limits = stage_names,
      expand = c(0.06, 0.02)
    ) +
    labs(
      title = paste0("Figure 6b. Ecotype transition trajectories: ", map_method_labels[preferred_method]),
      subtitle = "Rows show SSP245 and SSP585. Each flow tracks area from the predicted reference map through three future periods.",
      x = NULL,
      y = expression("Transition area (" * 10^6 * " km"^2 * ")")
    ) +
    theme_ms(base_size = 10) +
    theme(
      strip.text.y = element_text(angle = 0, face = "bold"),
      axis.text.x = element_text(face = "bold"),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.spacing.y = grid::unit(0.9, "lines")
    )
  
  save_gg(p, outfile, width = 12.8, height = 8.8)
}

# 12. Figure 6 =================================================================

run_step("Figure 6 | Area change and zone-transition trajectories", {
  remove_existing_outputs(c(
    file.path(fig_dir, "Figure6a_all_zone_area_change.png"),
    file.path(fig_dir, "Figure6_all_zone_area_change.png"),
    file.path(fig_dir, "Figure6b_zone_transition_sankey.png"),
    file.path(tab_dir, "Figure6_all_zone_area_by_scenario.csv"),
    file.path(tab_dir, "Figure6b_zone_transition_paths.csv"),
    file.path(tab_dir, "Figure6b_zone_transition_pairwise_area.csv")
  ))
  
  map_files <- sapply(
    scenario_levels,
    function(s) find_assigned_map(preferred_method, s),
    USE.NAMES = TRUE
  )
  
  missing_scenarios <- names(map_files)[is.na(map_files) | !file.exists(map_files)]
  if (length(missing_scenarios) > 0) {
    stop(
      "Figure 6 requires the predicted reference map and all six future maps. Missing: ",
      paste(missing_scenarios, collapse = ", ")
    )
  }
  validate_common_prediction_mask(
    map_files,
    paste0("Figure 6 assigned maps for ", preferred_method)
  )
  
  area_list <- list()
  
  for (s in names(map_files)) {
    cat0("  Calculating area: ", preferred_method, " | ", s)
    dt <- area_by_zone(map_files[[s]])
    dt[, scenario := s]
    dt[, scenario_label := scenario_label(s)]
    dt[, scenario_order := match(s, scenario_levels)]
    area_list[[length(area_list) + 1L]] <- dt
  }
  
  area_dt <- rbindlist(area_list, fill = TRUE)
  
  # Store explicit zeroes for absent classes so the table has a complete
  # scenario-by-zone grid and cannot silently omit extirpated/novel classes.
  area_grid <- CJ(
    scenario = scenario_levels,
    zoneID = c(model_zoneID, novel_value),
    unique = TRUE
  )
  area_dt <- merge(
    area_grid,
    area_dt[, .(scenario, zoneID, area_km2)],
    by = c("scenario", "zoneID"),
    all.x = TRUE,
    sort = FALSE
  )
  area_dt[is.na(area_km2), area_km2 := 0]
  area_dt[, `:=`(
    scenario_label = scenario_label(scenario),
    scenario_order = match(scenario, scenario_levels)
  )]
  setorder(area_dt, scenario_order, zoneID)
  
  area_totals <- area_dt[, .(total_area_km2 = sum(area_km2)), by = scenario]
  area_tolerance <- max(1e-6, max(area_totals$total_area_km2) * 1e-10)
  if (max(area_totals$total_area_km2) - min(area_totals$total_area_km2) > area_tolerance) {
    stop(
      "Figure 6 area totals differ across scenarios despite the required common mask: ",
      paste(capture.output(print(area_totals)), collapse = " ")
    )
  }
  
  area_file <- file.path(tab_dir, "Figure6_all_zone_area_by_scenario.csv")
  fwrite(area_dt, area_file)
  cat0("[SAVED] ", area_file)
  
  area_plot <- copy(area_dt)
  area_plot[, area_million_km2 := area_km2 / 1e6]
  area_plot[, scenario_label := factor(scenario_label, levels = scenario_label(names(map_files)))]
  
  zone_levels <- as.character(model_zoneID)
  if (novel_value %in% area_plot$zoneID) zone_levels <- c(zone_levels, "Novel")
  
  area_plot[, zone_label := as.character(zoneID)]
  area_plot[zoneID == novel_value, zone_label := "Novel"]
  area_plot[, zone_label := factor(zone_label, levels = zone_levels)]
  
  zone_cols <- zone_color_vector(model_zoneID)
  fill_cols <- zone_cols
  if (novel_value %in% area_plot$zoneID) fill_cols <- c(fill_cols, Novel = "#333333")
  
  p <- ggplot(area_plot, aes(x = scenario_label, y = area_million_km2, fill = zone_label)) +
    geom_col(width = 0.58) +
    scale_fill_manual(values = fill_cols, drop = FALSE) +
    scale_x_discrete(expand = expansion(add = 0.08)) +
    labs(
      title = paste0("Figure 6a. Projected ecotype area change: ", map_method_labels[preferred_method]),
      subtitle = "All modeled zones are shown separately; colors follow the vegetation-zone palette.",
      x = NULL,
      y = expression("Area (" * 10^6 * " km"^2 * ")"),
      fill = "Zone"
    ) +
    theme_ms() +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "right",
      legend.text = element_text(size = 6.2),
      legend.key.size = grid::unit(0.30, "cm")
    ) +
    guides(fill = guide_legend(ncol = 3, byrow = TRUE))
  
  save_gg(p, file.path(fig_dir, "Figure6a_all_zone_area_change.png"), width = 10.2, height = 6.4)
  save_gg(p, file.path(fig_dir, "Figure6_all_zone_area_change.png"), width = 10.2, height = 6.4)
  
  # Figure 6b: transition trajectories across stages.
  stage_names <- c("Predicted reference", "2011-2040", "2041-2070", "2071-2100")
  
  files_245 <- c(
    find_assigned_map(preferred_method, "normal"),
    find_assigned_map(preferred_method, "2011-2040SSP245"),
    find_assigned_map(preferred_method, "2041-2070SSP245"),
    find_assigned_map(preferred_method, "2071-2100SSP245")
  )
  
  files_585 <- c(
    find_assigned_map(preferred_method, "normal"),
    find_assigned_map(preferred_method, "2011-2040SSP585"),
    find_assigned_map(preferred_method, "2041-2070SSP585"),
    find_assigned_map(preferred_method, "2071-2100SSP585")
  )
  
  path_245 <- build_transition_paths(files_245, stage_names, "SSP245")
  path_585 <- build_transition_paths(files_585, stage_names, "SSP585")
  path_dt <- rbindlist(list(path_245, path_585), fill = TRUE)
  
  if (nrow(path_dt) == 0) {
    stop("No transition trajectories could be derived for Figure 6b.")
  }
  
  path_totals <- path_dt[, .(total_path_area_km2 = sum(area_km2)), by = family]
  expected_path_area <- area_totals[scenario == "normal", total_area_km2]
  if (length(expected_path_area) != 1L || any(
    abs(path_totals$total_path_area_km2 - expected_path_area) > area_tolerance
  )) {
    stop(
      "Figure 6b path areas do not equal the common mapped area. ",
      paste(capture.output(print(path_totals)), collapse = " ")
    )
  }
  
  fwrite(path_dt, file.path(tab_dir, "Figure6b_zone_transition_paths.csv"))
  
  pair_dt <- build_pairwise_transitions(path_dt, stage_names)
  fwrite(pair_dt, file.path(tab_dir, "Figure6b_zone_transition_pairwise_area.csv"))
  
  plot_transition_sankey(
    path_dt = path_dt,
    stage_names = stage_names,
    outfile = file.path(fig_dir, "Figure6b_zone_transition_sankey.png")
  )
})

# 13. Table 3 ==================================================================

shannon_H <- function(abundance) {
  abundance <- as.numeric(abundance)
  abundance <- abundance[is.finite(abundance) & abundance > 0]
  if (length(abundance) == 0) return(NA_real_)
  p <- abundance / sum(abundance)
  -sum(p * log(p))
}

run_step("Table 3 | Reference population abundance and diversity by species", {
  table3_outputs <- c(
    file.path(tab_dir, "Table3_ten_species_populations_ShannonH.csv"),
    file.path(tab_dir, "Table3_all_species_population_diversity.csv"),
    file.path(tab_dir, "Table3_species_population_abundance_long.csv")
  )
  remove_existing_outputs(table3_outputs)
  
  pop <- read_population_lookup()
  
  if (anyNA(pop$reference_abundance)) {
    stop(
      "Population lookup contains missing reference_abundance values. ",
      "Table 3 and Figure 7 require the code 6 abundance column."
    )
  }
  
  pop <- pop[is.finite(reference_abundance) & reference_abundance >= 0]
  if (nrow(pop) == 0) stop("No finite reference population abundances are available.")
  
  ref_summary <- pop[
    reference_abundance > 0,
    .(
      n_reference_populations = .N,
      reference_total_abundance = sum(reference_abundance),
      reference_shannon_H = shannon_H(reference_abundance)
    ),
    by = species
  ]
  
  eligible_summary <- pop[
    projected == TRUE &
      source_zone %in% model_zoneID &
      reference_abundance > 0,
    .(
      n_projection_eligible_populations = .N,
      projection_eligible_total_abundance = sum(reference_abundance),
      projection_eligible_shannon_H = shannon_H(reference_abundance)
    ),
    by = species
  ]
  
  tab3 <- merge(
    ref_summary,
    eligible_summary,
    by = "species",
    all = TRUE,
    sort = FALSE
  )
  
  count_cols <- c(
    "n_reference_populations",
    "n_projection_eligible_populations",
    "reference_total_abundance",
    "projection_eligible_total_abundance"
  )
  for (cc in count_cols) set(tab3, which(is.na(tab3[[cc]])), cc, 0)
  
  tab3[, `:=`(
    n_excluded_unmodeled_zone_populations =
      n_reference_populations - n_projection_eligible_populations,
    # Compatibility columns: these explicitly refer to the full reference data.
    n_populations = n_reference_populations,
    total_abundance = reference_total_abundance,
    shannon_H = reference_shannon_H,
    source_file = unique(pop$source_file)[1]
  )]
  
  setorder(tab3, -n_reference_populations, -reference_shannon_H, species)
  tab3_top <- tab3[seq_len(min(10L, .N))]
  
  out_file <- file.path(tab_dir, "Table3_ten_species_populations_ShannonH.csv")
  fwrite(tab3_top, out_file)
  
  all_file <- file.path(tab_dir, "Table3_all_species_population_diversity.csv")
  fwrite(tab3, all_file)
  
  pop_long <- pop[
    ,
    .(
      species,
      zoneID = source_zone,
      abundance = reference_abundance,
      projection_eligible = projected,
      source_file
    )
  ]
  setorder(pop_long, species, zoneID)
  
  pop_file <- file.path(tab_dir, "Table3_species_population_abundance_long.csv")
  fwrite(pop_long, pop_file)
  
  cat0("[SAVED] ", out_file)
  cat0("[SAVED] ", all_file)
  cat0("[SAVED] ", pop_file)
})

# 14. Figure 7 =================================================================

run_step("Figure 7 | Species x ecotype reference population-abundance heatmap", {
  figure7_file <- file.path(
    fig_dir,
    "Figure7_species_ecotype_population_abundance_heatmap.png"
  )
  remove_existing_outputs(figure7_file)
  
  pop_file <- file.path(tab_dir, "Table3_species_population_abundance_long.csv")
  tab3_file <- file.path(tab_dir, "Table3_ten_species_populations_ShannonH.csv")
  
  if (!file.exists(pop_file) || !file.exists(tab3_file)) {
    stop("Corrected Table 3 outputs were not found.")
  }
  
  pop <- fread(pop_file)
  tab3 <- fread(tab3_file)
  species_keep <- as.character(tab3$species)
  
  observed <- pop[species %in% species_keep]
  if (nrow(observed) == 0) stop("No population-abundance records for selected species.")
  
  zones_keep <- sort(unique(observed[abundance > 0, zoneID]))
  if (length(zones_keep) == 0) stop("Selected species have no positive abundance.")
  
  heat <- CJ(
    species = species_keep,
    zoneID = zones_keep,
    unique = TRUE
  )
  heat <- merge(
    heat,
    observed[, .(species, zoneID, abundance)],
    by = c("species", "zoneID"),
    all.x = TRUE,
    sort = FALSE
  )
  heat[, `:=`(
    log10_abundance = fifelse(
      is.finite(abundance),
      log10(abundance + 1),
      NA_real_
    ),
    zoneID_chr = as.character(zoneID),
    zone_modeled = zoneID %in% model_zoneID
  )]
  heat[, species := factor(species, levels = rev(species_keep))]
  heat[, zoneID_fac := factor(zoneID_chr, levels = as.character(zones_keep))]
  
  zone_axis_labels <- setNames(
    ifelse(
      zones_keep %in% unmodeled_zoneID,
      paste0(zones_keep, "*"),
      as.character(zones_keep)
    ),
    as.character(zones_keep)
  )
  
  p <- ggplot(heat, aes(x = zoneID_fac, y = species, fill = log10_abundance)) +
    geom_tile(colour = "white", linewidth = 0.18) +
    scale_x_discrete(labels = zone_axis_labels) +
    scale_fill_gradient(
      low = "#F5F5F2",
      high = "#3F4A4A",
      name = "log10(reference\nabundance + 1)",
      na.value = "white"
    ) +
    labs(
      title = "Species-by-ecotype reference population abundance",
      subtitle = paste0(
        "Abundance is read directly from population_projection_lookup.csv. Blank cells indicate ",
        "no retained population (fewer than ", min_reference_population_abundance,
        " occupied cells or no recorded occurrence)."
      ),
      x = "Ecotype / vegetation zone",
      y = NULL,
      caption = "* Zones 8 and 51 occur in the reference data but are excluded from future ecosystem projections."
    ) +
    theme_ms() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.5),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank()
    )
  
  save_gg(
    p,
    figure7_file,
    width = 10.5,
    height = 5.8
  )
})

# 15. Figure 8 =================================================================

run_step("Figure 8 | Tree population dual-niche projection maps", {
  pop_tbl <- get_species_population_table()
  
  expected_pages <- ceiling(
    uniqueN(pop_tbl$species) /
      (species_page_ncol * species_page_nrow)
  )
  
  expected_files <- unlist(
    lapply(
      tree_plot_scenarios,
      function(s) {
        file.path(
          fig_dir,
          paste0(
            "Figure8_tree_population_dual_niche_projection_",
            safe_name(tree_plot_method),
            "_",
            safe_name(s),
            "_page",
            sprintf("%02d", seq_len(expected_pages)),
            ".png"
          )
        )
      }
    ),
    use.names = FALSE
  )
  
  source_table_file <- file.path(
    tab_dir,
    "Figure8_species_population_source_zones.csv"
  )
  
  if (
    resume_existing_outputs &&
    length(expected_files) > 0L &&
    all(file.exists(expected_files)) &&
    file.exists(source_table_file)
  ) {
    cat0(
      "[REUSE] Figure 8 already has all ",
      length(expected_files),
      " expected pages. Expensive raster plotting was skipped."
    )
  } else {
    validate_dual_scenario_grid(
      tree_plot_method,
      tree_plot_scenarios,
      sort(unique(pop_tbl$source_zone))
    )
    
    fwrite(pop_tbl, source_table_file)
    
    for (s in tree_plot_scenarios) {
      scenario_expected <- file.path(
        fig_dir,
        paste0(
          "Figure8_tree_population_dual_niche_projection_",
          safe_name(tree_plot_method),
          "_",
          safe_name(s),
          "_page",
          sprintf("%02d", seq_len(expected_pages)),
          ".png"
        )
      )
      
      if (
        resume_existing_outputs &&
        all(file.exists(scenario_expected))
      ) {
        cat0(
          "[REUSE] Figure 8 | ",
          s,
          " already complete."
        )
      } else {
        remove_existing_outputs(scenario_expected)
        
        plot_species_population_dual_pages(
          pop_tbl = pop_tbl,
          method = tree_plot_method,
          scenario = s,
          outfile_prefix = file.path(
            fig_dir,
            "Figure8_tree_population_dual_niche_projection"
          )
        )
      }
    }
    
    missing_after <- expected_files[!file.exists(expected_files)]
    
    if (length(missing_after) > 0L) {
      stop(
        "Figure 8 is missing expected pages after generation: ",
        paste(missing_after, collapse = "; ")
      )
    }
  }
})



# 16. Figure 9 =================================================================

run_step("Figure 9 | Species-level dual-niche projection maps", {
  pop_tbl <- get_species_population_table()
  
  expected_pages <- ceiling(
    uniqueN(pop_tbl$species) /
      (species_page_ncol * species_page_nrow)
  )
  
  expected_files <- unlist(
    lapply(
      tree_plot_scenarios,
      function(s) {
        file.path(
          fig_dir,
          paste0(
            "Figure9_species_level_dual_niche_projection_",
            safe_name(tree_plot_method),
            "_",
            safe_name(s),
            "_page",
            sprintf("%02d", seq_len(expected_pages)),
            ".png"
          )
        )
      }
    ),
    use.names = FALSE
  )
  
  source_table_file <- file.path(
    tab_dir,
    "Figure9_species_population_source_zones.csv"
  )
  
  if (
    resume_existing_outputs &&
    length(expected_files) > 0L &&
    all(file.exists(expected_files)) &&
    file.exists(source_table_file)
  ) {
    cat0(
      "[REUSE] Figure 9 already has all ",
      length(expected_files),
      " expected pages. Expensive raster plotting was skipped."
    )
  } else {
    validate_dual_scenario_grid(
      tree_plot_method,
      tree_plot_scenarios,
      sort(unique(pop_tbl$source_zone))
    )
    
    fwrite(pop_tbl, source_table_file)
    
    for (s in tree_plot_scenarios) {
      scenario_expected <- file.path(
        fig_dir,
        paste0(
          "Figure9_species_level_dual_niche_projection_",
          safe_name(tree_plot_method),
          "_",
          safe_name(s),
          "_page",
          sprintf("%02d", seq_len(expected_pages)),
          ".png"
        )
      )
      
      if (
        resume_existing_outputs &&
        all(file.exists(scenario_expected))
      ) {
        cat0(
          "[REUSE] Figure 9 | ",
          s,
          " already complete."
        )
      } else {
        remove_existing_outputs(scenario_expected)
        
        plot_species_level_dual_pages(
          pop_tbl = pop_tbl,
          method = tree_plot_method,
          scenario = s,
          outfile_prefix = file.path(
            fig_dir,
            "Figure9_species_level_dual_niche_projection"
          )
        )
      }
    }
    
    missing_after <- expected_files[!file.exists(expected_files)]
    
    if (length(missing_after) > 0L) {
      stop(
        "Figure 9 is missing expected pages after generation: ",
        paste(missing_after, collapse = "; ")
      )
    }
  }
})



# 17. Table 4 ==================================================================

run_step("Table 4 | Optimized binary vs multiclass reference-map robustness check", {
  table4_outputs <- c(
    file.path(tab_dir, "Table4_map_metric_source_consistency_audit.csv"),
    file.path(tab_dir, "Table4_reference_map_mask_scope_audit.csv"),
    file.path(tab_dir, "Table4_binary_vs_multiclass_robustness_check.csv"),
    file.path(
      tab_dir,
      "Table4_optimized_binary_vs_multiclass_reference_map_robustness.csv"
    )
  )
  remove_existing_outputs(table4_outputs)
  
  audit <- audit_map_metric_consistency()
  fwrite(
    audit,
    file.path(tab_dir, "Table4_map_metric_source_consistency_audit.csv")
  )
  fwrite(
    audit,
    file.path(tab_dir, "Table4_reference_map_mask_scope_audit.csv")
  )
  
  overall <- read_map_comparison_overall()
  robust <- overall[method_key %in% c("optimized_rf", "multiclass_rf")]
  
  metric_labels <- c(
    individual_coverage = "Individual map coverage",
    coverage = "Common comparison coverage",
    exact_zone_accuracy = "Exact-zone accuracy",
    broad_category_accuracy = "Broad-category accuracy",
    macro_balanced_accuracy = "Macro balanced accuracy",
    macro_recall = "Macro recall",
    macro_specificity = "Macro specificity",
    macro_precision = "Macro precision",
    macro_f1 = "Macro F1",
    macro_tss = "Macro TSS"
  )
  
  metrics_keep <- names(metric_labels)[
    names(metric_labels) %in% names(robust)
  ]
  metrics_keep <- metrics_keep[
    vapply(
      metrics_keep,
      function(x) any(is.finite(robust[[x]])),
      logical(1)
    )
  ]
  
  if (length(metrics_keep) == 0) {
    stop("No finite, same-scope robustness metrics were found.")
  }
  
  long <- melt(
    robust,
    id.vars = c("method_key", "method_label", "source_file"),
    measure.vars = metrics_keep,
    variable.name = "metric",
    value.name = "value"
  )
  long[, metric_label := unname(metric_labels[as.character(metric)])]
  
  tab4 <- dcast(
    long,
    metric + metric_label ~ method_key,
    value.var = "value"
  )
  
  require_columns(
    tab4,
    c("optimized_rf", "multiclass_rf"),
    "Corrected Table 4"
  )
  
  tab4[, `:=`(
    multiclass_minus_optimized_binary = multiclass_rf - optimized_rf,
    higher_value = fifelse(
      multiclass_rf > optimized_rf,
      "Multiclass RF",
      fifelse(
        multiclass_rf < optimized_rf,
        "Optimized binary RF",
        "Equal"
      )
    ),
    comparison_scope = "Reference-map reconstruction against the original vegetation raster",
    source_file = "Recalculated directly from original and assigned-zone rasters"
  )]
  
  tab4[, metric_order := match(metric, metrics_keep)]
  setorder(tab4, metric_order)
  tab4[, metric_order := NULL]
  
  out_file <- file.path(
    tab_dir,
    "Table4_binary_vs_multiclass_robustness_check.csv"
  )
  descriptive_file <- file.path(
    tab_dir,
    "Table4_optimized_binary_vs_multiclass_reference_map_robustness.csv"
  )
  
  fwrite(tab4, out_file)
  fwrite(tab4, descriptive_file)
  cat0("[SAVED] ", out_file)
  cat0("[SAVED] ", descriptive_file)
})

# 19. Figure 10b ================================================================

run_step("Figure 10b | Overall reference-map reconstruction metrics for all workflows", {
  figure10b_png <- file.path(
    fig_dir,
    "Figure10b_all_workflows_reference_map_overall_metrics.png"
  )
  figure10b_long_csv <- file.path(
    tab_dir,
    "Figure10b_all_workflows_reference_map_overall_metrics_long.csv"
  )
  figure10b_csv <- file.path(
    tab_dir,
    "Figure10b_all_workflows_reference_map_overall_metrics.csv"
  )
  
  # Remove the obsolete testing-set comparison. Binary climate/soil one-vs-rest
  # tests and the multiclass test are different tasks and test samples.
  obsolete_files <- c(
    file.path(
      fig_dir,
      "Figure10b_all_binary_workflows_vs_multiclass_testing_set_F1_area_scaled.png"
    ),
    file.path(
      tab_dir,
      "Figure10b_all_binary_workflows_vs_multiclass_testing_set_F1.csv"
    ),
    figure10b_png,
    figure10b_long_csv,
    figure10b_csv
  )
  obsolete_files <- obsolete_files[file.exists(obsolete_files)]
  if (length(obsolete_files) > 0) {
    unlink(obsolete_files)
    cat0("[REMOVED OBSOLETE] ", paste(basename(obsolete_files), collapse = ", "))
  }
  
  invisible(audit_map_metric_consistency())
  overall <- read_map_comparison_overall()
  
  metric_labels <- c(
    individual_coverage = "Individual map coverage",
    coverage = "Common comparison coverage",
    exact_zone_accuracy = "Exact-zone accuracy",
    broad_category_accuracy = "Broad-category accuracy",
    macro_balanced_accuracy = "Macro balanced accuracy",
    macro_recall = "Macro recall",
    macro_specificity = "Macro specificity",
    macro_precision = "Macro precision",
    macro_f1 = "Macro F1"
  )
  
  metrics_keep <- names(metric_labels)[
    names(metric_labels) %in% names(overall)
  ]
  metrics_keep <- metrics_keep[
    vapply(
      metrics_keep,
      function(x) any(is.finite(overall[[x]])),
      logical(1)
    )
  ]
  
  if (length(metrics_keep) == 0) {
    stop("No finite overall map-reconstruction metrics are available.")
  }
  
  plot_dt <- melt(
    overall,
    id.vars = c("method_key", "method_label", "source_file"),
    measure.vars = metrics_keep,
    variable.name = "metric",
    value.name = "value"
  )
  plot_dt <- plot_dt[is.finite(value)]
  plot_dt[, `:=`(
    metric_label = factor(
      unname(metric_labels[as.character(metric)]),
      levels = unname(metric_labels[metrics_keep])
    ),
    method_label = factor(
      method_label,
      levels = rev(all_map_method_labels[all_map_method_order])
    ),
    model_family = fifelse(
      method_key == "multiclass_rf",
      "Multiclass robustness model",
      "Binary dual-niche workflow"
    )
  )]
  
  p <- ggplot(
    plot_dt,
    aes(x = value, y = method_label, shape = model_family)
  ) +
    geom_point(size = 2.8, stroke = 0.7) +
    facet_wrap(~ metric_label, ncol = 2) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = expansion(mult = c(0.01, 0.03))
    ) +
    labs(
      title = "Overall reference-map reconstruction performance",
      subtitle = paste0(
        "All values come from the same map-overlay assessment against the original ",
        "vegetation raster. This is reconstruction agreement, not independent validation."
      ),
      x = "Metric value",
      y = NULL,
      shape = "Model family"
    ) +
    theme_ms(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 8.2),
      legend.position = "bottom"
    )
  
  save_gg(
    p,
    figure10b_png,
    width = 10.5,
    height = 9.0
  )
  
  fwrite(
    plot_dt,
    figure10b_long_csv
  )
  fwrite(
    overall,
    figure10b_csv
  )
})

# 20. Save step log =============================================================

cat0("\n============================================================")
cat0("SAVE VISUALIZATION STEP LOG")
cat0("============================================================")

step_log_dt <- rbindlist(step_log, fill = TRUE)
log_file <- file.path(vis_dir, "visualization_step_log.csv")
fwrite(step_log_dt, log_file)

cat0("[SAVED] ", log_file)

failed_steps <- step_log_dt[status != "done"]
if (nrow(failed_steps) > 0) {
  cat0("\nVISUALIZATION INCOMPLETE")
  print(failed_steps[, .(step, status, message)])
  stop(
    "One or more manuscript visualization steps failed after all independent ",
    "steps were attempted. No COMPLETE status is reported. See: ",
    log_file,
    call. = FALSE
  )
}

cat0("\nCOMPLETE")
cat0(
  "All runnable visualization steps completed and passed validation. ",
  "Review Figure5_input_availability_audit.csv for workflows whose future ",
  "assigned maps have not yet been generated."
)
cat0("Figure folder: ", fig_dir)
cat0("Table folder: ", tab_dir)
cat0("Step log: ", log_file)
