# 11.2 Visualization for final plain RF and plain multi-Forest models
# ==============================================================================
# Rebuild the original manuscript visualization structure using only:
#   rf_var : ordinary RF using mcRFop-selected climate variables
#   mf_var : plain multi-Forest using mcRFop-selected climate variables
#
# Reproduced outputs:
#   Figure 1a-b  model performance
#   Figure 2a-b  predicted reference maps + original reference map
#   Figure 3a-b  zone-level climate/soil metrics
#   Figure 4     major ecotype confusion flows
#   Figure 5a-b  ecosystem niche maps (2 rows x 4 columns)
#   Figure 6a-b  projected ecotype area + transition Sankey
#   Figure 7     species x ecotype reference-population abundance heatmap
#   Figure 8     population dual-niche maps (10 species/page; 2 x 5; alpha 0.30)
#   Figure 9     species dual-niche maps (10 species/page; 2 x 5; alpha 0.30)
#   Figure 10a-b optional multiclass robustness comparison
#
# Also reproduces the old normal-map chord PDFs for rf_var and mf_var.
#
# Scientific rules:
#   - soil gate was applied upstream at 0.2 in script 3.21;
#   - final suitability threshold is 0.4;
#   - future pixel is novel only when max dual suitability < 0.4;
#   - population/species display therefore uses dual suitability >= 0.4.
#
# Script 8.1 is not required for the original Figure 1-10 design. Ranking
# products remain available for appendix/ranking analyses.
# ==============================================================================

library(terra)
library(data.table)
library(ggplot2)

rm(list = ls())
gc()

# 0. Settings ===================================================================
base_dir <- "H:/Jing/ecoChina2"
vis_dir <- file.path(base_dir, "visualization var threshold0.4")
fig_dir <- file.path(vis_dir, "figures")
tab_dir <- file.path(vis_dir, "tables")
assessment_dir <- file.path(base_dir, "assessment_var")
old_assessment_dir <- file.path(base_dir, "assessment")
chord_dir <- file.path(assessment_dir, "chord diagrams")
result_map_root <- file.path(base_dir, "result maps")
dual_root <- file.path(base_dir, "dual suit")
reference_file <- file.path(base_dir, "raster", "ecosys_ori.tif")
palette_file <- file.path(base_dir, "color_palette_China.csv")
population_lookup_file <- file.path(
  base_dir, "future tree niche var", "tables", "population_projection_lookup_var.csv"
)
dual_population_lookup_file <- file.path(
  base_dir, "future tree niche dual suitability var", "tables",
  "dual_population_projection_lookup_var.csv"
)
climate_train_file <- file.path(
  base_dir, "accuracy_climate_var", "climate_rf_var_accuracy_summary.csv"
)
soil_train_file <- file.path(base_dir, "accuracy_soil", "soil_rf_accuracy_summary.csv")
climate_test_file <- file.path(assessment_dir, "climate_test_zone_metrics_var.csv")
old_test_file <- file.path(old_assessment_dir, "rf_test_zone_metrics.csv")
saved_map_confusion_file <- file.path(assessment_dir, "normal_map_confusion_long_var.csv")
multiclass_reference_map_file <- file.path(
  result_map_root, "multiclass_rf", "assigned_zone_normal_multiclass_rf.tif"
)

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(chord_dir, recursive = TRUE, showWarnings = FALSE)

model_zoneID <- c(1:7, 9:50, 52:55)
unmodeled_zoneID <- c(8L, 51L)
reference_zoneID <- sort(unique(c(model_zoneID, unmodeled_zoneID)))
novel_value <- 99L
novel_threshold <- 0.4
tie_tol <- 1e-4

method_order <- c("rf_var", "mf_var")
method_labels <- c(rf_var = "Plain RF", mf_var = "Plain MF RF")
preferred_method <- "rf_var"
tree_plot_method <- preferred_method

future_order <- c(
  "2011-2040SSP245", "2041-2070SSP245", "2071-2100SSP245",
  "2011-2040SSP585", "2041-2070SSP585", "2071-2100SSP585"
)
scenario_order <- c("normal", future_order)
tree_plot_scenarios <- scenario_order

species_page_ncol <- 5L
species_page_nrow <- 2L
species_per_page <- species_page_ncol * species_page_nrow
display_max_cells <- 300000L
background_max_cells <- 120000L
fig_dpi <- 320
tree_point_alpha <- 0.30
tree_point_size <- 0.075
strict_mode <- TRUE
terraOptions(memfrac = 0.10)

# 1. Helpers ====================================================================
cat0 <- function(...) cat(..., "\n", sep = "")

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nchar(x) == 0, "unnamed", x)
}

require_file <- function(file, label = basename(file)) {
  if (!file.exists(file)) stop("Missing ", label, ": ", file)
  file
}

require_columns <- function(dt, cols, label) {
  missing_cols <- setdiff(cols, names(dt))
  if (length(missing_cols) > 0) {
    stop(label, " is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  invisible(TRUE)
}

div <- function(a, b) ifelse(is.finite(b) & b > 0, a / b, NA_real_)

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
      strip.background = element_rect(fill = "grey95", colour = "grey78"),
      strip.text = element_text(face = "bold", colour = "grey15"),
      axis.text = element_text(colour = "grey20"),
      axis.title = element_text(colour = "grey10", face = "bold"),
      plot.title = element_text(face = "bold"),
      panel.spacing.x = grid::unit(1.05, "lines"),
      panel.spacing.y = grid::unit(0.42, "lines")
    )
}

save_gg <- function(p, file, width = 8, height = 6) {
  ggsave(file, p, width = width, height = height, dpi = fig_dpi, bg = "white")
  cat0("[SAVED] ", file)
}

step_log <- list()
run_step <- function(step_name, code, optional = FALSE) {
  cat0("\n============================================================")
  cat0(step_name)
  cat0("============================================================")
  t0 <- Sys.time()
  tryCatch({
    force(code)
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
    step_log[[length(step_log) + 1L]] <<- data.table(
      step = step_name, status = "done", message = NA_character_, elapsed_sec = elapsed
    )
    cat0("[DONE] ", step_name, " | ", elapsed, " sec")
  }, error = function(e) {
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
    status_value <- if (optional) "optional_skipped" else "error"
    step_log[[length(step_log) + 1L]] <<- data.table(
      step = step_name, status = status_value,
      message = conditionMessage(e), elapsed_sec = elapsed
    )
    cat0("[", toupper(status_value), "] ", step_name)
    cat0("  ", conditionMessage(e))
  })
  invisible(NULL)
}

# 2. Palette ====================================================================
require_file(palette_file, "vegetation-zone palette")
palette <- fread(palette_file)
require_columns(palette, c("zoneID", "COLOR"), "color_palette_China.csv")
palette[, zoneID := as.integer(zoneID)]
palette[!grepl("^#", COLOR) & !is.na(COLOR), COLOR := paste0("#", COLOR)]

zone_color_vector <- function(values) {
  values <- sort(unique(as.integer(values)))
  values <- values[!is.na(values)]
  out <- setNames(rep(NA_character_, length(values)), as.character(values))
  matched <- match(values, palette$zoneID)
  good <- !is.na(matched)
  out[as.character(values[good])] <- palette$COLOR[matched[good]]
  missing_values <- values[is.na(out[as.character(values)])]
  if (length(missing_values) > 0) {
    fallback <- grDevices::hcl.colors(length(missing_values), palette = "Dark 3")
    names(fallback) <- as.character(missing_values)
    out[names(fallback)] <- fallback
  }
  if (as.character(novel_value) %in% names(out)) {
    out[as.character(novel_value)] <- "#333333"
  }
  out
}

category_lookup <- NULL
if ("category2" %in% names(palette)) {
  category_lookup <- unique(
    palette[zoneID %in% model_zoneID, .(zoneID, category2 = as.character(category2))],
    by = "zoneID"
  )
}

# 3. Raster helpers ==============================================================
assigned_map_file <- function(method, scenario) {
  file.path(
    result_map_root, method,
    paste0(
      "assigned_zone_", scenario, "_threshold", novel_threshold,
      "_tol", tie_tol, "_novel", novel_value, "_maskNA8_noNovelNormal.tif"
    )
  )
}

dual_file <- function(method, scenario, zone) {
  file.path(dual_root, method, scenario, paste0("dual_suitability_zone", zone, ".tif"))
}

scenario_label <- function(scenario) {
  out <- scenario
  out[out == "normal"] <- "Predicted reference"
  out <- gsub("SSP245$", " SSP245", out)
  out <- gsub("SSP585$", " SSP585", out)
  out
}

thin_raster_for_plot <- function(x, categorical = TRUE, max_cells = display_max_cells) {
  x <- x[[1]]
  if (ncell(x) <= max_cells) return(x)
  fact <- ceiling(sqrt(ncell(x) / max_cells))
  aggregate(x, fact = fact, fun = if (categorical) "modal" else "mean", na.rm = TRUE)
}

raster_to_plot_dt <- function(file, categorical = TRUE, allowed_values = NULL) {
  require_file(file, "raster panel")
  x <- thin_raster_for_plot(rast(file)[[1]], categorical = categorical)
  dt <- as.data.table(as.data.frame(x, xy = TRUE, na.rm = TRUE))
  if (nrow(dt) == 0) return(data.table())
  value_col <- setdiff(names(dt), c("x", "y"))[1]
  setnames(dt, value_col, "value")
  if (categorical) {
    dt[, value := as.integer(round(value))]
    if (!is.null(allowed_values)) dt <- dt[value %in% as.integer(allowed_values)]
  }
  dt[!is.na(value) & is.finite(value)]
}

plot_zone_panel_gg <- function(files, titles, outfile, ncol, width, height,
                               plot_title, allowed_values) {
  if (length(files) != length(titles)) stop("files and titles must have equal length.")
  plot_list <- list()
  for (panel_index in seq_along(files)) {
    dt <- raster_to_plot_dt(files[[panel_index]], TRUE, allowed_values)
    if (nrow(dt) == 0) stop("No valid raster cells for panel: ", titles[panel_index])
    dt[, panel := titles[panel_index]]
    plot_list[[length(plot_list) + 1L]] <- dt
  }
  plot_dt <- rbindlist(plot_list, fill = TRUE)
  plot_dt[, panel := factor(panel, levels = titles)]
  values <- sort(unique(plot_dt$value))
  colors <- zone_color_vector(values)
  plot_dt[, value_chr := factor(as.character(value), levels = as.character(values))]
  p <- ggplot(plot_dt, aes(x = x, y = y, fill = value_chr)) +
    geom_raster() + facet_wrap(~panel, ncol = ncol, drop = FALSE) +
    scale_fill_manual(values = colors, drop = FALSE, guide = "none") +
    coord_equal(expand = FALSE) + labs(title = plot_title, x = NULL, y = NULL) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 13),
          strip.text = element_text(face = "bold", size = 9.5))
  save_gg(p, outfile, width, height)
  invisible(plot_dt)
}

# 4. Performance data ============================================================
build_performance_long <- function() {
  require_file(climate_test_file, "final climate independent-test metrics")
  require_file(old_test_file, "soil independent-test metrics")
  
  climate_test <- fread(climate_test_file)
  climate_test <- climate_test[climate_test[["method"]] %in% method_order]
  require_columns(
    climate_test,
    c("method", "zone", "accuracy", "balanced_accuracy", "f1", "precision",
      "recall", "specificity", "tss", "auc"),
    "climate_test_zone_metrics_var.csv"
  )
  climate_test[, niche_type := "climate"]
  
  soil_test <- fread(old_test_file)
  require_columns(
    soil_test,
    c("method", "niche", "zone", "accuracy", "balanced_accuracy", "f1",
      "precision", "recall", "specificity", "tss", "auc"),
    "assessment/rf_test_zone_metrics.csv"
  )
  soil_test <- soil_test[
    tolower(as.character(soil_test[["niche"]])) == "soil" &
      soil_test[["method"]] %in% c("plain_rf", "plain_mf")
  ]
  soil_test[, method := fifelse(method == "plain_rf", "rf_var", "mf_var")]
  soil_test[, niche_type := "soil"]
  
  test_metric_names <- c(
    accuracy = "Test accuracy", balanced_accuracy = "Test balanced accuracy",
    auc = "Test AUC", f1 = "Test F1", precision = "Test precision",
    recall = "Test recall", specificity = "Test specificity", tss = "Test TSS"
  )
  make_test_long <- function(dt) {
    rbindlist(lapply(names(test_metric_names), function(metric_col) {
      data.table(
        niche_type = dt$niche_type, method = dt$method,
        zoneID = as.integer(dt$zone), metric = unname(test_metric_names[metric_col]),
        value = as.numeric(dt[[metric_col]])
      )
    }), fill = TRUE)
  }
  test_long <- rbindlist(list(make_test_long(climate_test), make_test_long(soil_test)), fill = TRUE)
  train_list <- list()
  
  if (file.exists(climate_train_file)) {
    climate_train <- fread(climate_train_file)
    require_columns(climate_train, c("zone", "model", "train_accuracy", "oob_accuracy"),
                    "climate_rf_var_accuracy_summary.csv")
    climate_train <- climate_train[climate_train[["model"]] %in% method_order]
    train_list[[length(train_list) + 1L]] <- rbindlist(list(
      data.table(niche_type = "climate", method = climate_train$model,
                 zoneID = as.integer(climate_train$zone), metric = "Train accuracy",
                 value = as.numeric(climate_train$train_accuracy)),
      data.table(niche_type = "climate", method = climate_train$model,
                 zoneID = as.integer(climate_train$zone), metric = "OOB accuracy",
                 value = as.numeric(climate_train$oob_accuracy))
    ), fill = TRUE)
  }
  
  if (file.exists(soil_train_file)) {
    soil_train <- fread(soil_train_file)
    require_columns(soil_train, c("zone", "model", "train_accuracy", "oob_accuracy"),
                    "soil_rf_accuracy_summary.csv")
    soil_train <- soil_train[soil_train[["model"]] %in% c("plain_rf", "plain_mf_rf")]
    soil_train[, final_method := fifelse(model == "plain_rf", "rf_var", "mf_var")]
    train_list[[length(train_list) + 1L]] <- rbindlist(list(
      data.table(niche_type = "soil", method = soil_train$final_method,
                 zoneID = as.integer(soil_train$zone), metric = "Train accuracy",
                 value = as.numeric(soil_train$train_accuracy)),
      data.table(niche_type = "soil", method = soil_train$final_method,
                 zoneID = as.integer(soil_train$zone), metric = "OOB accuracy",
                 value = as.numeric(soil_train$oob_accuracy))
    ), fill = TRUE)
  }
  
  out <- if (length(train_list) > 0) {
    rbindlist(c(train_list, list(test_long)), fill = TRUE)
  } else test_long
  out <- out[method %in% method_order & zoneID %in% model_zoneID & is.finite(value)]
  out[, method_label := unname(method_labels[method])]
  out[, order_tmp := match(method, method_order)]
  setorder(out, niche_type, order_tmp, zoneID, metric)
  out[, order_tmp := NULL]
  out[]
}

# 5. Direct reference-map comparison ============================================
get_map_file <- function(method) {
  if (method == "multiclass_rf") return(multiclass_reference_map_file)
  assigned_map_file(method, "normal")
}

build_reference_comparison <- function(methods, label) {
  require_file(reference_file, "original vegetation map")
  map_files <- setNames(vapply(methods, get_map_file, character(1)), methods)
  missing_methods <- names(map_files)[!file.exists(map_files)]
  if (length(missing_methods) > 0) {
    stop("Missing reference map(s) for ", label, ": ", paste(missing_methods, collapse = ", "))
  }
  
  original_raw <- rast(reference_file)[[1]]
  original <- subst(original_raw, from = model_zoneID, to = model_zoneID, others = NA)
  names(original) <- "original_zone"
  
  predictions <- lapply(methods, function(method_key) {
    x <- rast(map_files[[method_key]])[[1]]
    if (!compareGeom(original, x, stopOnError = FALSE)) {
      stop("Geometry mismatch for ", method_key)
    }
    x <- subst(x, from = model_zoneID, to = model_zoneID, others = NA)
    names(x) <- method_key
    x
  })
  names(predictions) <- methods
  
  common_valid <- !is.na(original)
  for (method_key in methods) common_valid <- common_valid & !is.na(predictions[[method_key]])
  compared_pixels <- as.numeric(global(common_valid, "sum", na.rm = TRUE)[1, 1])
  valid_original_pixels <- as.numeric(global(!is.na(original), "sum", na.rm = TRUE)[1, 1])
  if (!is.finite(compared_pixels) || compared_pixels <= 0) stop("No common valid pixels for ", label)
  
  original_common <- ifel(common_valid, original, NA)
  names(original_common) <- "original_zone"
  confusion_list <- list(); zone_list <- list(); overall_list <- list()
  
  for (method_key in methods) {
    pred_common <- ifel(common_valid, predictions[[method_key]], NA)
    names(pred_common) <- "predicted_zone"
    ct <- as.data.table(crosstab(c(original_common, pred_common), long = TRUE, useNA = FALSE))
    if (ncol(ct) != 3) stop("Unexpected crosstab structure for ", method_key)
    setnames(ct, names(ct), c("original_zone", "predicted_zone", "n"))
    ct[, `:=`(
      method = method_key,
      original_zone = as.integer(original_zone),
      predicted_zone = as.integer(predicted_zone),
      n = as.numeric(n)
    )]
    ct <- ct[
      original_zone %in% model_zoneID & predicted_zone %in% model_zoneID & n > 0
    ]
    total <- sum(ct$n)
    if (total != compared_pixels) stop("Crosstab total differs from common mask for ", method_key)
    
    zone_dt <- rbindlist(lapply(model_zoneID, function(zone_value) {
      tp <- ct[original_zone == zone_value & predicted_zone == zone_value, sum(n)]
      fn <- ct[original_zone == zone_value & predicted_zone != zone_value, sum(n)]
      fp <- ct[original_zone != zone_value & predicted_zone == zone_value, sum(n)]
      tn <- total - tp - fn - fp
      recall_value <- div(tp, tp + fn)
      specificity_value <- div(tn, tn + fp)
      precision_value <- div(tp, tp + fp)
      data.table(
        method = method_key, zone = zone_value,
        original_pixels = tp + fn, predicted_pixels = tp + fp,
        TP = tp, TN = tn, FP = fp, FN = fn,
        accuracy = div(tp + tn, total),
        balanced_accuracy = div(recall_value + specificity_value, 2),
        recall = recall_value, specificity = specificity_value,
        precision = precision_value,
        f1 = div(2 * precision_value * recall_value, precision_value + recall_value),
        tss = recall_value + specificity_value - 1
      )
    }), fill = TRUE)
    
    exact_accuracy <- ct[original_zone == predicted_zone, sum(n)] / total
    broad_accuracy <- NA_real_
    if (!is.null(category_lookup)) {
      category_map <- setNames(category_lookup$category2, as.character(category_lookup$zoneID))
      ct[, original_category := unname(category_map[as.character(original_zone)])]
      ct[, predicted_category := unname(category_map[as.character(predicted_zone)])]
      if (!anyNA(ct$original_category) && !anyNA(ct$predicted_category)) {
        broad_accuracy <- ct[original_category == predicted_category, sum(n)] / total
      }
    }
    
    overall_dt <- data.table(
      method = method_key,
      valid_original_pixels = valid_original_pixels,
      compared_pixels = compared_pixels,
      missing_predictions = valid_original_pixels - compared_pixels,
      coverage = compared_pixels / valid_original_pixels,
      exact_zone_accuracy = exact_accuracy,
      broad_category_accuracy = broad_accuracy,
      macro_balanced_accuracy = mean(zone_dt$balanced_accuracy, na.rm = TRUE),
      macro_recall = mean(zone_dt$recall, na.rm = TRUE),
      macro_specificity = mean(zone_dt$specificity, na.rm = TRUE),
      macro_precision = mean(zone_dt$precision, na.rm = TRUE),
      macro_f1 = mean(zone_dt$f1, na.rm = TRUE),
      macro_tss = mean(zone_dt$tss, na.rm = TRUE),
      source_map = map_files[[method_key]], comparison_scope = label
    )
    confusion_list[[method_key]] <- ct
    zone_list[[method_key]] <- zone_dt
    overall_list[[method_key]] <- overall_dt
  }
  
  list(
    overall = rbindlist(overall_list, fill = TRUE),
    zone = rbindlist(zone_list, fill = TRUE),
    confusion = rbindlist(confusion_list, fill = TRUE),
    files = map_files
  )
}

# 6. Table 1 and Figure 1 ========================================================
run_step("Table 1 and Figure 1 | Final plain RF/MF performance", {
  performance_long <- build_performance_long()
  fwrite(performance_long, file.path(tab_dir, "Table1_binary_RF_zone_level_metrics_long.csv"))
  performance_summary <- performance_long[, {
    n_value <- .N
    sd_value <- if (n_value > 1) sd(value) else 0
    .(n_zones = uniqueN(zoneID), mean = mean(value), sd = sd_value, se = sd_value / sqrt(n_value))
  }, by = .(niche_type, method, method_label, metric)]
  performance_summary[, `:=`(ci95_low = mean - 1.96 * se, ci95_high = mean + 1.96 * se)]
  fwrite(performance_summary, file.path(tab_dir, "Table1_binary_RF_workflow_summary_long.csv"))
  
  metric_label_map <- c(
    "accuracy" = "Accuracy", "balanced accuracy" = "Balanced accuracy", "AUC" = "AUC",
    "F1" = "F1", "precision" = "Precision", "recall" = "Recall",
    "specificity" = "Specificity", "TSS" = "TSS"
  )
  metric_levels <- c("Accuracy", "Balanced accuracy", "AUC", "F1", "Precision", "Recall", "Specificity", "TSS")
  
  summary_plot <- copy(performance_summary)
  summary_plot[, eval_set := sub(" .*", "", metric)]
  summary_plot[, metric_name := sub("^(OOB|Train|Test) ", "", metric)]
  summary_plot[, metric_label := unname(metric_label_map[metric_name])]
  summary_plot <- summary_plot[!is.na(metric_label)]
  
  zone_plot <- copy(performance_long)
  zone_plot[, eval_set := sub(" .*", "", metric)]
  zone_plot[, metric_name := sub("^(OOB|Train|Test) ", "", metric)]
  zone_plot[, metric_label := unname(metric_label_map[metric_name])]
  zone_plot <- zone_plot[!is.na(metric_label)]
  
  summary_plot[, `:=`(
    metric_label = factor(metric_label, levels = rev(metric_levels)),
    eval_set = factor(eval_set, levels = c("OOB", "Train", "Test")),
    method_label = factor(method_label, levels = method_labels[method_order]),
    niche_type = factor(niche_type, levels = c("climate", "soil"))
  )]
  zone_plot[, `:=`(
    metric_label = factor(metric_label, levels = rev(metric_levels)),
    eval_set = factor(eval_set, levels = c("OOB", "Train", "Test")),
    method_label = factor(method_label, levels = method_labels[method_order]),
    niche_type = factor(niche_type, levels = c("climate", "soil")),
    zoneID_chr = factor(as.character(zoneID), levels = as.character(model_zoneID))
  )]
  
  finite_values <- c(zone_plot$value, summary_plot$ci95_low, summary_plot$ci95_high)
  finite_values <- finite_values[is.finite(finite_values)]
  x_lower <- max(-1, min(0.60, floor(min(finite_values) / 0.05) * 0.05))
  
  p1a <- ggplot(summary_plot, aes(x = mean, y = metric_label)) +
    geom_errorbarh(aes(xmin = ci95_low, xmax = ci95_high), height = 0.22,
                   linewidth = 0.75, colour = "grey35") +
    geom_point(aes(fill = niche_type), shape = 21, size = 2.8, colour = "grey12") +
    facet_grid(niche_type + eval_set ~ method_label, scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = c(climate = "#2E5E8C", soil = "#9A7B32"), guide = "none") +
    scale_x_continuous(limits = c(x_lower, 1), breaks = seq(x_lower, 1, by = 0.1)) +
    labs(
      title = "Figure 1a. One-hot binary RF performance across vegetation zones",
      subtitle = "Points are across-zone means; horizontal lines show mean +/- 1.96 SE.",
      x = "Metric value", y = NULL
    ) + theme_fig1(base_size = 10.2)
  save_gg(p1a, file.path(fig_dir, "Figure1a_onehot_binary_RF_summary.png"), 10.5, 9.6)
  save_gg(p1a, file.path(fig_dir, "Figure1_onehot_binary_RF_OOB_train_test_performance_dotrange.png"), 10.5, 9.6)
  
  zone_colors <- zone_color_vector(model_zoneID)
  p1b <- ggplot(zone_plot, aes(x = value, y = metric_label)) +
    geom_point(aes(colour = zoneID_chr), position = position_jitter(width = 0, height = 0.14),
               size = 1.75, alpha = 0.75, show.legend = FALSE) +
    geom_errorbarh(data = summary_plot, aes(xmin = ci95_low, xmax = ci95_high, y = metric_label),
                   inherit.aes = FALSE, height = 0.22, linewidth = 0.72, colour = "grey15") +
    geom_point(data = summary_plot, aes(x = mean, y = metric_label, fill = niche_type),
               inherit.aes = FALSE, shape = 23, size = 2.8, colour = "grey10", show.legend = FALSE) +
    facet_grid(niche_type + eval_set ~ method_label, scales = "free_y", space = "free_y", switch = "y") +
    scale_colour_manual(values = zone_colors, guide = "none") +
    scale_fill_manual(values = c(climate = "#2E5E8C", soil = "#9A7B32"), guide = "none") +
    scale_x_continuous(limits = c(x_lower, 1), breaks = seq(x_lower, 1, by = 0.1)) +
    labs(
      title = "Figure 1b. Zone-level variation in one-hot binary RF performance",
      subtitle = "Small points are vegetation zones; diamonds and lines show the across-zone mean and mean +/- 1.96 SE.",
      x = "Metric value", y = NULL
    ) + theme_fig1(base_size = 10)
  save_gg(p1b, file.path(fig_dir, "Figure1b_onehot_binary_RF_zone_level_metrics.png"), 10.5, 10.4)
})

# 7. Recalculate binary reference comparison ====================================
binary_comparison <- NULL
run_step("Reference-map recalculation | RF and MF common valid mask", {
  binary_comparison <<- build_reference_comparison(
    method_order, "Common valid mask of final RF and MF reference maps"
  )
  fwrite(binary_comparison$overall, file.path(tab_dir, "Binary_reference_map_recalculated_overall_metrics.csv"))
  fwrite(binary_comparison$zone, file.path(tab_dir, "Binary_reference_map_recalculated_zone_metrics.csv"))
  fwrite(binary_comparison$confusion, file.path(tab_dir, "Binary_reference_map_recalculated_confusion_long.csv"))
})

# 8. Figure 2 ===================================================================
run_step("Figure 2 | Predicted reference maps", {
  letters_final <- c("a", "b")
  for (method_index in seq_along(method_order)) {
    method_key <- method_order[method_index]
    map_file <- assigned_map_file(method_key, "normal")
    require_file(map_file, paste0(method_key, " predicted reference map"))
    plot_zone_panel_gg(
      files = list(map_file), titles = method_labels[method_key],
      outfile = file.path(fig_dir, paste0("Figure2", letters_final[method_index],
                                          "_predicted_reference_map_", method_key, ".png")),
      ncol = 1, width = 6.2, height = 4.8,
      plot_title = paste0("Predicted reference map: ", method_labels[method_key]),
      allowed_values = model_zoneID
    )
  }
  plot_zone_panel_gg(
    files = list(reference_file), titles = "Original vegetation map",
    outfile = file.path(fig_dir, "Figure2_reference_vegetation_map.png"),
    ncol = 1, width = 6.2, height = 4.8, plot_title = "Original vegetation map",
    allowed_values = reference_zoneID
  )
})

# 9. Table 2 ====================================================================
run_step("Table 2 | Reference-map reconstruction performance", {
  if (is.null(binary_comparison)) stop("Binary reference comparison was not generated.")
  table2 <- copy(binary_comparison$overall)
  table2[, method_label := unname(method_labels[method])]
  table2[, `:=`(
    coverage_percent = 100 * coverage,
    exact_zone_agreement_percent = 100 * exact_zone_accuracy,
    broad_category_agreement_percent = 100 * broad_category_accuracy
  )]
  fwrite(table2, file.path(tab_dir, "Table2_reference_map_accuracy_binary_workflows.csv"))
})

# 10. Figure 3 ==================================================================
run_step("Figure 3 | Climate and soil zone-level metrics", {
  performance_long <- build_performance_long()
  metrics_keep <- c("Test F1", "Test precision", "Test recall")
  plot_dt <- performance_long[metric %in% metrics_keep]
  if (nrow(plot_dt) == 0) stop("No Test F1/precision/recall records are available.")
  plot_dt[, `:=`(
    metric = factor(metric, levels = metrics_keep),
    method_label = factor(method_label, levels = method_labels[method_order]),
    zoneID_fac = factor(as.character(zoneID), levels = as.character(rev(model_zoneID))),
    zoneID_chr = as.character(zoneID)
  )]
  zone_colors <- zone_color_vector(model_zoneID)
  for (niche_value in c("climate", "soil")) {
    dt_sub <- plot_dt[niche_type == niche_value]
    p <- ggplot(dt_sub, aes(x = value, y = zoneID_fac, colour = zoneID_chr)) +
      geom_point(size = 1.5, alpha = 0.75) +
      facet_grid(metric ~ method_label, scales = "free_x") +
      scale_colour_manual(values = zone_colors, guide = "none") +
      labs(
        title = paste0(tools::toTitleCase(niche_value), " binary RF zone-level performance"),
        subtitle = "Each point is one vegetation zone; colors follow the vegetation-zone palette.",
        x = "Metric value", y = "Vegetation zone"
      ) + theme_ms(base_size = 9.5) +
      theme(axis.text.y = element_text(size = 6.2), axis.text.x = element_text(size = 7),
            panel.grid.major.y = element_blank())
    save_gg(
      p,
      file.path(fig_dir, paste0("Figure3", ifelse(niche_value == "climate", "a", "b"),
                                "_", niche_value, "_binary_RF_zone_metrics_all_workflows.png")),
      8.5, 8.5
    )
  }
})

# 11. Figure 4 ==================================================================
run_step("Figure 4 | Major ecotype confusion flows", {
  if (is.null(binary_comparison)) stop("Binary reference comparison was not generated.")
  confusion_dt <- binary_comparison$confusion[
    binary_comparison$confusion[["method"]] == preferred_method &
      binary_comparison$confusion[["original_zone"]] != binary_comparison$confusion[["predicted_zone"]]
  ]
  flow <- confusion_dt[, .(count = sum(n)), by = .(from = original_zone, to = predicted_zone)]
  setorder(flow, -count, from, to)
  if (nrow(flow) == 0) stop("No off-diagonal confusion flows found.")
  flow <- flow[seq_len(min(20L, nrow(flow)))]
  flow[, `:=`(
    from_chr = as.character(from), to_chr = as.character(to),
    flow_label = paste0("Zone ", from, " -> Zone ", to)
  )]
  fwrite(flow, file.path(tab_dir, "Figure4_major_ecotype_confusion_flows.csv"))
  colors <- zone_color_vector(flow$from)
  
  if (requireNamespace("ggalluvial", quietly = TRUE)) {
    p <- ggplot(flow, aes(y = count, axis1 = from_chr, axis2 = to_chr)) +
      ggalluvial::geom_alluvium(aes(fill = from_chr), width = 1/12, alpha = 0.78, show.legend = FALSE) +
      ggalluvial::geom_stratum(width = 1/8, fill = "grey94", colour = "grey40") +
      ggalluvial::stat_stratum(
        geom = "text",
        aes(label = after_stat(stratum)),
        size = 2.6
      ) +
      scale_fill_manual(values = colors) +
      scale_x_discrete(limits = c("Original", "Predicted"), expand = c(0.08, 0.08)) +
      labs(title = paste0("Major reference-map confusion flows: ", method_labels[preferred_method]),
           x = NULL, y = "Pixel count") + theme_ms()
  } else {
    flow[, flow_label := factor(flow_label, levels = rev(flow_label))]
    p <- ggplot(flow, aes(x = count, y = flow_label, fill = from_chr)) +
      geom_col(width = 0.72, alpha = 0.80) +
      scale_fill_manual(values = colors, guide = "none") +
      labs(title = paste0("Major reference-map confusion flows: ", method_labels[preferred_method]),
           x = "Pixel count", y = NULL) + theme_ms(base_size = 10)
  }
  save_gg(p, file.path(fig_dir, "Figure4_major_ecotype_confusion_flows.png"), 8.8, 6.6)
})

# 12. Figure 5 ==================================================================
run_step("Figure 5 | Ecosystem niche maps, 2 x 4", {
  panel_titles <- c(
    "Original", "SSP245\n2011-2040", "SSP245\n2041-2070", "SSP245\n2071-2100",
    "Predicted\nreference", "SSP585\n2011-2040", "SSP585\n2041-2070", "SSP585\n2071-2100"
  )
  figure_letters <- c("a", "b")
  for (method_index in seq_along(method_order)) {
    method_key <- method_order[method_index]
    panel_files <- c(
      reference_file,
      assigned_map_file(method_key, "2011-2040SSP245"),
      assigned_map_file(method_key, "2041-2070SSP245"),
      assigned_map_file(method_key, "2071-2100SSP245"),
      assigned_map_file(method_key, "normal"),
      assigned_map_file(method_key, "2011-2040SSP585"),
      assigned_map_file(method_key, "2041-2070SSP585"),
      assigned_map_file(method_key, "2071-2100SSP585")
    )
    missing_files <- panel_files[!file.exists(panel_files)]
    if (length(missing_files) > 0) {
      stop("Figure 5 missing maps for ", method_key, ":\n", paste(missing_files, collapse = "\n"))
    }
    plot_zone_panel_gg(
      files = as.list(panel_files), titles = panel_titles,
      outfile = file.path(fig_dir, paste0("Figure5", figure_letters[method_index],
                                          "_ecosystem_niche_maps_", method_key, "_2x4.png")),
      ncol = 4, width = 13.2, height = 6.8,
      plot_title = paste0("Ecosystem niche projection: ", method_labels[method_key]),
      allowed_values = c(reference_zoneID, novel_value)
    )
  }
})

# 13. Figure 6 ==================================================================
area_by_zone <- function(map_file, cell_area) {
  require_file(map_file, "assigned map")
  z <- rast(map_file)[[1]]
  area_table <- as.data.table(zonal(cell_area, z, fun = "sum", na.rm = TRUE))
  if (ncol(area_table) < 2) stop("Unexpected zonal-area output for ", map_file)
  setnames(area_table, names(area_table)[1:2], c("zoneID", "area_km2"))
  area_table[, `:=`(zoneID = as.integer(zoneID), area_km2 = as.numeric(area_km2))]
  area_table[zoneID %in% c(model_zoneID, novel_value)]
}

encode_transition_stack <- function(stack) {
  n_layer <- nlyr(stack)
  weights <- 100^(0:(n_layer - 1))
  app(stack, fun = function(x) {
    if (is.null(dim(x))) {
      if (anyNA(x)) return(NA_real_)
      return(sum(x * weights))
    }
    out <- rep(NA_real_, nrow(x))
    good <- complete.cases(x)
    if (any(good)) out[good] <- as.numeric(x[good, , drop = FALSE] %*% weights)
    out
  })
}

decode_transition_code <- function(code, n_stage) {
  out <- matrix(NA_integer_, nrow = length(code), ncol = n_stage)
  for (stage_index in seq_len(n_stage)) {
    out[, stage_index] <- as.integer(floor(code / (100^(stage_index - 1))) %% 100)
  }
  out
}

build_transition_paths <- function(files, stage_names, family, cell_area) {
  stack <- rast(files)
  encoded <- encode_transition_stack(stack)
  area_table <- as.data.table(zonal(cell_area, encoded, fun = "sum", na.rm = TRUE))
  if (ncol(area_table) < 2) stop("Unexpected transition-area table.")
  setnames(area_table, names(area_table)[1:2], c("transition_code", "area_km2"))
  decoded <- decode_transition_code(as.numeric(area_table$transition_code), length(stage_names))
  decoded_dt <- as.data.table(decoded)
  setnames(decoded_dt, paste0("stage", seq_along(stage_names)))
  out <- cbind(decoded_dt, data.table(area_km2 = as.numeric(area_table$area_km2), family = family))
  valid_values <- c(model_zoneID, novel_value)
  keep <- rep(TRUE, nrow(out))
  for (stage_col in paste0("stage", seq_along(stage_names))) {
    keep <- keep & out[[stage_col]] %in% valid_values
  }
  out[keep]
}

run_step("Figure 6 | Area change and zone-transition trajectories", {
  reference_map <- rast(reference_file)[[1]]
  cell_area <- cellSize(reference_map, unit = "km")
  map_files <- setNames(vapply(scenario_order, function(scenario_value) {
    assigned_map_file(preferred_method, scenario_value)
  }, character(1)), scenario_order)
  missing_files <- map_files[!file.exists(map_files)]
  if (length(missing_files) > 0) {
    stop("Figure 6 missing assigned maps:\n", paste(missing_files, collapse = "\n"))
  }
  
  area_list <- lapply(names(map_files), function(scenario_value) {
    out <- area_by_zone(map_files[[scenario_value]], cell_area)
    out[, scenario := scenario_value]
    out
  })
  area_dt <- rbindlist(area_list, fill = TRUE)
  area_grid <- CJ(scenario = scenario_order, zoneID = c(model_zoneID, novel_value), unique = TRUE)
  area_dt <- merge(area_grid, area_dt, by = c("scenario", "zoneID"), all.x = TRUE, sort = FALSE)
  area_dt[is.na(area_km2), area_km2 := 0]
  area_dt[, scenario_label := scenario_label(scenario)]
  fwrite(area_dt, file.path(tab_dir, "Figure6_all_zone_area_by_scenario.csv"))
  
  area_plot <- copy(area_dt)
  area_plot[, `:=`(
    area_million_km2 = area_km2 / 1e6,
    scenario_label = factor(scenario_label, levels = scenario_label(scenario_order)),
    zone_label = fifelse(zoneID == novel_value, "Novel", as.character(zoneID))
  )]
  area_plot[, zone_label := factor(zone_label, levels = c(as.character(model_zoneID), "Novel"))]
  fill_colors <- zone_color_vector(c(model_zoneID, novel_value))
  names(fill_colors)[names(fill_colors) == as.character(novel_value)] <- "Novel"
  
  p6a <- ggplot(area_plot, aes(x = scenario_label, y = area_million_km2, fill = zone_label)) +
    geom_col(width = 0.58) + scale_fill_manual(values = fill_colors, drop = FALSE) +
    labs(
      title = paste0("Figure 6a. Projected ecotype area change: ", method_labels[preferred_method]),
      x = NULL, y = expression("Area (" * 10^6 * " km"^2 * ")"), fill = "Zone"
    ) + theme_ms() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.text = element_text(size = 6.2),
          legend.key.size = grid::unit(0.30, "cm")) +
    guides(fill = guide_legend(ncol = 3, byrow = TRUE))
  save_gg(p6a, file.path(fig_dir, "Figure6a_all_zone_area_change.png"), 10.2, 6.4)
  save_gg(p6a, file.path(fig_dir, "Figure6_all_zone_area_change.png"), 10.2, 6.4)
  
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    stop("Package 'ggalluvial' is required for Figure 6b.")
  }
  stage_names <- c("Predicted reference", "2011-2040", "2041-2070", "2071-2100")
  files_245 <- c(
    assigned_map_file(preferred_method, "normal"),
    assigned_map_file(preferred_method, "2011-2040SSP245"),
    assigned_map_file(preferred_method, "2041-2070SSP245"),
    assigned_map_file(preferred_method, "2071-2100SSP245")
  )
  files_585 <- c(
    assigned_map_file(preferred_method, "normal"),
    assigned_map_file(preferred_method, "2011-2040SSP585"),
    assigned_map_file(preferred_method, "2041-2070SSP585"),
    assigned_map_file(preferred_method, "2071-2100SSP585")
  )
  path_dt <- rbindlist(list(
    build_transition_paths(files_245, stage_names, "SSP245", cell_area),
    build_transition_paths(files_585, stage_names, "SSP585", cell_area)
  ), fill = TRUE)
  fwrite(path_dt, file.path(tab_dir, "Figure6b_zone_transition_paths.csv"))
  path_dt[, path_id := .I]
  long_path <- melt(
    path_dt,
    id.vars = c("path_id", "area_km2", "family"),
    measure.vars = paste0("stage", seq_along(stage_names)),
    variable.name = "stage_var", value.name = "zoneID"
  )
  long_path[, stage_index := as.integer(sub("stage", "", stage_var))]
  long_path[, stage := factor(stage_names[stage_index], levels = stage_names)]
  long_path[, zone_label := fifelse(zoneID == novel_value, "Novel", as.character(zoneID))]
  long_path[, area_million_km2 := area_km2 / 1e6]
  trajectory_colors <- zone_color_vector(c(model_zoneID, novel_value))
  names(trajectory_colors)[names(trajectory_colors) == as.character(novel_value)] <- "Novel"
  
  p6b <- ggplot(long_path, aes(
    x = stage, stratum = zone_label, alluvium = path_id,
    y = area_million_km2, fill = zone_label
  )) +
    ggalluvial::geom_flow(alpha = 0.55, colour = NA) +
    ggalluvial::geom_stratum(width = 0.22, colour = "grey35", linewidth = 0.15) +
    facet_grid(family ~ .) + scale_fill_manual(values = trajectory_colors, guide = "none") +
    labs(
      title = paste0("Figure 6b. Ecotype transition trajectories: ", method_labels[preferred_method]),
      subtitle = "Rows show SSP245 and SSP585. Flows track mapped area from the predicted reference through three future periods.",
      x = NULL, y = expression("Transition area (" * 10^6 * " km"^2 * ")")
    ) + theme_ms(base_size = 10) +
    theme(strip.text.y = element_text(angle = 0, face = "bold"),
          axis.text.x = element_text(face = "bold"))
  save_gg(p6b, file.path(fig_dir, "Figure6b_zone_transition_sankey.png"), 12.8, 8.8)
})

# 14. Table 3 and Figure 7 =======================================================
shannon_H <- function(abundance) {
  abundance <- as.numeric(abundance)
  abundance <- abundance[is.finite(abundance) & abundance > 0]
  if (length(abundance) == 0) return(NA_real_)
  p <- abundance / sum(abundance)
  -sum(p * log(p))
}

run_step("Table 3 and Figure 7 | Reference population abundance", {
  require_file(population_lookup_file, "population projection lookup from script 6.1")
  population <- fread(population_lookup_file)
  require_columns(
    population,
    c("Species", "source_zone", "reference_abundance", "projected"),
    "population_projection_lookup_var.csv"
  )
  population[, `:=`(
    species = as.character(Species), zoneID = as.integer(source_zone),
    abundance = as.numeric(reference_abundance)
  )]
  table3 <- population[is.finite(abundance) & abundance > 0, .(
    n_reference_populations = .N,
    reference_total_abundance = sum(abundance),
    reference_shannon_H = shannon_H(abundance)
  ), by = species]
  eligible <- population[
    projected == TRUE & zoneID %in% model_zoneID & is.finite(abundance) & abundance > 0,
    .(
      n_projection_eligible_populations = .N,
      projection_eligible_total_abundance = sum(abundance),
      projection_eligible_shannon_H = shannon_H(abundance)
    ), by = species
  ]
  table3 <- merge(table3, eligible, by = "species", all = TRUE, sort = FALSE)
  setorder(table3, -n_reference_populations, -reference_shannon_H, species)
  table3_top <- table3[seq_len(min(10L, nrow(table3)))]
  fwrite(table3_top, file.path(tab_dir, "Table3_ten_species_populations_ShannonH.csv"))
  fwrite(table3, file.path(tab_dir, "Table3_all_species_population_diversity.csv"))
  pop_long <- population[, .(species, zoneID, abundance, projection_eligible = projected)]
  fwrite(pop_long, file.path(tab_dir, "Table3_species_population_abundance_long.csv"))
  
  species_keep <- table3_top$species
  observed <- pop_long[species %in% species_keep]
  zones_keep <- sort(unique(observed[abundance > 0, zoneID]))
  heat <- CJ(species = species_keep, zoneID = zones_keep, unique = TRUE)
  heat <- merge(heat, observed[, .(species, zoneID, abundance)],
                by = c("species", "zoneID"), all.x = TRUE, sort = FALSE)
  heat[, `:=`(
    log10_abundance = fifelse(is.finite(abundance), log10(abundance + 1), NA_real_),
    species = factor(species, levels = rev(species_keep)),
    zoneID_fac = factor(as.character(zoneID), levels = as.character(zones_keep))
  )]
  axis_labels <- setNames(
    ifelse(zones_keep %in% unmodeled_zoneID, paste0(zones_keep, "*"), as.character(zones_keep)),
    as.character(zones_keep)
  )
  p7 <- ggplot(heat, aes(x = zoneID_fac, y = species, fill = log10_abundance)) +
    geom_tile(colour = "white", linewidth = 0.18) +
    scale_x_discrete(labels = axis_labels) +
    scale_fill_gradient(low = "#F5F5F2", high = "#3F4A4A",
                        name = "log10(reference\nabundance + 1)", na.value = "white") +
    labs(
      title = "Species-by-ecotype reference population abundance",
      x = "Ecotype / vegetation zone", y = NULL,
      caption = "* Zones 8 and 51 occur in reference data but are excluded from future ecosystem projections."
    ) + theme_ms() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.5),
          axis.text.y = element_text(size = 8), panel.grid = element_blank())
  save_gg(p7, file.path(fig_dir, "Figure7_species_ecotype_population_abundance_heatmap.png"), 10.5, 5.8)
})

# 15. Figure 8 and Figure 9 helpers =============================================
reference_plot_background <- function(max_cells = background_max_cells) {
  ref <- rast(reference_file)[[1]]
  ref_plot <- thin_raster_for_plot(ref, categorical = TRUE, max_cells = max_cells)
  bg <- as.data.table(as.data.frame(ref_plot, xy = TRUE, na.rm = TRUE))
  if (nrow(bg) == 0) stop("Reference raster has no cells for map background.")
  list(
    data = bg[, .(x, y)],
    xlim = c(xmin(ref), xmax(ref)),
    ylim = c(ymin(ref), ymax(ref))
  )
}

dual_to_points <- function(file, source_zone, threshold = novel_threshold) {
  require_file(file, "dual-suitability raster")
  x <- rast(file)[[1]]
  ref <- rast(reference_file)[[1]]
  if (!compareGeom(x, ref, stopOnError = FALSE)) x <- resample(x, ref, method = "bilinear")
  x <- mask(x, ref)
  x <- ifel(!is.na(x) & x >= threshold, x, NA)
  if (ncell(x) > display_max_cells) {
    fact <- ceiling(sqrt(ncell(x) / display_max_cells))
    x <- aggregate(x, fact = fact, fun = "max", na.rm = TRUE)
    x <- ifel(!is.na(x) & x >= threshold, x, NA)
  }
  dt <- as.data.table(as.data.frame(x, xy = TRUE, na.rm = TRUE))
  if (nrow(dt) == 0) return(data.table())
  value_col <- setdiff(names(dt), c("x", "y"))[1]
  setnames(dt, value_col, "dual_suitability")
  dt <- dt[is.finite(dual_suitability) & dual_suitability >= threshold]
  dt[, source_zone := as.integer(source_zone)]
  dt[]
}

build_dual_cache <- function(method, scenario, zones) {
  zones <- sort(unique(as.integer(zones)))
  cache <- lapply(zones, function(zone_value) {
    dual_to_points(dual_file(method, scenario, zone_value), zone_value, novel_threshold)
  })
  names(cache) <- as.character(zones)
  cache
}

read_projection_population <- function() {
  require_file(dual_population_lookup_file, "6.21 dual population lookup")
  population <- fread(dual_population_lookup_file)
  require_columns(population, c("Species", "source_zone"),
                  "dual_population_projection_lookup_var.csv")
  population[, `:=`(species = as.character(Species), source_zone = as.integer(source_zone))]
  population <- population[source_zone %in% model_zoneID]
  population[]
}

plot_species_population_dual_pages <- function(population, method, scenario, outfile_prefix) {
  species_order <- sort(unique(population$species))
  zones_needed <- sort(unique(population$source_zone))
  cache <- build_dual_cache(method, scenario, zones_needed)
  map_bg <- reference_plot_background()
  page_id <- ceiling(seq_along(species_order) / species_per_page)
  zone_colors <- zone_color_vector(zones_needed)
  
  for (page_value in sort(unique(page_id))) {
    species_selected <- species_order[page_id == page_value]
    plot_list <- list()
    for (species_value in species_selected) {
      species_zones <- population[species == species_value, sort(unique(source_zone))]
      dts <- cache[as.character(species_zones)]
      dts <- dts[lengths(dts) > 0]
      if (length(dts) > 0) {
        dt <- rbindlist(dts, fill = TRUE)
        if (nrow(dt) > 0) {
          dt[, species := species_value]
          plot_list[[length(plot_list) + 1L]] <- dt
        }
      }
    }
    plot_dt <- if (length(plot_list) > 0) {
      rbindlist(plot_list, fill = TRUE)
    } else {
      data.table(x = numeric(), y = numeric(), dual_suitability = numeric(),
                 source_zone = integer(), species = character())
    }
    plot_dt[, `:=`(
      species = factor(species, levels = species_selected),
      source_zone_chr = factor(as.character(source_zone), levels = as.character(zones_needed))
    )]
    panel_frame <- data.table(
      species = factor(species_selected, levels = species_selected),
      x = mean(map_bg$xlim), y = mean(map_bg$ylim)
    )
    p <- ggplot() +
      geom_raster(data = map_bg$data, aes(x = x, y = y), fill = "grey95") +
      geom_blank(data = panel_frame, aes(x = x, y = y)) +
      geom_point(data = plot_dt, aes(x = x, y = y, colour = source_zone_chr),
                 size = tree_point_size, alpha = tree_point_alpha) +
      facet_wrap(~species, ncol = species_page_ncol, drop = FALSE) +
      coord_equal(xlim = map_bg$xlim, ylim = map_bg$ylim, expand = FALSE) +
      scale_colour_manual(values = zone_colors, drop = FALSE) +
      labs(
        title = paste0("Tree population dual-niche projections: ", method_labels[method], " | ", scenario),
        subtitle = paste0(
          "Full China extent. Suitable cells have dual suitability >= ", novel_threshold,
          "; color indicates source zone; point opacity = 30%."
        ),
        colour = "Source zone"
      ) + theme_void(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9), strip.text = element_text(face = "bold", size = 9.2),
        legend.position = "right", legend.text = element_text(size = 6.6),
        legend.key.size = grid::unit(0.30, "cm")
      ) +
      guides(colour = guide_legend(ncol = 2, byrow = TRUE,
                                   override.aes = list(size = 2, alpha = 1)))
    outfile <- paste0(
      outfile_prefix, "_", safe_name(method), "_", safe_name(scenario),
      "_page", sprintf("%02d", page_value), ".png"
    )
    save_gg(p, outfile, 15.5, 7.0)
  }
}

plot_species_level_dual_pages <- function(population, method, scenario, outfile_prefix) {
  species_order <- sort(unique(population$species))
  zones_needed <- sort(unique(population$source_zone))
  cache <- build_dual_cache(method, scenario, zones_needed)
  map_bg <- reference_plot_background()
  page_id <- ceiling(seq_along(species_order) / species_per_page)
  
  for (page_value in sort(unique(page_id))) {
    species_selected <- species_order[page_id == page_value]
    plot_list <- list()
    for (species_value in species_selected) {
      species_zones <- population[species == species_value, sort(unique(source_zone))]
      dts <- cache[as.character(species_zones)]
      dts <- dts[lengths(dts) > 0]
      if (length(dts) > 0) {
        dt <- rbindlist(dts, fill = TRUE)
        if (nrow(dt) > 0) {
          dt <- dt[, .(dual_suitability = max(dual_suitability, na.rm = TRUE)), by = .(x, y)]
          dt[, species := species_value]
          plot_list[[length(plot_list) + 1L]] <- dt
        }
      }
    }
    plot_dt <- if (length(plot_list) > 0) {
      rbindlist(plot_list, fill = TRUE)
    } else {
      data.table(x = numeric(), y = numeric(), dual_suitability = numeric(), species = character())
    }
    plot_dt[, species := factor(species, levels = species_selected)]
    panel_frame <- data.table(
      species = factor(species_selected, levels = species_selected),
      x = mean(map_bg$xlim), y = mean(map_bg$ylim)
    )
    p <- ggplot() +
      geom_raster(data = map_bg$data, aes(x = x, y = y), fill = "grey95") +
      geom_blank(data = panel_frame, aes(x = x, y = y)) +
      geom_point(data = plot_dt, aes(x = x, y = y, colour = dual_suitability),
                 size = tree_point_size, alpha = tree_point_alpha) +
      facet_wrap(~species, ncol = species_page_ncol, drop = FALSE) +
      coord_equal(xlim = map_bg$xlim, ylim = map_bg$ylim, expand = FALSE) +
      scale_colour_gradientn(
        colours = grDevices::hcl.colors(60, palette = "YlGnBu"),
        limits = c(novel_threshold, 1), oob = scales::squish,
        name = "Max dual\nsuitability"
      ) +
      labs(
        title = paste0("Species-level dual-niche projections: ", method_labels[method], " | ", scenario),
        subtitle = paste0(
          "Full China extent. Species niche is the pixel-wise maximum across population dual niches; suitable cells >= ",
          novel_threshold, " are shown with 30% opacity."
        )
      ) + theme_void(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9), strip.text = element_text(face = "bold", size = 9.2),
        legend.position = "right", legend.text = element_text(size = 7),
        legend.key.size = grid::unit(0.35, "cm")
      )
    outfile <- paste0(
      outfile_prefix, "_", safe_name(method), "_", safe_name(scenario),
      "_page", sprintf("%02d", page_value), ".png"
    )
    save_gg(p, outfile, 15.5, 7.0)
  }
}

# 16. Figure 8 ==================================================================
run_step("Figure 8 | Tree population dual-niche projection maps", {
  population <- read_projection_population()
  for (scenario_value in tree_plot_scenarios) {
    plot_species_population_dual_pages(
      population, tree_plot_method, scenario_value,
      file.path(fig_dir, "Figure8_tree_population_dual_niche_projection")
    )
  }
  fwrite(population, file.path(tab_dir, "Figure8_species_population_source_zones.csv"))
})

# 17. Figure 9 ==================================================================
run_step("Figure 9 | Species-level dual-niche projection maps", {
  population <- read_projection_population()
  for (scenario_value in tree_plot_scenarios) {
    plot_species_level_dual_pages(
      population, tree_plot_method, scenario_value,
      file.path(fig_dir, "Figure9_species_level_dual_niche_projection")
    )
  }
  fwrite(population, file.path(tab_dir, "Figure9_species_population_source_zones.csv"))
})

# 18. Chord diagrams =============================================================
plot_chord <- function(flow, item_order, item_color, item_label, out_file, main,
                       label_cex = 0.65, label_track_height = 0.13) {
  flow <- as.data.table(copy(flow))
  flow[, `:=`(from = as.character(from), to = as.character(to), n = as.numeric(n))]
  flow <- flow[n > 0]
  from_id <- item_order[item_order %in% unique(flow$from)]
  to_id <- rev(item_order[item_order %in% unique(flow$to)])
  from_sector <- paste0("O_", from_id)
  to_sector <- paste0("P_", to_id)
  sector_order <- c(from_sector, to_sector)
  flow[, `:=`(from_sector = paste0("O_", from), to_sector = paste0("P_", to))]
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
  on.exit({circlize::circos.clear(); dev.off()}, add = TRUE)
  par(mar = c(0.5, 0.5, 2.8, 0.5), xpd = NA)
  circlize::circos.par(
    start.degree = 90, clock.wise = FALSE, cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.002, 0.002), canvas.xlim = c(-1.38, 1.38),
    canvas.ylim = c(-1.28, 1.28), points.overflow.warning = FALSE
  )
  circlize::chordDiagram(
    x = as.data.frame(flow[, .(from_sector, to_sector, n)]),
    order = sector_order, grid.col = sector_color, grid.border = NA,
    col = link_color, transparency = 0.72, directional = 1,
    direction.type = "diffHeight", diffHeight = circlize::mm_h(1),
    link.target.prop = FALSE, link.sort = "default", link.decreasing = TRUE,
    link.largest.ontop = TRUE, annotationTrack = "grid",
    annotationTrackHeight = circlize::mm_h(2),
    preAllocateTracks = list(track.height = label_track_height),
    big.gap = 14, small.gap = 0.15, reduce = -1
  )
  circlize::circos.trackPlotRegion(
    track.index = 1, bg.border = NA,
    panel.fun = function(x, y) {
      sector <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      ylim <- circlize::get.cell.meta.data("ylim")
      circlize::circos.text(
        x = mean(xlim), y = mean(ylim), labels = unname(sector_label[sector]),
        facing = "clockwise", niceFacing = TRUE, adj = c(0.5, 0.5), cex = label_cex
      )
    }
  )
  mtext(main, side = 3, line = 0.5, font = 2, cex = 1.25)
  text(-1.27, 0, labels = "Original", srt = 90, font = 2, cex = 1.15)
  text(1.27, 0, labels = "Assigned", srt = 270, font = 2, cex = 1.15)
}

run_step("Chord diagrams | Original to assigned pixel transitions", {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required for chord diagrams.")
  }
  require_file(saved_map_confusion_file, "5.1 normal-map confusion table")
  chord_dt <- fread(saved_map_confusion_file)
  require_columns(
    chord_dt, c("method", "original_zone", "predicted_zone", "n"),
    "normal_map_confusion_long_var.csv"
  )
  chord_dt <- chord_dt[
    chord_dt[["method"]] %in% method_order &
      chord_dt[["original_zone"]] %in% model_zoneID &
      chord_dt[["predicted_zone"]] %in% model_zoneID & chord_dt[["n"]] > 0
  ]
  zone_order <- as.character(model_zoneID)
  zone_colors <- zone_color_vector(model_zoneID)
  zone_labels <- setNames(as.character(model_zoneID), as.character(model_zoneID))
  if (is.null(category_lookup)) stop("Palette must contain category2 for category chord diagrams.")
  
  category_names <- c(
    cropland = "Cropland", orchard = "Orchard", plantation = "Tree\nplantation",
    forest = "Forest", grassland_meadow_steppe = "Grassland,\nmeadow & steppe",
    wetland = "Wetland", scrub = "Scrub", alpine_vegetation = "Alpine\nvegetation",
    desert = "Desert", no_vegetation = "No\nvegetation",
    rare_or_error = "Rare / error", novel = "Novel"
  )
  category_palette <- palette[zoneID %in% model_zoneID, {
    first_index <- which.min(zoneID)
    .(first_zone = zoneID[first_index], COLOR = COLOR[first_index])
  }, by = category2][order(first_zone)]
  category_palette[, display_label := unname(category_names[as.character(category2)])]
  category_palette[is.na(display_label), display_label := tools::toTitleCase(
    gsub("_", " ", as.character(category2), fixed = TRUE)
  )]
  category_order <- as.character(category_palette$category2)
  category_colors <- setNames(category_palette$COLOR, category_order)
  category_labels <- setNames(category_palette$display_label, category_order)
  zone_to_category <- setNames(as.character(category_lookup$category2),
                               as.character(category_lookup$zoneID))
  category_flow_list <- list()
  
  for (method_key in method_order) {
    method_dt <- chord_dt[chord_dt[["method"]] == method_key]
    zone_flow <- method_dt[, .(n = sum(n)), by = .(
      from = as.character(original_zone), to = as.character(predicted_zone)
    )]
    plot_chord(
      zone_flow, zone_order, zone_colors, zone_labels,
      file.path(chord_dir, paste0("normal_map_zone_chord_", method_key, ".pdf")),
      paste0(method_labels[method_key], ": Original zone to assigned zone"),
      label_cex = 0.62, label_track_height = 0.13
    )
    
    category_flow <- copy(zone_flow)
    category_flow[, from := unname(zone_to_category[from])]
    category_flow[, to := unname(zone_to_category[to])]
    category_flow <- category_flow[!is.na(from) & !is.na(to), .(n = sum(n)), by = .(from, to)]
    category_flow_list[[method_key]] <- data.table(
      method = method_key, original_category2 = category_flow$from,
      assigned_category2 = category_flow$to, n = category_flow$n
    )
    plot_chord(
      category_flow, category_order, category_colors, category_labels,
      file.path(chord_dir, paste0("normal_map_category_chord_", method_key, ".pdf")),
      paste0(method_labels[method_key], ": Original category to assigned category"),
      label_cex = 0.90, label_track_height = 0.18
    )
  }
  fwrite(rbindlist(category_flow_list, fill = TRUE),
         file.path(chord_dir, "normal_map_category_confusion_long_var.csv"))
  fwrite(category_palette, file.path(chord_dir, "category_chord_legend_var.csv"))
})

# 19. Figure 10 | Optional multiclass robustness =================================
run_step("Figure 10 | Optional multiclass robustness comparison", {
  require_file(multiclass_reference_map_file, "multiclass robustness reference map")
  robustness_methods <- c(method_order, "multiclass_rf")
  robustness_labels <- c(method_labels, multiclass_rf = "Multiclass RF")
  robustness <- build_reference_comparison(
    robustness_methods, "Common valid mask of final RF, MF and multiclass maps"
  )
  
  binary_f1 <- robustness$zone[
    robustness$zone[["method"]] %in% method_order,
    .(method, zone, binary_f1 = f1)
  ]
  multiclass_f1 <- robustness$zone[
    robustness$zone[["method"]] == "multiclass_rf",
    .(zone, multiclass_f1 = f1)
  ]
  compare_f1 <- merge(binary_f1, multiclass_f1, by = "zone", all.x = TRUE, sort = FALSE)
  ref <- rast(reference_file)[[1]]
  ref_area <- cellSize(ref, unit = "km")
  ref_modeled <- subst(ref, from = model_zoneID, to = model_zoneID, others = NA)
  area_dt <- as.data.table(zonal(ref_area, ref_modeled, fun = "sum", na.rm = TRUE))
  setnames(area_dt, names(area_dt)[1:2], c("zone", "area_km2"))
  area_dt[, `:=`(zone = as.integer(zone), area_km2 = as.numeric(area_km2))]
  compare_f1 <- merge(compare_f1, area_dt, by = "zone", all.x = TRUE, sort = FALSE)
  compare_f1[, `:=`(
    method_label = factor(robustness_labels[method], levels = robustness_labels[method_order]),
    zone_chr = as.character(zone)
  )]
  zone_colors <- zone_color_vector(model_zoneID)
  
  p10a <- ggplot(compare_f1, aes(
    x = binary_f1, y = multiclass_f1, colour = zone_chr, size = area_km2
  )) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.45) +
    geom_point(alpha = 0.80) +
    geom_text(aes(label = zone), check_overlap = TRUE, size = 2.2,
              nudge_y = 0.012, show.legend = FALSE) +
    facet_wrap(~method_label, nrow = 1) + scale_colour_manual(values = zone_colors, guide = "none") +
    scale_size_continuous(range = c(1.2, 6.0), name = expression("Original area (km"^2 * ")")) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(
      title = "Binary workflows vs multiclass reference-map F1",
      subtitle = "Above the dashed line favors multiclass; below the line favors the binary dual-niche workflow.",
      x = "Binary workflow reference-map F1", y = "Multiclass reference-map F1"
    ) + theme_ms(base_size = 10)
  save_gg(
    p10a,
    file.path(fig_dir, "Figure10a_all_binary_workflows_vs_multiclass_reference_map_F1_area_scaled.png"),
    8.8, 5.5
  )
  fwrite(compare_f1, file.path(tab_dir, "Figure10a_all_binary_workflows_vs_multiclass_reference_map_F1.csv"))
  
  overall <- copy(robustness$overall)
  overall[, method_label := robustness_labels[method]]
  metric_labels <- c(
    coverage = "Common comparison coverage", exact_zone_accuracy = "Exact-zone accuracy",
    broad_category_accuracy = "Broad-category accuracy",
    macro_balanced_accuracy = "Macro balanced accuracy", macro_recall = "Macro recall",
    macro_specificity = "Macro specificity", macro_precision = "Macro precision", macro_f1 = "Macro F1"
  )
  metrics_keep <- names(metric_labels)[names(metric_labels) %in% names(overall)]
  overall_long <- melt(
    overall, id.vars = c("method", "method_label"), measure.vars = metrics_keep,
    variable.name = "metric", value.name = "value"
  )
  overall_long[, `:=`(
    metric_label = factor(unname(metric_labels[metric]), levels = unname(metric_labels[metrics_keep])),
    method_label = factor(method_label, levels = rev(robustness_labels[robustness_methods]))
  )]
  p10b <- ggplot(overall_long, aes(x = value, y = method_label)) +
    geom_point(size = 2.8) + facet_wrap(~metric_label, ncol = 2) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    labs(
      title = "Overall reference-map reconstruction performance",
      subtitle = "All workflows are recalculated against the original vegetation raster on one common valid mask.",
      x = "Metric value", y = NULL
    ) + theme_ms(base_size = 10) + theme(panel.grid.major.y = element_blank())
  save_gg(p10b, file.path(fig_dir, "Figure10b_all_workflows_reference_map_overall_metrics.png"), 9.5, 8.5)
  fwrite(overall, file.path(tab_dir, "Figure10b_all_workflows_reference_map_overall_metrics.csv"))
}, optional = TRUE)

# 20. Save log ==================================================================
step_log_dt <- rbindlist(step_log, fill = TRUE)
log_file <- file.path(vis_dir, "visualization_step_log.csv")
fwrite(step_log_dt, log_file)
cat0("[SAVED] ", log_file)
failed_steps <- step_log_dt[status == "error"]
if (strict_mode && nrow(failed_steps) > 0) {
  cat0("\nVISUALIZATION INCOMPLETE")
  print(failed_steps[, .(step, message)])
  stop("One or more required visualization steps failed. See: ", log_file, call. = FALSE)
}
cat0("\nCOMPLETE")
cat0("Main methods: ", paste(method_labels[method_order], collapse = ", "))
cat0("Suitability threshold: >= ", novel_threshold)
cat0("Figure folder: ", fig_dir)
cat0("Chord folder: ", chord_dir)
cat0("Step log: ", log_file)
