# 11. visualization.R
# ============================================================
# Figure 1 only.
#
# Purpose:
#   Rebuild Table 1 support files and Figure 1 from existing results.
#   This script is intentionally limited to Figure 1, so "Run All"
#   will not regenerate Figures 2-10.
#
# Required input files:
#   accuracy_climate/climate_rf_accuracy_summary.csv
#   accuracy_soil/soil_rf_accuracy_summary.csv
#   assessment/rf_test_zone_metrics.csv
#
# Key outputs:
#   visualization/tables/Table1_*.csv
#   visualization/figures/Figure1a_onehot_binary_RF_summary.png
#   visualization/figures/Figure1b_onehot_binary_RF_zone_level_metrics.png
#   visualization/figures/Figure1_onehot_binary_RF_OOB_train_test_performance_dotrange.png
#   visualization/visualization_step_log.csv
# ============================================================

library(data.table)
library(ggplot2)
library(grid)

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

# Modeled zones. Zone 8 and zone 51 were not modeled.
model_zoneID <- c(1:7, 9:50, 52:55)

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

# Non-zone colors.
niche_cols <- c(
  climate = "#2E5E8C",
  soil = "#9A7B32"
)

fig_dpi <- 320

# 1. General helpers ============================================================

cat0 <- function(...) {
  cat(..., "\n", sep = "")
}

norm_name <- function(x) {
  tolower(gsub("[^a-z0-9]+", "_", x))
}

relative_path <- function(x) {
  base_norm <- normalizePath(base_dir, winslash = "/", mustWork = FALSE)
  x_norm <- normalizePath(x, winslash = "/", mustWork = FALSE)
  gsub(paste0("^", base_norm, "/?"), "", x_norm)
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

find_files <- function(pattern, root = base_dir) {
  if (!dir.exists(root)) return(character(0))
  list.files(
    root,
    pattern = pattern,
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
}

pick_file <- function(files, prefer = character()) {
  files <- unique(files)
  files <- files[file.exists(files)]
  if (length(files) == 0) return(NA_character_)
  score <- rep(0L, length(files))
  for (p in prefer) {
    score <- score + as.integer(grepl(p, files, ignore.case = TRUE, fixed = TRUE))
  }
  files[order(-score, nchar(files), files)][1]
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

infer_workflow <- function(x) {
  x <- tolower(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("optimized.*mf|optimized_mf", x)] <- "optimized_mf_rf"
  out[grepl("optimized.*rf|optimized_rf", x) & is.na(out)] <- "optimized_rf"
  out[grepl("plain.*mf|plain_mf", x) & is.na(out)] <- "plain_mf_rf"
  out[grepl("plain.*rf|plain_rf", x) & is.na(out)] <- "plain_rf"
  out
}

metric_from_col <- function(x) {
  x0 <- norm_name(x)
  
  # OOB metrics.
  if (grepl("oob.*acc|acc.*oob|oob.*accuracy|accuracy.*oob", x0)) return("OOB accuracy")
  if (grepl("oob.*auc|auc.*oob", x0)) return("OOB AUC")
  if (grepl("oob.*f1|f1.*oob", x0)) return("OOB F1")
  if (grepl("oob.*precision|precision.*oob", x0)) return("OOB precision")
  if (grepl("oob.*recall|recall.*oob|oob.*sensitivity|sensitivity.*oob", x0)) return("OOB recall")
  if (grepl("oob.*specificity|specificity.*oob", x0)) return("OOB specificity")
  
  # Training metrics.
  if (grepl("train.*acc|acc.*train|train.*accuracy|accuracy.*train", x0)) return("Train accuracy")
  if (grepl("train.*auc|auc.*train", x0)) return("Train AUC")
  if (grepl("train.*f1|f1.*train", x0)) return("Train F1")
  if (grepl("train.*precision|precision.*train", x0)) return("Train precision")
  if (grepl("train.*recall|recall.*train|train.*sensitivity|sensitivity.*train", x0)) return("Train recall")
  if (grepl("train.*specificity|specificity.*train", x0)) return("Train specificity")
  
  NA_character_
}

default_zone_colors <- function(vals) {
  vals <- sort(unique(as.integer(vals)))
  vals <- vals[!is.na(vals)]
  cols <- grDevices::hcl.colors(length(vals), palette = "Dark 3")
  names(cols) <- as.character(vals)
  cols
}

zone_color_vector <- function(vals = model_zoneID) {
  cols <- default_zone_colors(vals)
  cols[as.character(vals)]
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
      panel.spacing.x = unit(0.42, "lines"),
      panel.spacing.y = unit(0.42, "lines"),
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
      cat0("[SKIP/ERROR] ", step_name)
      cat0("  ", conditionMessage(e))
    }
  )
  invisible(NULL)
}

# 2. Binary RF metric builders ==================================================

find_binary_summary_file <- function(niche_type) {
  if (niche_type == "climate") {
    exact <- c(
      file.path(base_dir, "accuracy_climate", "climate_rf_accuracy_summary.csv"),
      file.path(base_dir, "climate_rf_accuracy_summary.csv")
    )
    prefer <- c("accuracy_climate", "climate_rf_accuracy_summary")
    pattern <- "(clim|climate).*rf.*(accuracy|summary|metric).*\\.csv$"
  } else {
    exact <- c(
      file.path(base_dir, "accuracy_soil", "soil_rf_accuracy_summary.csv"),
      file.path(base_dir, "soil_rf_accuracy_summary.csv")
    )
    prefer <- c("accuracy_soil", "soil_rf_accuracy_summary")
    pattern <- "soil.*rf.*(accuracy|summary|metric).*\\.csv$"
  }
  
  exact <- exact[file.exists(exact)]
  if (length(exact) > 0) return(exact[1])
  
  files <- find_files(pattern)
  files <- files[
    !grepl(
      "(normal_map|confusion|assigned|future tree niche|population|species|ranking|visualization|area|lookup|index|multiclass|multi_class|robust)",
      files,
      ignore.case = TRUE
    )
  ]
  
  pick_file(files, prefer = prefer)
}

read_train_oob_metrics <- function(niche_type, file) {
  if (is.na(file) || !file.exists(file)) {
    stop("Missing binary RF summary file for ", niche_type)
  }
  
  dt <- fread(file)
  
  z_col <- pick_col(
    dt,
    c("zoneID", "zone_id", "zone", "ecotype", "class", "veg_zone")
  )
  if (is.na(z_col)) stop("Cannot find zone column in: ", file)
  
  workflow_col <- pick_col(
    dt,
    c("workflow", "method", "model", "model_version", "version", "model_type", "rf_type")
  )
  
  dt[, zoneID_tmp := extract_zone_id(get(z_col))]
  dt <- dt[zoneID_tmp %in% model_zoneID]
  if (nrow(dt) == 0) stop("No modeled zones found in: ", file)
  
  if (!is.na(workflow_col)) {
    dt[, workflow := infer_workflow(get(workflow_col))]
    dt[, workflow_source := workflow_col]
  } else {
    # Used only when a summary file stores exactly four rows per zone but no method column.
    dt[, row_in_zone := seq_len(.N), by = zoneID_tmp]
    dt[, workflow := workflow_order[((row_in_zone - 1L) %% length(workflow_order)) + 1L]]
    dt[, workflow_source := "inferred_from_row_order"]
  }
  
  numeric_cols <- names(dt)[sapply(dt, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, c("zoneID_tmp", "row_in_zone"))
  
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
  
  if (length(out) == 0) {
    stop("No usable OOB/Train metric columns found in: ", file)
  }
  
  out <- rbindlist(out, fill = TRUE)
  out <- out[
    workflow %in% workflow_order &
      !is.na(value) &
      is.finite(value)
  ]
  
  if (nrow(out) == 0) {
    stop("No valid OOB/Train workflow metrics extracted from: ", file)
  }
  
  out[]
}

find_test_metrics_file <- function() {
  exact <- c(
    file.path(base_dir, "assessment", "rf_test_zone_metrics.csv"),
    file.path(base_dir, "rf_test_zone_metrics.csv")
  )
  
  exact <- exact[file.exists(exact)]
  if (length(exact) > 0) return(exact[1])
  
  files <- find_files("rf_test_zone_metrics.*\\.csv$")
  files <- files[!grepl("visualization", files, ignore.case = TRUE)]
  
  pick_file(files, prefer = c("assessment", "rf_test_zone_metrics"))
}

read_test_metrics <- function(file) {
  if (is.na(file) || !file.exists(file)) {
    cat0("[WARNING] Missing rf_test_zone_metrics.csv. Figure 1 will use OOB/Train metrics only.")
    return(data.table())
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
  
  if (nrow(dt) == 0) {
    cat0("[WARNING] No valid rows found in rf_test_zone_metrics.csv.")
    return(data.table())
  }
  
  test_metric_cols <- c(
    accuracy = "Test accuracy",
    balanced_accuracy = "Test balanced accuracy",
    auc = "Test AUC",
    f1 = "Test F1",
    precision = "Test precision",
    recall = "Test recall",
    specificity = "Test specificity",
    tss = "Test TSS"
  )
  
  test_metric_cols <- test_metric_cols[names(test_metric_cols) %in% names(dt)]
  
  if (length(test_metric_cols) == 0) {
    cat0("[WARNING] No supported test metric columns found in rf_test_zone_metrics.csv.")
    return(data.table())
  }
  
  out <- list()
  
  for (mc in names(test_metric_cols)) {
    tmp <- data.table(
      niche_type = dt$niche_type,
      zoneID = dt$zoneID,
      workflow = dt$workflow,
      metric = unname(test_metric_cols[mc]),
      value = as.numeric(dt[[mc]]),
      metric_raw = mc,
      workflow_source = "assessment/rf_test_zone_metrics.csv",
      source_file = relative_path(file)
    )
    
    out[[length(out) + 1L]] <- tmp
  }
  
  out <- rbindlist(out, fill = TRUE)
  out <- out[!is.na(value) & is.finite(value)]
  
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
  
  long <- long[
    ,
    .(
      value = mean(value, na.rm = TRUE),
      workflow_source = paste(unique(workflow_source), collapse = "; "),
      metric_raw = paste(unique(metric_raw), collapse = "; "),
      source_file = paste(unique(source_file), collapse = "; ")
    ),
    by = .(niche_type, workflow, zoneID, metric)
  ]
  
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
  
  # Most Figure 1 metrics are bounded by 0 and 1. TSS can be negative.
  summary_long[!grepl("TSS", metric), `:=`(
    ci95_low = pmax(ci95_low, 0),
    ci95_high = pmin(ci95_high, 1)
  )]
  
  summary_long[, workflow_label := workflow_labels[as.character(workflow)]]
  
  setorder(summary_long, niche_type, workflow, metric)
  setorder(long, niche_type, workflow, zoneID, metric)
  
  list(long = long, summary_long = summary_long)
}

# 3. Table 1 support files + Figure 1 ==========================================

run_step("Figure 1 | One-hot binary RF summary and zone-level panels", {
  tbl <- make_binary_rf_table1()
  long <- tbl$long
  summary_long <- tbl$summary_long
  
  long_file <- file.path(tab_dir, "Table1_binary_RF_zone_level_metrics_long.csv")
  fwrite(long, long_file)
  cat0("[SAVED] ", long_file)
  
  zone_wide <- dcast(
    long,
    niche_type + workflow + zoneID ~ metric,
    value.var = "value",
    fun.aggregate = mean
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
  
  # -----------------------------------------------------------------------------
  # Shared plotting data preparation
  # -----------------------------------------------------------------------------
  
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
  
  if (nrow(plot_summary) == 0) {
    stop("No OOB/Train/Test metrics found.")
  }
  
  plot_summary[, eval_set := sub(" .*", "", metric)]
  plot_summary[, metric_name := sub("^(OOB|Train|Test) ", "", metric)]
  plot_summary[, metric_label := unname(metric_label_map[metric_name])]
  plot_summary <- plot_summary[!is.na(metric_label)]
  
  if (nrow(plot_summary) == 0) {
    stop("No supported Figure 1 metrics found.")
  }
  
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
  
  x_lower <- min(x_all, na.rm = TRUE)
  x_upper <- max(x_all, na.rm = TRUE)
  x_lower <- min(0, floor(x_lower / 0.1) * 0.1)
  x_upper <- max(1, ceiling(x_upper / 0.1) * 0.1)
  if (x_upper <= x_lower) x_upper <- x_lower + 1
  x_breaks <- pretty(c(x_lower, x_upper), n = 6)
  
  # -----------------------------------------------------------------------------
  # Figure 1a: workflow-level summary
  # -----------------------------------------------------------------------------
  
  p1 <- ggplot(plot_summary, aes(x = mean, y = metric_label)) +
    geom_segment(
      aes(x = ci95_low, xend = ci95_high, yend = metric_label),
      linewidth = 0.9,
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
      limits = c(x_lower, x_upper),
      breaks = x_breaks,
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = "Figure 1a. One-hot binary RF performance across vegetation zones",
      subtitle = "Points are means across modeled vegetation zones; horizontal lines show approximate 95% confidence intervals. OOB/Train metrics come from RF training summaries, and Test metrics come from independent balanced test-set assessment.",
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
  
  fig1a_file <- file.path(fig_dir, "Figure1a_onehot_binary_RF_summary.png")
  save_gg(p1, fig1a_file, width = 13.5, height = 9.6)
  
  # Keep the original Figure 1 filename for backward compatibility.
  save_gg(
    p1,
    file.path(fig_dir, "Figure1_onehot_binary_RF_OOB_train_test_performance_dotrange.png"),
    width = 13.5,
    height = 9.6
  )
  
  # -----------------------------------------------------------------------------
  # Figure 1b: zone-level variation
  # -----------------------------------------------------------------------------
  
  p2 <- ggplot(plot_zone, aes(x = value, y = metric_label)) +
    geom_point(
      aes(colour = zoneID_chr),
      position = position_jitter(width = 0, height = 0.13),
      size = 1.1,
      alpha = 0.60,
      stroke = 0,
      show.legend = FALSE
    ) +
    geom_segment(
      data = plot_summary,
      aes(x = ci95_low, xend = ci95_high, y = metric_label, yend = metric_label),
      inherit.aes = FALSE,
      linewidth = 0.8,
      colour = "grey15",
      alpha = 0.95
    ) +
    geom_point(
      data = plot_summary,
      aes(x = mean, y = metric_label, fill = niche_type),
      inherit.aes = FALSE,
      shape = 23,
      size = 2.35,
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
      limits = c(x_lower, x_upper),
      breaks = x_breaks,
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = "Figure 1b. Zone-level variation in one-hot binary RF performance",
      subtitle = "Small points are individual vegetation zones. Diamonds and horizontal lines show the corresponding mean and approximate 95% confidence interval across zones.",
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
  
  fig1b_file <- file.path(fig_dir, "Figure1b_onehot_binary_RF_zone_level_metrics.png")
  save_gg(p2, fig1b_file, width = 13.5, height = 10.4)
})

# 4. Save step log ============================================================== 

cat0("\n============================================================")
cat0("SAVE VISUALIZATION STEP LOG")
cat0("============================================================")

step_log_dt <- rbindlist(step_log, fill = TRUE)
log_file <- file.path(vis_dir, "visualization_step_log.csv")
fwrite(step_log_dt, log_file)

cat0("[SAVED] ", log_file)
cat0("\nCOMPLETE")
cat0("Only Figure 1 was regenerated.")
cat0("Figure folder: ", fig_dir)
cat0("Table folder: ", tab_dir)
cat0("Step log: ", log_file)
