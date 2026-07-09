# 11. visualization.R
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
#   - Each table/figure is an independent step.
#   - Missing files or incompatible columns are skipped without stopping.
#   - Zone-related figures use color_palette_China.csv whenever possible.
#   - Figure 1 has two panels:
#       Figure 1a = mean performance across modeled zones.
#       Figure 1b = zone-level performance values.
#   - Figure 6 now has Figure 6a (area change) and Figure 6b (zone-transition sankey).
#   - Figures 8 and 9 plot all tree species, 10 species per page (2 rows x 5 columns), for all scenarios.
#   - Figure 10 point transparency is 10% (alpha = 0.90); zone labels are fully opaque.
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

# Modeled zones. Zone 8 and zone 51 were not modeled.
model_zoneID <- c(1:7, 9:50, 52:55)
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

scenario_order <- c("normal", future_order)

fig_dpi <- 320
display_max_cells <- 300000
species_page_ncol <- 5
species_page_nrow <- 2

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
      cat0("[SKIP/ERROR] ", step_name)
      cat0("  ", conditionMessage(e))
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
  vals <- sort(unique(c(model_zoneID, 8, novel_value)))
  cols <- default_zone_colors(vals)
  lookup <- data.table(
    zoneID = vals,
    zone_label = ifelse(vals == novel_value, "Novel ecotype", paste0("Zone ", vals)),
    color = unname(cols[as.character(vals)])
  )
  
  palette_files <- unique(c(
    file.path(base_dir, "color_palette_China.csv"),
    find_files("(color_palette|palette).*\\.csv$")
  ))
  pal_file <- pick_file(palette_files, prefer = c("color_palette_China.csv", "palette"))
  
  if (!is.na(pal_file)) {
    pal <- tryCatch(fread(pal_file), error = function(e) NULL)
    if (!is.null(pal)) {
      z_col <- pick_col(pal, c("zoneID", "zone_id", "zoneid", "zone", "id", "value"))
      c_col <- pick_col(pal, c("COLOR", "color", "colour", "hex", "hex_color"))
      n_col <- pick_col(pal, c("zone_name", "zone_label", "name", "vegetation", "ecosystem", "type"))
      
      if (!is.na(z_col)) {
        tmp <- data.table(zoneID = as.integer(pal[[z_col]]))
        
        if (!is.na(c_col)) {
          tmp[, color_new := as.character(pal[[c_col]])]
          tmp[!grepl("^#", color_new) & !is.na(color_new), color_new := paste0("#", color_new)]
        }
        
        if (!is.na(n_col)) {
          tmp[, label_new := as.character(pal[[n_col]])]
        }
        
        lookup <- merge(lookup, tmp, by = "zoneID", all.x = TRUE)
        
        if ("color_new" %in% names(lookup)) {
          lookup[!is.na(color_new), color := color_new]
          lookup[, color_new := NULL]
        }
        
        if ("label_new" %in% names(lookup)) {
          lookup[!is.na(label_new), zone_label := label_new]
          lookup[, label_new := NULL]
        }
      }
    }
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
  
  z_col <- pick_col(dt, c("zoneID", "zone_id", "zone", "ecotype", "class", "veg_zone"))
  if (is.na(z_col)) stop("Cannot find zone column in: ", file)
  
  workflow_col <- pick_col(dt, c("workflow", "method", "model", "model_version", "version", "model_type", "rf_type"))
  
  dt[, zoneID_tmp := extract_zone_id(get(z_col))]
  dt <- dt[zoneID_tmp %in% model_zoneID]
  if (nrow(dt) == 0) stop("No modeled zones found in: ", file)
  
  if (!is.na(workflow_col)) {
    dt[, workflow := infer_workflow(get(workflow_col))]
    dt[, workflow_source := workflow_col]
  } else {
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
  
  if (length(out) == 0) stop("No usable OOB/Train metric columns found in: ", file)
  
  out <- rbindlist(out, fill = TRUE)
  out <- out[workflow %in% workflow_order & !is.na(value) & is.finite(value)]
  
  if (nrow(out) == 0) stop("No valid OOB/Train workflow metrics extracted from: ", file)
  
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
  
  summary_long[!grepl("TSS", metric), `:=`(
    ci95_low = pmax(ci95_low, 0),
    ci95_high = pmin(ci95_high, 1)
  )]
  
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

raster_to_plot_dt <- function(file, categorical = TRUE, population_mode = FALSE) {
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
      dt <- dt[value %in% c(model_zoneID, novel_value)]
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
    plot_title = NULL) {
  
  plot_list <- list()
  
  for (i in seq_along(files)) {
    cat0("  Preparing raster panel: ", titles[i])
    
    dt <- tryCatch(
      raster_to_plot_dt(
        files[[i]],
        categorical = categorical,
        population_mode = population_mode
      ),
      error = function(e) {
        cat0("  [SKIP PANEL] ", titles[i], " | ", conditionMessage(e))
        data.table()
      }
    )
    
    if (nrow(dt) == 0) next
    dt[, panel := titles[i]]
    plot_list[[length(plot_list) + 1L]] <- dt
  }
  
  if (length(plot_list) == 0) stop("No valid raster cells were available for plotting.")
  
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
    "Reference",
    paste0(sub("SSP[0-9]+$", "", s), "\n", sub("^.*(SSP[0-9]+)$", "\\1", s))
  )
}

scenario_period <- function(s) {
  ifelse(s == "normal", "Reference", sub("SSP[0-9]+$", "", s))
}

scenario_ssp <- function(s) {
  ifelse(s == "normal", "Reference", sub("^.*(SSP[0-9]+)$", "\\1", s))
}

find_assigned_map <- function(method, scenario) {
  roots <- unique(c(file.path(result_map_root, method), result_map_root))
  roots <- roots[dir.exists(roots)]
  
  files <- unique(unlist(lapply(
    roots,
    function(root) {
      list.files(
        root,
        pattern = "\\.tif$",
        recursive = TRUE,
        full.names = TRUE,
        ignore.case = TRUE
      )
    }
  )))
  
  files <- files[
    grepl("assigned_zone", basename(files), ignore.case = TRUE) &
      grepl(method, files, ignore.case = TRUE, fixed = TRUE) &
      grepl(scenario, basename(files), ignore.case = TRUE, fixed = TRUE)
  ]
  
  files <- files[!grepl("(color|rgb|plot|legend)", basename(files), ignore.case = TRUE)]
  
  pick_file(
    files,
    prefer = c("threshold0.2", "threshold0.1", "maskNA8", "noNovelNormal", method, scenario)
  )
}

find_dual_file <- function(method, scenario, zone) {
  f <- file.path(
    dual_root,
    method,
    scenario,
    paste0("dual_suitability_zone", zone, ".tif")
  )
  
  if (file.exists(f)) return(f)
  
  files <- find_files(paste0("dual_suitability_zone", zone, "\\.tif$"), root = dual_root)
  files <- files[
    grepl(method, files, ignore.case = TRUE, fixed = TRUE) &
      grepl(scenario, files, ignore.case = TRUE, fixed = TRUE)
  ]
  
  pick_file(files, prefer = c(method, scenario, paste0("zone", zone)))
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
  area_by_zone(reference_file)[zoneID %in% model_zoneID]
}

get_species_population_table <- function() {
  files <- unique(c(
    file.path(base_dir, "future tree niche", "tables", "population_projection_lookup.csv"),
    file.path(base_dir, "future tree niche dual suitability", "tables", "population_projection_lookup.csv"),
    find_files("species_zone_population_long.*\\.csv$")
  ))
  
  files <- files[file.exists(files)]
  files <- files[!grepl("visualization", files, ignore.case = TRUE)]
  f <- pick_file(files, prefer = c("population_projection_lookup", "species_zone_population_long"))
  
  if (is.na(f)) stop("No species-population lookup file found.")
  
  dt <- fread(f)
  
  species_col <- pick_col(dt, c("species", "species_name", "taxon", "tree_species"))
  zone_col <- pick_col(dt, c("source_zone", "population_zone", "zoneID", "zone_id", "zone", "ecotype"))
  
  if (is.na(species_col) || is.na(zone_col)) {
    stop("Cannot identify species/source-zone columns in: ", f)
  }
  
  out <- unique(data.table(
    species = as.character(dt[[species_col]]),
    source_zone = extract_zone_id(dt[[zone_col]])
  ))
  
  out <- out[!is.na(species) & source_zone %in% model_zoneID]
  if (nrow(out) == 0) stop("No valid species-source-zone records found in: ", f)
  
  out[]
}

dual_to_points <- function(file, source_zone, threshold = dual_plot_threshold) {
  if (is.na(file) || !file.exists(file)) return(data.table())
  
  x <- rast(file)[[1]]
  x <- thin_raster_for_plot(x, categorical = FALSE, max_cells = display_max_cells)
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

build_dual_cache <- function(method, scenario, zones) {
  zones <- sort(unique(as.integer(zones)))
  cache <- list()
  
  for (z in zones) {
    f <- find_dual_file(method, scenario, z)
    
    if (is.na(f) || !file.exists(f)) {
      cat0("  [MISSING DUAL] ", method, " | ", scenario, " | zone ", z)
      next
    }
    
    cache[[as.character(z)]] <- dual_to_points(f, z, threshold = dual_plot_threshold)
  }
  
  cache
}

plot_species_population_dual_pages <- function(pop_tbl, method, scenario, outfile_prefix) {
  species_order <- sort(unique(pop_tbl$species))
  zones_needed <- sort(unique(pop_tbl$source_zone))
  cache <- build_dual_cache(method, scenario, zones_needed)
  
  if (length(cache) == 0) stop("No dual suitability rasters could be read for ", method, " | ", scenario)
  
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
      
      if (length(dts) == 0) next
      
      dt <- rbindlist(dts, fill = TRUE)
      if (nrow(dt) == 0) next
      
      dt[, species := sp]
      plot_list[[length(plot_list) + 1L]] <- dt
    }
    
    if (length(plot_list) == 0) next
    
    plot_dt <- rbindlist(plot_list, fill = TRUE)
    plot_dt[, species := factor(species, levels = sp_sel)]
    plot_dt[, source_zone_chr := factor(as.character(source_zone), levels = as.character(zones_needed))]
    
    p <- ggplot(plot_dt, aes(x = x, y = y, colour = source_zone_chr)) +
      geom_point(size = 0.075, alpha = 0.30) +
      facet_wrap(~ species, ncol = species_page_ncol) +
      coord_equal(expand = FALSE) +
      scale_colour_manual(values = cols, drop = FALSE) +
      labs(
        title = paste0("Tree population dual-niche projections: ", method, " | ", scenario),
        subtitle = paste0("Each population niche is drawn from the dual niche of its corresponding source zone. Cells with dual suitability > ", dual_plot_threshold, " are colored by source zone with 30% opacity."),
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
      guides(colour = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(size = 2, alpha = 1)))
    
    outfile <- paste0(
      outfile_prefix,
      "_",
      safe_name(method),
      "_",
      safe_name(scenario),
      "_page",
      sprintf("%02d", pg),
      ".png"
    )
    
    ggsave(outfile, p, width = 15.5, height = 7.0, dpi = fig_dpi, bg = "white")
    cat0("[SAVED] ", outfile)
  }
}

plot_species_level_dual_pages <- function(pop_tbl, method, scenario, outfile_prefix) {
  species_order <- sort(unique(pop_tbl$species))
  zones_needed <- sort(unique(pop_tbl$source_zone))
  cache <- build_dual_cache(method, scenario, zones_needed)
  
  if (length(cache) == 0) stop("No dual suitability rasters could be read for ", method, " | ", scenario)
  
  page_size <- species_page_ncol * species_page_nrow
  page_id <- ceiling(seq_along(species_order) / page_size)
  
  for (pg in sort(unique(page_id))) {
    sp_sel <- species_order[page_id == pg]
    plot_list <- list()
    
    for (sp in sp_sel) {
      z_sp <- pop_tbl[species == sp, sort(unique(source_zone))]
      dts <- cache[as.character(z_sp)]
      dts <- dts[lengths(dts) > 0]
      
      if (length(dts) == 0) next
      
      dt <- rbindlist(dts, fill = TRUE)
      if (nrow(dt) == 0) next
      
      dt <- dt[
        ,
        .(dual_suitability = max(dual_suitability, na.rm = TRUE)),
        by = .(x, y)
      ]
      
      dt[, species := sp]
      plot_list[[length(plot_list) + 1L]] <- dt
    }
    
    if (length(plot_list) == 0) next
    
    plot_dt <- rbindlist(plot_list, fill = TRUE)
    plot_dt[, species := factor(species, levels = sp_sel)]
    
    p <- ggplot(plot_dt, aes(x = x, y = y, colour = dual_suitability)) +
      geom_point(size = 0.075, alpha = 0.30) +
      facet_wrap(~ species, ncol = species_page_ncol) +
      coord_equal(expand = FALSE) +
      scale_colour_gradientn(
        colours = grDevices::hcl.colors(60, palette = "YlGnBu"),
        name = "Max dual\nsuitability"
      ) +
      labs(
        title = paste0("Species-level dual-niche projections: ", method, " | ", scenario),
        subtitle = paste0("Each species map is the pixel-wise maximum across its population dual niches; cells with dual suitability > ", dual_plot_threshold, " are shown with 30% opacity.")
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
      outfile_prefix,
      "_",
      safe_name(method),
      "_",
      safe_name(scenario),
      "_page",
      sprintf("%02d", pg),
      ".png"
    )
    
    ggsave(outfile, p, width = 15.5, height = 7.0, dpi = fig_dpi, bg = "white")
    cat0("[SAVED] ", outfile)
  }
}

# 6. Table 1 + Figure 1 =========================================================

run_step("Table 1 and Figure 1 | One-hot climate/soil binary RF performance", {
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
  
  # Zoomed range for readability. Current metrics are not expected below 0.4.
  x_lower <- 0.4
  x_upper <- max(1, ceiling(max(x_all, na.rm = TRUE) / 0.05) * 0.05)
  if (x_upper <= x_lower) x_upper <- 1
  x_breaks <- seq(x_lower, x_upper, by = 0.1)
  
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
      limits = c(x_lower, x_upper),
      breaks = x_breaks,
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = "Figure 1a. One-hot binary RF performance across vegetation zones",
      subtitle = "Points are means across modeled vegetation zones; horizontal lines show approximate 95% confidence intervals. OOB/Train metrics come from RF training summaries, and Test metrics come from independent balanced test-set assessment. The x-axis is zoomed to the observed performance range for readability.",
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
      limits = c(x_lower, x_upper),
      breaks = x_breaks,
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = "Figure 1b. Zone-level variation in one-hot binary RF performance",
      subtitle = "Small points are individual vegetation zones. Diamonds and horizontal lines show the corresponding mean and approximate 95% confidence interval across zones. The x-axis is zoomed to the observed performance range for readability.",
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

run_step("Figure 2 | Original vs predicted reference maps for all workflows", {
  if (!file.exists(reference_file)) stop("Missing reference raster: ", reference_file)
  
  figure_letters <- letters[seq_along(map_method_order)]
  
  for (i in seq_along(map_method_order)) {
    m <- map_method_order[i]
    pred_file <- find_assigned_map(m, "normal")
    
    if (is.na(pred_file)) {
      cat0("  [SKIP] missing predicted reference map for method: ", m)
      next
    }
    
    plot_zone_panel_gg(
      files = list(reference_file, pred_file),
      titles = c("Original vegetation map", "Predicted reference map"),
      outfile = file.path(
        fig_dir,
        paste0("Figure2", figure_letters[i], "_original_vs_predicted_reference_map_", m, ".png")
      ),
      categorical = TRUE,
      show_legend = FALSE,
      ncol = 2,
      width = 10,
      height = 4.8,
      plot_title = paste0("Reference-map reconstruction: ", map_method_labels[m])
    )
  }
})

# 8. Table 2 ===================================================================

run_step("Table 2 | Reference-map accuracy for four binary workflows", {
  files <- find_files("normal_map_overall_metrics.*\\.csv$")
  files <- files[!grepl("visualization", files, ignore.case = TRUE)]
  f <- pick_file(files, prefer = c("normal_map_overall_metrics"))
  
  if (is.na(f)) stop("No normal_map_overall_metrics.csv found.")
  
  dt <- fread(f)
  dt[, source_file := relative_path(f)]
  
  out_file <- file.path(tab_dir, "Table2_reference_map_accuracy_binary_workflows.csv")
  fwrite(dt, out_file)
  cat0("[SAVED] ", out_file)
})

# 9. Figure 3 ==================================================================

run_step("Figure 3 | Climate and soil binary RF zone-level metrics for all workflows", {
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
      cat0("  [SKIP] no records for ", nt)
      next
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
  files <- find_files("normal_map_confusion_long.*\\.csv$")
  files <- files[!grepl("visualization", files, ignore.case = TRUE)]
  f <- pick_file(files, prefer = c(preferred_method, "normal_map_confusion_long"))
  
  if (is.na(f)) stop("No normal_map_confusion_long.csv found.")
  
  dt <- fread(f)
  
  from_col <- pick_col(dt, c("original_zone", "reference_zone", "actual_zone", "true_zone", "from", "truth", "observed", "ori"))
  to_col <- pick_col(dt, c("assigned_zone", "predicted_zone", "pred_zone", "prediction", "to", "pred"))
  count_col <- pick_col(dt, c("n", "count", "freq", "frequency", "pixels", "pixel_count"))
  
  if (is.na(from_col) || is.na(to_col) || is.na(count_col)) {
    stop("Cannot identify from/to/count columns in confusion long table.")
  }
  
  flow <- data.table(
    from = extract_zone_id(dt[[from_col]]),
    to = extract_zone_id(dt[[to_col]]),
    count = as.numeric(dt[[count_col]])
  )
  
  if ("method" %in% names(dt)) {
    flow[, method := normalise_map_method(dt$method)]
    flow <- flow[is.na(method) | method == preferred_method]
  }
  
  flow <- flow[
    from %in% model_zoneID &
      to %in% model_zoneID &
      from != to &
      !is.na(count) &
      is.finite(count) &
      count > 0
  ]
  
  if (nrow(flow) == 0) stop("No off-diagonal confusion flows found.")
  
  flow <- flow[, .(count = sum(count, na.rm = TRUE)), by = .(from, to)][order(-count)]
  top_n <- min(20L, nrow(flow))
  flow_top <- flow[seq_len(top_n)]
  
  flow_top[, from_chr := as.character(from)]
  flow_top[, to_chr := as.character(to)]
  flow_top[, flow_label := paste0("Zone ", from, " -> Zone ", to)]
  flow_top[, flow_label := factor(flow_label, levels = rev(flow_label))]
  
  cols <- zone_color_vector(unique(flow_top$from))
  
  out_file <- file.path(fig_dir, "Figure4_major_ecotype_confusion_flows.png")
  
  if (requireNamespace("ggalluvial", quietly = TRUE)) {
    p <- ggplot(flow_top, aes(y = count, axis1 = from_chr, axis2 = to_chr)) +
      ggalluvial::geom_alluvium(
        aes(fill = from_chr),
        width = 1 / 12,
        alpha = 0.78,
        show.legend = FALSE
      ) +
      ggalluvial::geom_stratum(width = 1 / 8, fill = "grey94", colour = "grey40") +
      ggalluvial::geom_text(
        stat = "stratum",
        aes(label = paste0("Zone ", after_stat(stratum))),
        size = 2.6
      ) +
      scale_fill_manual(values = cols) +
      scale_x_discrete(limits = c("Original", "Predicted"), expand = c(0.08, 0.08)) +
      labs(
        title = paste0("Major reference-map confusion flows: ", preferred_method),
        subtitle = "Top off-diagonal original-to-predicted ecotype transitions.",
        x = NULL,
        y = "Pixel count"
      ) +
      theme_ms()
  } else {
    p <- ggplot(flow_top, aes(x = count, y = flow_label, fill = from_chr)) +
      geom_col(width = 0.72, alpha = 0.80) +
      scale_fill_manual(values = cols, guide = "none") +
      labs(
        title = paste0("Major reference-map confusion flows: ", preferred_method),
        subtitle = "Top off-diagonal original-to-predicted ecotype transitions.",
        x = "Pixel count",
        y = NULL
      ) +
      theme_ms(base_size = 10) +
      theme(panel.grid.major.y = element_blank())
  }
  
  save_gg(p, out_file, width = 8.8, height = 6.6)
})

# 11. Figure 5 =================================================================

run_step("Figure 5 | Ecosystem niche maps for all workflows and periods", {
  if (!file.exists(reference_file)) stop("Missing reference raster: ", reference_file)
  
  panel_scenarios <- c(
    "normal",
    "2011-2040SSP245",
    "2041-2070SSP245",
    "2071-2100SSP245",
    "normal_pred",
    "2011-2040SSP585",
    "2041-2070SSP585",
    "2071-2100SSP585"
  )
  
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
  
  figure_letters <- letters[seq_along(map_method_order)]
  index_list <- list()
  
  for (i in seq_along(map_method_order)) {
    m <- map_method_order[i]
    
    panel_files <- list(
      reference_file,
      find_assigned_map(m, "2011-2040SSP245"),
      find_assigned_map(m, "2041-2070SSP245"),
      find_assigned_map(m, "2071-2100SSP245"),
      find_assigned_map(m, "normal"),
      find_assigned_map(m, "2011-2040SSP585"),
      find_assigned_map(m, "2041-2070SSP585"),
      find_assigned_map(m, "2071-2100SSP585")
    )
    
    exists_vec <- !is.na(unlist(panel_files)) & file.exists(unlist(panel_files))
    
    if (!any(exists_vec)) {
      cat0("  [SKIP MODEL] no maps found for: ", m)
      next
    }
    
    if (!all(exists_vec)) {
      cat0("  [WARNING] missing one or more panels for model: ", m)
      print(data.table(panel = panel_titles, file = unlist(panel_files), exists = exists_vec))
    }
    
    files_keep <- panel_files[exists_vec]
    titles_keep <- panel_titles[exists_vec]
    
    index_list[[length(index_list) + 1L]] <- data.table(
      figure = paste0("Figure5", figure_letters[i]),
      method = m,
      method_label = map_method_labels[m],
      panel = panel_titles,
      file = unlist(panel_files),
      exists = exists_vec
    )
    
    plot_zone_panel_gg(
      files = files_keep,
      titles = titles_keep,
      outfile = file.path(
        fig_dir,
        paste0("Figure5", figure_letters[i], "_ecosystem_niche_maps_", m, "_2x4.png")
      ),
      categorical = TRUE,
      population_mode = FALSE,
      show_legend = FALSE,
      ncol = 4,
      width = 13.2,
      height = 6.8,
      draw_as_points = FALSE,
      point_alpha = 1,
      plot_title = paste0("Ecosystem niche projection: ", map_method_labels[m])
    )
  }
  
  if (length(index_list) == 0) stop("No ecosystem niche maps were found for any model.")
  
  idx <- rbindlist(index_list, fill = TRUE)
  fwrite(idx, file.path(tab_dir, "Figure5_ecosystem_niche_map_index_all_models_periods.csv"))
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
  out_list <- list()
  for (i in seq_len(length(stage_names) - 1L)) {
    from_col <- paste0("stage", i)
    to_col <- paste0("stage", i + 1L)
    tmp <- path_dt[
      ,
      .(area_km2 = sum(area_km2, na.rm = TRUE)),
      by = .(
        family,
        interval = paste0(stage_names[i], " -> ", stage_names[i + 1L]),
        from_zone = get(from_col),
        to_zone = get(to_col)
      )
    ]
    out_list[[length(out_list) + 1L]] <- tmp
  }
  rbindlist(out_list, fill = TRUE)
}

plot_transition_sankey <- function(path_dt, stage_names, outfile) {
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    stop("Package 'ggalluvial' is required for Figure 6b.")
  }
  if (nrow(path_dt) == 0) stop("No transition paths are available for Figure 6b.")
  
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
    ggalluvial::geom_text(
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
      title = paste0("Figure 6b. Ecotype transition trajectories: ", preferred_method),
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
  map_files <- sapply(
    scenario_order,
    function(s) find_assigned_map(preferred_method, s),
    USE.NAMES = TRUE
  )
  
  keep <- !is.na(map_files) & file.exists(map_files)
  map_files <- map_files[keep]
  
  if (length(map_files) == 0) stop("No assigned maps found for area calculation.")
  
  area_list <- list()
  
  for (s in names(map_files)) {
    cat0("  Calculating area: ", preferred_method, " | ", s)
    dt <- area_by_zone(map_files[[s]])
    dt[, scenario := s]
    dt[, scenario_label := scenario_label(s)]
    dt[, scenario_order := match(s, scenario_order)]
    area_list[[length(area_list) + 1L]] <- dt
  }
  
  area_dt <- rbindlist(area_list, fill = TRUE)
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
      title = paste0("Figure 6a. Projected ecotype area change: ", preferred_method),
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
  stage_names <- c("Reference", "2011-2040", "2041-2070", "2071-2100")
  
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

run_step("Table 3 | Ten species, number of populations, and Shannon H", {
  pop_tbl <- get_species_population_table()
  pop_agg <- pop_tbl[, .(abundance = .N), by = .(species, zoneID = source_zone)]
  
  tab3 <- pop_agg[
    ,
    {
      total <- sum(abundance, na.rm = TRUE)
      p <- abundance / total
      H <- -sum(p[p > 0] * log(p[p > 0]))
      .(n_populations = .N, total_abundance = total, shannon_H = H)
    },
    by = species
  ][order(-n_populations, -shannon_H)]
  
  tab3_top <- tab3[seq_len(min(10L, .N))]
  
  out_file <- file.path(tab_dir, "Table3_ten_species_populations_ShannonH.csv")
  fwrite(tab3_top, out_file)
  
  pop_file <- file.path(tab_dir, "Table3_species_population_abundance_long.csv")
  fwrite(pop_agg, pop_file)
  
  cat0("[SAVED] ", out_file)
  cat0("[SAVED] ", pop_file)
})

# 14. Figure 7 =================================================================

run_step("Figure 7 | Species x ecotype population abundance heatmap", {
  pop_file <- file.path(tab_dir, "Table3_species_population_abundance_long.csv")
  tab3_file <- file.path(tab_dir, "Table3_ten_species_populations_ShannonH.csv")
  
  if (!file.exists(pop_file) || !file.exists(tab3_file)) stop("Table 3 outputs not found.")
  
  pop <- fread(pop_file)
  tab3 <- fread(tab3_file)
  species_keep <- tab3$species
  
  heat <- pop[species %in% species_keep]
  if (nrow(heat) == 0) stop("No records for selected species.")
  
  heat[, zoneID_chr := as.character(zoneID)]
  heat[, log_abundance := log1p(abundance)]
  heat[, species := factor(species, levels = rev(species_keep))]
  heat[, zoneID_fac := factor(zoneID_chr, levels = as.character(sort(unique(zoneID))))]
  
  p <- ggplot(heat, aes(x = zoneID_fac, y = species, fill = log_abundance)) +
    geom_tile(colour = "white", linewidth = 0.18) +
    scale_fill_gradient(low = "#F5F5F2", high = "#3F4A4A", name = "log(1 + abundance)") +
    labs(
      title = "Species-by-ecotype population abundance",
      x = "Ecotype / vegetation zone",
      y = NULL
    ) +
    theme_ms() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.5),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank()
    )
  
  save_gg(
    p,
    file.path(fig_dir, "Figure7_species_ecotype_population_abundance_heatmap.png"),
    width = 10.5,
    height = 5.8
  )
})

# 15. Figure 8 =================================================================

run_step("Figure 8 | Tree population dual-niche projection maps", {
  pop_tbl <- get_species_population_table()
  fwrite(pop_tbl, file.path(tab_dir, "Figure8_species_population_source_zones.csv"))
  
  for (s in tree_plot_scenarios) {
    plot_species_population_dual_pages(
      pop_tbl = pop_tbl,
      method = tree_plot_method,
      scenario = s,
      outfile_prefix = file.path(fig_dir, "Figure8_tree_population_dual_niche_projection")
    )
  }
})

# 16. Figure 9 =================================================================

run_step("Figure 9 | Species-level dual-niche projection maps", {
  pop_tbl <- get_species_population_table()
  fwrite(pop_tbl, file.path(tab_dir, "Figure9_species_population_source_zones.csv"))
  
  for (s in tree_plot_scenarios) {
    plot_species_level_dual_pages(
      pop_tbl = pop_tbl,
      method = tree_plot_method,
      scenario = s,
      outfile_prefix = file.path(fig_dir, "Figure9_species_level_dual_niche_projection")
    )
  }
})

# 17. Table 4 ==================================================================

run_step("Table 4 | Binary optimized workflow vs multiclass robustness check", {
  files <- find_files("(robust|multiclass|multi_class).*\\.csv$")
  files <- files[!grepl("visualization", files, ignore.case = TRUE)]
  
  if (length(files) == 0) stop("No robustness or multiclass CSV files found.")
  
  tab_list <- list()
  
  for (f in files) {
    dt <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0) next
    dt[, source_file := relative_path(f)]
    tab_list[[length(tab_list) + 1L]] <- dt
  }
  
  if (length(tab_list) == 0) {
    stop("Robustness candidate files were found, but none could be read.")
  }
  
  tab4 <- rbindlist(tab_list, fill = TRUE)
  out_file <- file.path(tab_dir, "Table4_binary_vs_multiclass_robustness_check.csv")
  fwrite(tab4, out_file)
  cat0("[SAVED] ", out_file)
})

# 18. Figure 10a ================================================================

run_step("Figure 10a | Binary workflows vs multiclass reference-map F1 comparison", {
  bin_files <- find_files("normal_map_zone_metrics.*\\.csv$")
  bin_files <- bin_files[!grepl("visualization", bin_files, ignore.case = TRUE)]
  
  if (length(bin_files) == 0) stop("No binary normal_map_zone_metrics.csv files found.")
  
  bin_list <- list()
  
  for (f in bin_files) {
    dt <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0) next
    
    method_col <- pick_col(dt, c("method", "model", "workflow", "model_version", "version"))
    zone_col <- pick_col(dt, c("zoneID", "zone_id", "zone", "original_zone", "reference_zone", "class"))
    f1_col <- pick_col(dt, c("f1", "F1", "f1_score"))
    
    if (is.na(zone_col) || is.na(f1_col)) next
    
    if (!is.na(method_col)) {
      method_raw <- as.character(dt[[method_col]])
    } else {
      method_raw <- rep(normalise_map_method(f), nrow(dt))
    }
    
    tmp <- data.table(
      method_raw = method_raw,
      zoneID = extract_zone_id(dt[[zone_col]]),
      binary_f1 = as.numeric(dt[[f1_col]]),
      source_file = relative_path(f)
    )
    
    tmp[, method_key := normalise_map_method(method_raw)]
    tmp[is.na(method_key), method_key := normalise_map_method(source_file)]
    tmp <- tmp[method_key %in% map_method_order]
    
    bin_list[[length(bin_list) + 1L]] <- tmp
  }
  
  if (length(bin_list) == 0) stop("No usable binary zone-level F1 rows found.")
  
  bin <- rbindlist(bin_list, fill = TRUE)
  bin <- bin[zoneID %in% model_zoneID & !is.na(binary_f1) & is.finite(binary_f1)]
  bin <- bin[, .(binary_f1 = mean(binary_f1, na.rm = TRUE)), by = .(method_key, zoneID)]
  
  multi_files <- find_files("(multiclass|multi_class).*zone.*metrics.*\\.csv$")
  multi_files <- multi_files[!grepl("visualization", multi_files, ignore.case = TRUE)]
  multi_file <- pick_file(multi_files, prefer = c("multiclass", "zone", "metrics"))
  
  if (is.na(multi_file)) stop("No multiclass zone-metric CSV found.")
  
  multi <- fread(multi_file)
  multi_zone <- pick_col(multi, c("zoneID", "zone_id", "zone", "original_zone", "reference_zone", "class"))
  multi_f1 <- pick_col(multi, c("f1", "F1", "f1_score"))
  
  if (is.na(multi_zone) || is.na(multi_f1)) {
    stop("Cannot identify zone/F1 columns in multiclass file.")
  }
  
  multi_dt <- data.table(
    zoneID = extract_zone_id(multi[[multi_zone]]),
    multiclass_f1 = as.numeric(multi[[multi_f1]])
  )
  
  multi_dt <- multi_dt[
    zoneID %in% model_zoneID &
      !is.na(multiclass_f1) &
      is.finite(multiclass_f1)
  ]
  multi_dt <- multi_dt[, .(multiclass_f1 = mean(multiclass_f1, na.rm = TRUE)), by = zoneID]
  
  cmp <- merge(bin, multi_dt, by = "zoneID", all = FALSE)
  if (nrow(cmp) == 0) stop("No matched binary/multiclass zone-level F1 rows found.")
  
  area_dt <- zone_area_from_reference()[, .(zoneID, area_km2)]
  cmp <- merge(cmp, area_dt, by = "zoneID", all.x = TRUE)
  
  cmp[, method_label := map_method_labels[method_key]]
  cmp[, method_label := factor(method_label, levels = map_method_labels[map_method_order])]
  cmp[, zoneID_chr := as.character(zoneID)]
  
  cols <- zone_color_vector(cmp$zoneID)
  axis_lim <- c(0, 1)
  
  p <- ggplot(cmp, aes(x = binary_f1, y = multiclass_f1, colour = zoneID_chr, size = area_km2)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      colour = "grey45",
      linewidth = 0.45
    ) +
    geom_point(alpha = 0.90) +
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
      subtitle = "Map-overlay assessment against the original reference map, not testing-set F1. Point size is original zone area.",
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
    file.path(fig_dir, "Figure10a_all_binary_workflows_vs_multiclass_reference_map_F1_area_scaled.png"),
    width = 9.5,
    height = 8
  )
  
  fwrite(
    cmp,
    file.path(tab_dir, "Figure10a_all_binary_workflows_vs_multiclass_reference_map_F1.csv")
  )
})

# 19. Figure 10b ================================================================

find_multiclass_test_f1 <- function() {
  files <- find_files("(multiclass|multi_class|rf_test).*\\.csv$")
  files <- files[!grepl("visualization", files, ignore.case = TRUE)]
  
  candidates <- list()
  
  for (f in files) {
    dt <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0) next
    
    z_col <- pick_col(dt, c("zoneID", "zone_id", "zone", "class", "reference_zone"))
    if (is.na(z_col)) next
    
    f1_candidates <- names(dt)[grepl("test.*f1|f1.*test|^f1$|f1_score", norm_name(names(dt)))]
    f1_candidates <- f1_candidates[sapply(dt[, ..f1_candidates], is.numeric)]
    if (length(f1_candidates) == 0) next
    
    # Prefer explicitly testing-set F1; fall back to F1 if the file name itself is a test-metrics file.
    f1_col <- pick_col(dt, c("test_f1", "f1_test", "test_f1_score", "f1", "f1_score"))
    if (is.na(f1_col)) f1_col <- f1_candidates[1]
    
    tmp <- data.table(
      zoneID = extract_zone_id(dt[[z_col]]),
      multiclass_test_f1 = as.numeric(dt[[f1_col]]),
      multiclass_source_file = relative_path(f),
      multiclass_f1_col = f1_col
    )
    
    tmp <- tmp[
      zoneID %in% model_zoneID &
        !is.na(multiclass_test_f1) &
        is.finite(multiclass_test_f1)
    ]
    
    if (nrow(tmp) > 0) {
      candidates[[length(candidates) + 1L]] <- tmp
    }
  }
  
  if (length(candidates) == 0) return(data.table())
  
  out <- rbindlist(candidates, fill = TRUE)
  out <- out[
    ,
    .(
      multiclass_test_f1 = mean(multiclass_test_f1, na.rm = TRUE),
      multiclass_source_file = paste(unique(multiclass_source_file), collapse = "; "),
      multiclass_f1_col = paste(unique(multiclass_f1_col), collapse = "; ")
    ),
    by = zoneID
  ]
  
  out[]
}

run_step("Figure 10b | Binary workflows vs multiclass testing-set F1 comparison", {
  long_file <- file.path(tab_dir, "Table1_binary_RF_zone_level_metrics_long.csv")
  
  if (file.exists(long_file)) {
    long <- fread(long_file)
  } else {
    tbl <- make_binary_rf_table1()
    long <- tbl$long
    fwrite(long, long_file)
  }
  
  bin <- long[metric == "Test F1"]
  if (nrow(bin) == 0) {
    stop("No binary Test F1 found in Table1_binary_RF_zone_level_metrics_long.csv.")
  }
  
  # Binary climate/soil testing-set F1 are shown separately because these are binary climate/soil RFs.
  bin <- bin[
    ,
    .(binary_test_f1 = mean(value, na.rm = TRUE)),
    by = .(niche_type, workflow, zoneID)
  ]
  
  bin[, method_key := normalise_map_method(as.character(workflow))]
  bin <- bin[method_key %in% map_method_order]
  
  multi <- find_multiclass_test_f1()
  if (nrow(multi) == 0) {
    stop("No multiclass testing-set F1 could be found. Check multiclass/rf_test zone metrics files.")
  }
  
  cmp <- merge(bin, multi, by = "zoneID", all = FALSE)
  if (nrow(cmp) == 0) stop("No matched binary/multiclass testing-set F1 rows found.")
  
  area_dt <- zone_area_from_reference()[, .(zoneID, area_km2)]
  cmp <- merge(cmp, area_dt, by = "zoneID", all.x = TRUE)
  
  cmp[, method_label := map_method_labels[method_key]]
  cmp[, method_label := factor(method_label, levels = map_method_labels[map_method_order])]
  cmp[, niche_type := factor(niche_type, levels = c("climate", "soil"))]
  cmp[, zoneID_chr := as.character(zoneID)]
  
  cols <- zone_color_vector(cmp$zoneID)
  axis_lim <- c(0, 1)
  
  p <- ggplot(cmp, aes(x = binary_test_f1, y = multiclass_test_f1, colour = zoneID_chr, size = area_km2)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      colour = "grey45",
      linewidth = 0.45
    ) +
    geom_point(alpha = 0.90) +
    geom_text(
      aes(label = zoneID),
      check_overlap = TRUE,
      size = 2.0,
      nudge_y = 0.012,
      alpha = 1,
      show.legend = FALSE
    ) +
    facet_grid(niche_type ~ method_label) +
    scale_colour_manual(values = cols, guide = "none") +
    scale_size_continuous(
      range = c(1.0, 5.5),
      name = expression("Original area (km"^2*")")
    ) +
    scale_x_continuous(
      limits = axis_lim,
      breaks = seq(0, 1, 0.2),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      limits = axis_lim,
      breaks = seq(0, 1, 0.2),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    coord_equal(xlim = axis_lim, ylim = axis_lim, expand = FALSE) +
    labs(
      title = "Binary workflows vs multiclass testing-set F1",
      subtitle = "Testing-set comparison where available. Binary climate and soil RFs are shown separately.",
      x = "Binary RF testing-set F1",
      y = "Multiclass testing-set F1"
    ) +
    theme_ms(base_size = 9.5) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  
  save_gg(
    p,
    file.path(fig_dir, "Figure10b_all_binary_workflows_vs_multiclass_testing_set_F1_area_scaled.png"),
    width = 12.5,
    height = 7.5
  )
  
  fwrite(
    cmp,
    file.path(tab_dir, "Figure10b_all_binary_workflows_vs_multiclass_testing_set_F1.csv")
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

cat0("\nCOMPLETE")
cat0("Figure folder: ", fig_dir)
cat0("Table folder: ", tab_dir)
cat0("Step log: ", log_file)
