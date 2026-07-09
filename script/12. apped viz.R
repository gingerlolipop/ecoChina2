# 12. appendix visualization.R
# ============================================================
# Rich appendix / supplementary visualizations from existing results.
#
# Outputs:
#   visualization_appendix/tables/
#   visualization_appendix/figures/
#   visualization_appendix/appendix_step_log.csv
#
# Design:
#   - Each appendix table/figure is an independent step.
#   - Missing files are skipped without stopping the script.
#   - Zone-related figures use color_palette_China.csv whenever possible.
#   - Raster panels are drawn through ggplot to avoid empty base-plot panels.

library(terra)
library(data.table)
library(ggplot2)

rm(list = ls())
gc()

# 0. Settings ==================================================================

base_dir <- "H:/Jing/ecoChina2"

app_dir <- file.path(base_dir, "visualization_appendix")
fig_dir <- file.path(app_dir, "figures")
tab_dir <- file.path(app_dir, "tables")

dir.create(app_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

reference_file <- file.path(base_dir, "raster/ecosys_ori.tif")
result_map_root <- file.path(base_dir, "result maps")

model_zoneID <- c(1:7, 9:50, 52:55)
novel_value <- 99

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

binary_workflow_order <- c(
  "plain_rf",
  "plain_mf_rf",
  "optimized_rf",
  "optimized_mf_rf"
)

binary_workflow_labels <- c(
  plain_rf = "Plain RF",
  plain_mf_rf = "Plain MF RF",
  optimized_rf = "Optimized RF",
  optimized_mf_rf = "Optimized MF RF"
)

preferred_method <- "optimized_rf"

future_order <- c(
  "2011-2040SSP245",
  "2041-2070SSP245",
  "2071-2100SSP245",
  "2011-2040SSP585",
  "2041-2070SSP585",
  "2071-2100SSP585"
)

scenario_order <- c("normal", future_order)

display_max_cells <- 350000
fig_dpi <- 320

terraOptions(memfrac = 0.10)

niche_cols <- c(
  climate = "#2F5D7C",
  soil = "#6F7D45"
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

theme_app <- function(base_size = 10.5) {
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

save_gg <- function(p, file, width = 8, height = 5.5) {
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
      
      elapsed <- round(
        as.numeric(difftime(Sys.time(), t0, units = "secs")),
        2
      )
      
      step_log[[length(step_log) + 1L]] <<- data.table(
        step = step_name,
        status = "done",
        message = NA_character_,
        elapsed_sec = elapsed
      )
      
      cat0("[DONE] ", step_name, " | ", elapsed, " sec")
    },
    error = function(e) {
      elapsed <- round(
        as.numeric(difftime(Sys.time(), t0, units = "secs")),
        2
      )
      
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

# 2. Zone palette ==============================================================

default_zone_colors <- function(vals) {
  vals <- sort(unique(as.integer(vals)))
  vals <- vals[!is.na(vals)]
  
  cols <- grDevices::hcl.colors(
    n = length(vals),
    palette = "Dark 3"
  )
  
  names(cols) <- vals
  
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
  
  pal_files <- unique(c(
    file.path(base_dir, "color_palette_China.csv"),
    find_files("(color_palette|palette).*\\.csv$")
  ))
  
  pal_file <- pick_file(
    pal_files,
    prefer = c("color_palette_China.csv", "palette")
  )
  
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
    extra <- default_zone_colors(lk$zoneID[miss])
    lk[miss, color := unname(extra[as.character(zoneID)])]
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

# 3. Raster helpers =============================================================

thin_raster_for_plot <- function(x, categorical = TRUE, max_cells = display_max_cells) {
  x <- x[[1]]
  
  if (ncell(x) <= max_cells) return(x)
  
  fact <- ceiling(sqrt(ncell(x) / max_cells))
  fact <- max(1L, fact)
  
  tryCatch(
    {
      if (categorical) {
        aggregate(x, fact = fact, fun = "modal", na.rm = TRUE)
      } else {
        aggregate(x, fact = fact, fun = "mean", na.rm = TRUE)
      }
    },
    error = function(e) {
      spatSample(
        x,
        size = max_cells,
        method = "regular",
        as.raster = TRUE,
        na.rm = TRUE
      )
    }
  )
}

infer_source_zone_from_filename <- function(file) {
  b <- basename(file)
  
  patterns <- c(
    "source[_-]?zone[_-]?(\\d+)",
    "zone[_-]?(\\d+)",
    "z[_-]?(\\d+)",
    "population[_-]?(\\d+)",
    "pop[_-]?(\\d+)"
  )
  
  for (p in patterns) {
    if (grepl(p, b, ignore.case = TRUE)) {
      z <- sub(paste0(".*", p, ".*"), "\\1", b, ignore.case = TRUE)
      z <- suppressWarnings(as.integer(z))
      if (!is.na(z)) return(z)
    }
  }
  
  NA_integer_
}

raster_to_plot_dt <- function(file, categorical = TRUE, population_mode = FALSE) {
  
  if (inherits(file, "SpatRaster")) {
    x <- file[[1]]
    source_file <- NA_character_
  } else {
    if (is.na(file) || !file.exists(file)) stop("Missing raster: ", file)
    x <- rast(file)[[1]]
    source_file <- file
  }
  
  x <- thin_raster_for_plot(x, categorical = categorical)
  
  dt <- as.data.table(as.data.frame(
    x,
    xy = TRUE,
    na.rm = TRUE
  ))
  
  if (nrow(dt) == 0) return(data.table())
  
  value_col <- setdiff(names(dt), c("x", "y"))[1]
  setnames(dt, value_col, "value")
  
  dt <- dt[
    !is.na(value) &
      is.finite(value)
  ]
  
  if (nrow(dt) == 0) return(data.table())
  
  if (categorical) {
    dt[, value := as.integer(round(value))]
    
    if (population_mode && !is.na(source_file)) {
      vals <- sort(unique(dt$value))
      src_zone <- infer_source_zone_from_filename(source_file)
      
      if (all(vals %in% c(0L, 1L)) && !is.na(src_zone)) {
        dt <- dt[value == 1L]
        dt[, value := src_zone]
      } else {
        dt <- dt[value != 0L]
      }
    }
  }
  
  dt[]
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
    height = 6
) {
  
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
  
  if (length(plot_list) == 0) {
    stop("No valid raster cells were available for plotting.")
  }
  
  plot_dt <- rbindlist(plot_list, fill = TRUE)
  plot_dt[, panel := factor(panel, levels = titles[titles %in% unique(panel)])]
  
  if (categorical) {
    
    vals <- sort(unique(plot_dt$value))
    vals <- vals[!is.na(vals)]
    
    cols <- zone_color_vector(vals)
    
    if (novel_value %in% vals) {
      cols[as.character(novel_value)] <- "#333333"
    }
    
    plot_dt[, value_chr := factor(as.character(value), levels = as.character(vals))]
    
    p <- ggplot(
      plot_dt,
      aes(x = x, y = y, fill = value_chr)
    ) +
      geom_raster() +
      facet_wrap(~ panel, ncol = ncol) +
      coord_equal(expand = FALSE) +
      scale_fill_manual(values = cols, drop = FALSE) +
      labs(fill = "Zone") +
      theme_void(base_size = 10) +
      theme(
        strip.text = element_text(face = "bold", size = 9.5),
        legend.position = ifelse(show_legend, "right", "none"),
        legend.text = element_text(size = 7),
        legend.key.size = grid::unit(0.35, "cm")
      ) +
      guides(
        fill = guide_legend(ncol = 2, byrow = TRUE)
      )
    
  } else {
    
    p <- ggplot(
      plot_dt,
      aes(x = x, y = y, fill = value)
    ) +
      geom_raster() +
      facet_wrap(~ panel, ncol = ncol) +
      coord_equal(expand = FALSE) +
      scale_fill_gradientn(
        colours = grDevices::hcl.colors(60, palette = "YlGnBu"),
        na.value = NA
      ) +
      labs(fill = "Suitability") +
      theme_void(base_size = 10) +
      theme(
        strip.text = element_text(face = "bold", size = 9.5),
        legend.position = ifelse(show_legend, "right", "none"),
        legend.text = element_text(size = 7),
        legend.key.size = grid::unit(0.35, "cm")
      )
  }
  
  ggsave(
    outfile,
    p,
    width = width,
    height = height,
    dpi = fig_dpi,
    bg = "white"
  )
  
  cat0("[SAVED] ", outfile)
}

# 4. Scenario/map helpers =======================================================

scenario_label <- function(s) {
  ifelse(
    s == "normal",
    "Reference",
    paste0(
      sub("SSP[0-9]+$", "", s),
      "\n",
      sub("^.*(SSP[0-9]+)$", "\\1", s)
    )
  )
}

scenario_period <- function(s) {
  ifelse(s == "normal", "Reference", sub("SSP[0-9]+$", "", s))
}

scenario_ssp <- function(s) {
  ifelse(s == "normal", "Reference", sub("^.*(SSP[0-9]+)$", "\\1", s))
}

find_assigned_map <- function(method, scenario) {
  
  roots <- unique(c(
    file.path(result_map_root, method),
    result_map_root
  ))
  
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
  
  files <- files[
    !grepl("(color|rgb|plot|legend)", basename(files), ignore.case = TRUE)
  ]
  
  pick_file(
    files,
    prefer = c(
      "threshold0.2",
      "threshold0.1",
      "maskNA8",
      "noNovelNormal",
      method,
      scenario
    )
  )
}

area_by_zone <- function(file) {
  x <- rast(file)[[1]]
  a <- cellSize(x, unit = "km")
  z <- zonal(a, x, fun = "sum", na.rm = TRUE)
  
  dt <- as.data.table(z)
  
  if (ncol(dt) < 2) {
    stop("Cannot calculate area by zone: ", file)
  }
  
  setnames(dt, names(dt)[1], "zoneID")
  setnames(dt, names(dt)[2], "area_km2")
  
  dt[, zoneID := as.integer(round(zoneID))]
  dt <- dt[!is.na(zoneID)]
  dt <- dt[zoneID %in% c(model_zoneID, novel_value)]
  
  dt[]
}

# 5. Binary RF helpers ==========================================================

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
  
  if (grepl("oob.*acc|acc.*oob", x0)) return("OOB accuracy")
  if (grepl("train.*acc|acc.*train", x0)) return("Train accuracy")
  if (grepl("balanced.*acc|acc.*balanced", x0)) return("Balanced accuracy")
  if (grepl("auc|roc", x0)) return("AUC")
  if (grepl("f1", x0)) return("F1")
  if (grepl("precision|positive_predictive", x0)) return("Precision")
  if (grepl("recall|sensitivity|true_positive_rate", x0)) return("Recall")
  if (grepl("specificity|true_negative_rate", x0)) return("Specificity")
  if (grepl("threshold|cutoff", x0)) return("Threshold")
  if (grepl("(^|_)accuracy($|_)|(^|_)acc($|_)", x0)) return("Accuracy")
  
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

read_one_binary_rf_summary <- function(niche_type, file) {
  
  if (is.na(file) || !file.exists(file)) {
    stop("Missing binary RF summary file for ", niche_type)
  }
  
  dt <- fread(file)
  
  z_col <- pick_col(dt, c("zoneID", "zone_id", "zone", "ecotype", "class", "veg_zone"))
  
  if (is.na(z_col)) {
    stop("Cannot find zone column in: ", file)
  }
  
  dt[, zoneID_tmp := extract_zone_id(get(z_col))]
  dt <- dt[zoneID_tmp %in% model_zoneID]
  
  if (nrow(dt) == 0) {
    stop("No modeled zones found in: ", file)
  }
  
  workflow_col <- pick_col(dt, c(
    "workflow",
    "method",
    "model",
    "model_version",
    "version",
    "model_type",
    "rf_type"
  ))
  
  if (!is.na(workflow_col)) {
    dt[, workflow_row := infer_workflow(get(workflow_col))]
    dt[, workflow_source := workflow_col]
  } else {
    dt[, row_in_zone := seq_len(.N), by = zoneID_tmp]
    dt[, workflow_row := binary_workflow_order[((row_in_zone - 1L) %% length(binary_workflow_order)) + 1L]]
    dt[, workflow_source := "inferred_from_row_order"]
  }
  
  path_workflow <- infer_workflow(file)
  
  numeric_cols <- names(dt)[sapply(dt, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, c("zoneID_tmp", "row_in_zone"))
  
  out <- list()
  
  for (mc in numeric_cols) {
    
    metric <- metric_from_col(mc)
    if (is.na(metric)) next
    
    col_workflow <- infer_workflow(mc)
    
    tmp <- data.table(
      niche_type = niche_type,
      zoneID = dt$zoneID_tmp,
      metric = metric,
      value = as.numeric(dt[[mc]]),
      metric_raw = mc,
      source_file = relative_path(file)
    )
    
    if (all(is.na(col_workflow))) {
      if (!all(is.na(dt$workflow_row))) {
        tmp[, workflow := dt$workflow_row]
        tmp[, workflow_source := dt$workflow_source]
      } else if (!all(is.na(path_workflow))) {
        tmp[, workflow := path_workflow[1]]
        tmp[, workflow_source := "inferred_from_path"]
      } else {
        tmp[, workflow := NA_character_]
        tmp[, workflow_source := "unknown"]
      }
    } else {
      tmp[, workflow := col_workflow[1]]
      tmp[, workflow_source := "inferred_from_metric_column"]
    }
    
    out[[length(out) + 1L]] <- tmp
  }
  
  if (length(out) == 0) {
    stop("No usable metric columns found in: ", file)
  }
  
  out <- rbindlist(out, fill = TRUE)
  
  out <- out[
    workflow %in% binary_workflow_order &
      !is.na(value) &
      is.finite(value)
  ]
  
  out[]
}

make_binary_rf_metrics <- function() {
  
  climate_file <- find_binary_summary_file("climate")
  soil_file <- find_binary_summary_file("soil")
  
  cat0("[BINARY RF SOURCE] climate: ", climate_file)
  cat0("[BINARY RF SOURCE] soil: ", soil_file)
  
  long <- rbindlist(
    list(
      read_one_binary_rf_summary("climate", climate_file),
      read_one_binary_rf_summary("soil", soil_file)
    ),
    fill = TRUE
  )
  
  long <- long[
    ,
    .(
      value = mean(value, na.rm = TRUE),
      metric_raw = paste(unique(metric_raw), collapse = "; "),
      source_file = paste(unique(source_file), collapse = "; ")
    ),
    by = .(niche_type, workflow, zoneID, metric)
  ]
  
  long[, workflow := factor(workflow, levels = binary_workflow_order)]
  
  long[]
}

# 6. Appendix Table S1 ==========================================================

run_step("Appendix Table S1 | Result file inventory", {
  
  files <- unique(c(
    find_files("\\.csv$"),
    find_files("\\.tif$")
  ))
  
  files <- files[!grepl("visualization_appendix", files, ignore.case = TRUE)]
  
  inv <- data.table(
    file = relative_path(files),
    extension = tools::file_ext(files),
    size_mb = round(file.info(files)$size / 1024^2, 3),
    modified_time = as.character(file.info(files)$mtime)
  )
  
  inv[, category := fifelse(
    grepl("dual suit", file, ignore.case = TRUE), "dual suitability",
    fifelse(
      grepl("result maps|assigned_zone", file, ignore.case = TRUE), "assigned maps",
      fifelse(
        grepl("accuracy_climate|climate", file, ignore.case = TRUE), "climate RF",
        fifelse(
          grepl("accuracy_soil|soil", file, ignore.case = TRUE), "soil RF",
          fifelse(
            grepl("future tree niche", file, ignore.case = TRUE), "tree niche",
            fifelse(
              grepl("confusion|normal_map", file, ignore.case = TRUE), "reference assessment",
              "other"
            )
          )
        )
      )
    )
  )]
  
  setorder(inv, category, file)
  
  out_file <- file.path(tab_dir, "Appendix_TableS1_result_file_inventory.csv")
  fwrite(inv, out_file)
  
  cat0("[SAVED] ", out_file)
})

# 7. Appendix Figure S1 =========================================================

run_step("Appendix Figure S1 | Binary RF zone-level metric heatmaps", {
  
  dt <- make_binary_rf_metrics()
  
  out_file <- file.path(tab_dir, "Appendix_binary_RF_zone_level_metrics_long.csv")
  fwrite(dt, out_file)
  cat0("[SAVED] ", out_file)
  
  metrics_keep <- c("OOB accuracy", "Train accuracy", "Accuracy", "Balanced accuracy", "AUC", "F1")
  dt <- dt[metric %in% metrics_keep]
  
  if (nrow(dt) == 0) {
    stop("No selected binary RF metrics available.")
  }
  
  dt[, workflow_label := binary_workflow_labels[as.character(workflow)]]
  dt[, workflow_label := factor(workflow_label, levels = binary_workflow_labels[binary_workflow_order])]
  dt[, zoneID_chr := factor(as.character(zoneID), levels = as.character(model_zoneID))]
  dt[, metric := factor(metric, levels = metrics_keep)]
  dt[, niche_type := factor(niche_type, levels = c("climate", "soil"))]
  
  p <- ggplot(
    dt,
    aes(x = zoneID_chr, y = workflow_label, fill = value)
  ) +
    geom_tile(colour = "white", linewidth = 0.15) +
    facet_grid(niche_type ~ metric, scales = "free_x") +
    scale_fill_gradientn(
      colours = grDevices::hcl.colors(60, palette = "YlGnBu"),
      limits = c(0, 1),
      oob = scales::squish,
      name = "Metric"
    ) +
    labs(
      title = "Zone-level binary RF performance",
      x = "Vegetation zone",
      y = NULL
    ) +
    theme_app(base_size = 9.5) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
      axis.text.y = element_text(size = 7),
      panel.grid = element_blank()
    )
  
  save_gg(
    p,
    file.path(fig_dir, "Appendix_FigureS1_binary_RF_zone_metric_heatmaps.png"),
    width = 15,
    height = 6.5
  )
})

# 8. Appendix Figure S2 =========================================================

run_step("Appendix Figure S2 | Climate vs soil zone-level metric scatter", {
  
  metric_file <- file.path(tab_dir, "Appendix_binary_RF_zone_level_metrics_long.csv")
  
  if (file.exists(metric_file)) {
    dt <- fread(metric_file)
  } else {
    dt <- make_binary_rf_metrics()
  }
  
  metrics_keep <- c("AUC", "F1", "Accuracy", "Balanced accuracy")
  dt <- dt[metric %in% metrics_keep]
  
  wide <- dcast(
    dt,
    workflow + zoneID + metric ~ niche_type,
    value.var = "value",
    fun.aggregate = mean
  )
  
  wide <- wide[
    !is.na(climate) &
      !is.na(soil) &
      is.finite(climate) &
      is.finite(soil)
  ]
  
  if (nrow(wide) == 0) {
    stop("No matched climate/soil metric pairs found.")
  }
  
  wide[, workflow_label := binary_workflow_labels[as.character(workflow)]]
  wide[, workflow_label := factor(workflow_label, levels = binary_workflow_labels[binary_workflow_order])]
  wide[, zoneID_chr := as.character(zoneID)]
  
  cols <- zone_color_vector(unique(wide$zoneID))
  
  p <- ggplot(
    wide,
    aes(x = climate, y = soil, colour = zoneID_chr)
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
    geom_point(size = 1.7, alpha = 0.85) +
    facet_grid(metric ~ workflow_label) +
    scale_colour_manual(values = cols, guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    coord_equal(expand = FALSE) +
    labs(
      title = "Climate vs soil binary RF performance by zone",
      subtitle = "Each point is one vegetation zone; colors follow the zone palette.",
      x = "Climate RF metric",
      y = "Soil RF metric"
    ) +
    theme_app(base_size = 9.5)
  
  save_gg(
    p,
    file.path(fig_dir, "Appendix_FigureS2_climate_vs_soil_metric_scatter.png"),
    width = 12,
    height = 8.5
  )
})

# 9. Appendix Figure S3 =========================================================

run_step("Appendix Figure S3 | Reference maps for all workflows", {
  
  if (!file.exists(reference_file)) {
    stop("Missing reference raster: ", reference_file)
  }
  
  files <- list(reference_file)
  titles <- c("Original map")
  
  for (m in map_method_order) {
    f <- find_assigned_map(m, "normal")
    
    if (!is.na(f)) {
      files[[length(files) + 1L]] <- f
      titles <- c(titles, paste0(map_method_labels[m], "\nreference prediction"))
    }
  }
  
  if (length(files) <= 1) {
    stop("No predicted reference maps found.")
  }
  
  plot_zone_panel_gg(
    files = files,
    titles = titles,
    outfile = file.path(fig_dir, "Appendix_FigureS3_reference_maps_all_workflows.png"),
    categorical = TRUE,
    population_mode = FALSE,
    show_legend = FALSE,
    ncol = 3,
    width = 12.5,
    height = 7.8
  )
})

# 10. Appendix Figure S4 ========================================================

run_step("Appendix Figure S4 | Reference-map confusion matrix heatmap", {
  
  files <- find_files("normal_map_confusion_long.*\\.csv$")
  f <- pick_file(files, prefer = c(preferred_method, "normal_map_confusion_long"))
  
  if (is.na(f)) {
    stop("No normal_map_confusion_long.csv found.")
  }
  
  dt <- fread(f)
  
  from_col <- pick_col(dt, c("original_zone", "reference_zone", "actual_zone", "true_zone", "from", "truth", "observed"))
  to_col <- pick_col(dt, c("assigned_zone", "predicted_zone", "pred_zone", "prediction", "to"))
  count_col <- pick_col(dt, c("n", "count", "freq", "frequency", "pixels", "pixel_count"))
  
  if (is.na(from_col) || is.na(to_col) || is.na(count_col)) {
    stop("Cannot identify confusion columns.")
  }
  
  cm <- data.table(
    original_zone = extract_zone_id(dt[[from_col]]),
    predicted_zone = extract_zone_id(dt[[to_col]]),
    count = as.numeric(dt[[count_col]])
  )
  
  cm <- cm[
    original_zone %in% model_zoneID &
      predicted_zone %in% c(model_zoneID, novel_value) &
      !is.na(count) &
      is.finite(count)
  ]
  
  cm <- cm[
    ,
    .(count = sum(count, na.rm = TRUE)),
    by = .(original_zone, predicted_zone)
  ]
  
  cm[, row_total := sum(count), by = original_zone]
  cm[, prop := count / row_total]
  
  cm[, original_zone := factor(as.character(original_zone), levels = as.character(model_zoneID))]
  cm[, predicted_zone := factor(as.character(predicted_zone), levels = as.character(c(model_zoneID, novel_value)))]
  
  p <- ggplot(
    cm,
    aes(x = predicted_zone, y = original_zone, fill = prop)
  ) +
    geom_tile(colour = "white", linewidth = 0.15) +
    scale_fill_gradientn(
      colours = grDevices::hcl.colors(60, palette = "YlGnBu"),
      name = "Row proportion"
    ) +
    labs(
      title = paste0("Reference-map confusion matrix: ", preferred_method),
      x = "Predicted zone",
      y = "Original zone"
    ) +
    theme_app(base_size = 9.5) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      axis.text.y = element_text(size = 6),
      panel.grid = element_blank()
    )
  
  save_gg(
    p,
    file.path(fig_dir, "Appendix_FigureS4_reference_confusion_matrix_heatmap.png"),
    width = 9.5,
    height = 8.5
  )
})

# 11. Appendix Figure S5 ========================================================

run_step("Appendix Figure S5 | Zone metrics for all binary workflows", {
  
  files <- find_files("normal_map_zone_metrics.*\\.csv$")
  
  if (length(files) == 0) {
    stop("No normal_map_zone_metrics.csv found.")
  }
  
  all_dt <- list()
  
  for (f in files) {
    dt <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0) next
    
    method_col <- pick_col(dt, c("method", "model", "workflow", "model_version", "version"))
    
    if (!is.na(method_col)) {
      dt[, method := as.character(get(method_col))]
    } else {
      dt[, method := infer_workflow(f)]
    }
    
    z_col <- pick_col(dt, c("zoneID", "zone_id", "zone", "original_zone", "reference_zone", "class"))
    if (is.na(z_col)) next
    
    dt[, zoneID := extract_zone_id(get(z_col))]
    dt <- dt[zoneID %in% model_zoneID]
    
    metric_cols <- names(dt)[norm_name(names(dt)) %in% c("precision", "recall", "f1", "f1_score")]
    metric_cols <- metric_cols[sapply(dt[, ..metric_cols], is.numeric)]
    
    if (length(metric_cols) == 0) next
    
    mm <- melt(
      dt,
      id.vars = c("method", "zoneID"),
      measure.vars = metric_cols,
      variable.name = "metric",
      value.name = "value"
    )
    
    mm[, metric := fifelse(norm_name(metric) == "f1_score", "F1", tools::toTitleCase(metric))]
    all_dt[[length(all_dt) + 1L]] <- mm
  }
  
  if (length(all_dt) == 0) {
    stop("No zone-level metric rows found.")
  }
  
  m <- rbindlist(all_dt, fill = TRUE)
  m <- m[!is.na(value) & is.finite(value)]
  
  m[, method_key := gsub("_mf_rf", "_mf", infer_workflow(method))]
  m[is.na(method_key), method_key := method]
  m[, method_label := map_method_labels[method_key]]
  m[is.na(method_label), method_label := method_key]
  
  m[, method_label := factor(method_label, levels = map_method_labels[map_method_order])]
  m[, zoneID_chr := as.character(zoneID)]
  m[, zoneID_fac := factor(zoneID_chr, levels = as.character(model_zoneID))]
  
  cols <- zone_color_vector(model_zoneID)
  
  p <- ggplot(
    m,
    aes(x = zoneID_fac, y = value, colour = zoneID_chr)
  ) +
    geom_point(size = 1.4, alpha = 0.85) +
    facet_grid(metric ~ method_label) +
    scale_colour_manual(values = cols, guide = "none") +
    labs(
      title = "Zone-level reference-map metrics across workflows",
      subtitle = "Point colors follow the vegetation-zone palette.",
      x = "Vegetation zone",
      y = "Metric value"
    ) +
    theme_app(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
      panel.grid.major.x = element_blank()
    )
  
  save_gg(
    p,
    file.path(fig_dir, "Appendix_FigureS5_zone_metrics_all_workflows.png"),
    width = 13,
    height = 7.5
  )
})

# 12. Appendix Figure S6 ========================================================

run_step("Appendix Figure S6 | Future ecotype maps by workflow", {
  
  for (m in map_method_order) {
    
    future_files <- sapply(
      future_order,
      function(s) find_assigned_map(m, s),
      USE.NAMES = TRUE
    )
    
    keep <- !is.na(future_files) & file.exists(future_files)
    future_files <- future_files[keep]
    
    if (length(future_files) == 0) {
      cat0("  [SKIP METHOD] no future maps: ", m)
      next
    }
    
    plot_zone_panel_gg(
      files = as.list(future_files),
      titles = scenario_label(names(future_files)),
      outfile = file.path(
        fig_dir,
        paste0("Appendix_FigureS6_future_ecotype_maps_", m, ".png")
      ),
      categorical = TRUE,
      population_mode = FALSE,
      show_legend = FALSE,
      ncol = 3,
      width = 11.5,
      height = 7.2
    )
  }
})

# 13. Appendix Figure S7 ========================================================

run_step("Appendix Figure S7 | Novel area trajectory across workflows", {
  
  area_list <- list()
  
  for (m in map_method_order) {
    for (s in scenario_order) {
      f <- find_assigned_map(m, s)
      
      if (is.na(f) || !file.exists(f)) next
      
      dt <- area_by_zone(f)
      dt[, method := m]
      dt[, method_label := map_method_labels[m]]
      dt[, scenario := s]
      dt[, scenario_label := scenario_label(s)]
      dt[, scenario_order := match(s, scenario_order)]
      
      area_list[[length(area_list) + 1L]] <- dt
    }
  }
  
  if (length(area_list) == 0) {
    stop("No assigned maps available for novel-area calculation.")
  }
  
  area_dt <- rbindlist(area_list, fill = TRUE)
  
  out_file <- file.path(tab_dir, "Appendix_area_by_method_scenario_zone.csv")
  fwrite(area_dt, out_file)
  cat0("[SAVED] ", out_file)
  
  novel_dt <- area_dt[
    zoneID == novel_value,
    .(area_km2 = sum(area_km2, na.rm = TRUE)),
    by = .(method, method_label, scenario, scenario_label, scenario_order)
  ]
  
  all_combo <- unique(area_dt[, .(method, method_label, scenario, scenario_label, scenario_order)])
  novel_dt <- merge(
    all_combo,
    novel_dt,
    by = c("method", "method_label", "scenario", "scenario_label", "scenario_order"),
    all.x = TRUE
  )
  
  novel_dt[is.na(area_km2), area_km2 := 0]
  novel_dt[, area_million_km2 := area_km2 / 1e6]
  novel_dt[, scenario_label := factor(scenario_label, levels = scenario_label(scenario_order))]
  novel_dt[, method_label := factor(method_label, levels = map_method_labels[map_method_order])]
  
  p <- ggplot(
    novel_dt,
    aes(x = scenario_label, y = area_million_km2, group = method_label, colour = method_label)
  ) +
    geom_line(linewidth = 0.65) +
    geom_point(size = 2) +
    labs(
      title = "Novel ecotype area trajectory across workflows",
      x = NULL,
      y = expression("Novel area (" * 10^6 * " km"^2 * ")"),
      colour = "Workflow"
    ) +
    theme_app() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  
  save_gg(
    p,
    file.path(fig_dir, "Appendix_FigureS7_novel_area_trajectory_all_workflows.png"),
    width = 9,
    height = 5.2
  )
})

# 14. Appendix Figure S8 ========================================================

run_step("Appendix Figure S8 | All-zone area trajectories for optimized workflow", {
  
  area_file <- file.path(tab_dir, "Appendix_area_by_method_scenario_zone.csv")
  
  if (file.exists(area_file)) {
    area_dt <- fread(area_file)
  } else {
    stop("Area table not found. Run Appendix Figure S7 step first.")
  }
  
  dt <- area_dt[method == preferred_method]
  
  if (nrow(dt) == 0) {
    stop("No area rows for preferred method: ", preferred_method)
  }
  
  dt[, area_million_km2 := area_km2 / 1e6]
  dt[, zoneID_chr := as.character(zoneID)]
  dt[zoneID == novel_value, zoneID_chr := "Novel"]
  dt[, scenario_label := factor(scenario_label, levels = scenario_label(scenario_order))]
  
  cols <- zone_color_vector(unique(dt$zoneID))
  if (novel_value %in% dt$zoneID) cols[as.character(novel_value)] <- "#333333"
  names(cols)[names(cols) == as.character(novel_value)] <- "Novel"
  
  p <- ggplot(
    dt,
    aes(x = scenario_label, y = area_million_km2, group = zoneID_chr, colour = zoneID_chr)
  ) +
    geom_line(linewidth = 0.4, alpha = 0.8) +
    geom_point(size = 1.1, alpha = 0.8) +
    scale_colour_manual(values = cols, guide = "none") +
    labs(
      title = paste0("All-zone area trajectories: ", preferred_method),
      subtitle = "Each line is one projected zone; colors follow the zone palette.",
      x = NULL,
      y = expression("Area (" * 10^6 * " km"^2 * ")")
    ) +
    theme_app() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  
  save_gg(
    p,
    file.path(fig_dir, "Appendix_FigureS8_all_zone_area_trajectories_optimized_rf.png"),
    width = 9.5,
    height = 5.5
  )
})

# 15. Appendix Figure S9 ========================================================

run_step("Appendix Figure S9 | Dual-suitability ranking summary maps", {
  
  ranking_root <- file.path(base_dir, "dual suit ranking")
  
  if (!dir.exists(ranking_root)) {
    stop("Ranking folder not found: ", ranking_root)
  }
  
  scenarios <- c("normal", "2071-2100SSP245", "2071-2100SSP585")
  
  for (s in scenarios) {
    
    rank_dir <- file.path(ranking_root, preferred_method, s)
    
    summary_file <- file.path(rank_dir, "ranked_summary.tif")
    zone_file <- file.path(rank_dir, "ranked_zone.tif")
    
    files <- list()
    titles <- character()
    cat_flags <- logical()
    
    if (file.exists(zone_file)) {
      rz <- rast(zone_file)[[1]]
      files[[length(files) + 1L]] <- rz
      titles <- c(titles, paste0(scenario_label(s), "\nRank-1 zone"))
      cat_flags <- c(cat_flags, TRUE)
    }
    
    if (file.exists(summary_file)) {
      rs <- rast(summary_file)
      
      if ("n_zone_above_threshold" %in% names(rs)) {
        files[[length(files) + 1L]] <- rs[["n_zone_above_threshold"]]
        titles <- c(titles, paste0(scenario_label(s), "\nNo. zones above threshold"))
        cat_flags <- c(cat_flags, FALSE)
      }
      
      if ("top1_minus_top2" %in% names(rs)) {
        files[[length(files) + 1L]] <- rs[["top1_minus_top2"]]
        titles <- c(titles, paste0(scenario_label(s), "\nTop1 - Top2 margin"))
        cat_flags <- c(cat_flags, FALSE)
      }
    }
    
    if (length(files) == 0) {
      cat0("  [SKIP SCENARIO] no ranking rasters: ", s)
      next
    }
    
    # Plot rank-1 zone separately if present.
    if (any(cat_flags)) {
      plot_zone_panel_gg(
        files = files[cat_flags],
        titles = titles[cat_flags],
        outfile = file.path(
          fig_dir,
          paste0("Appendix_FigureS9_rank1_zone_", safe_name(s), ".png")
        ),
        categorical = TRUE,
        population_mode = FALSE,
        show_legend = TRUE,
        ncol = 1,
        width = 7,
        height = 5.8
      )
    }
    
    if (any(!cat_flags)) {
      plot_zone_panel_gg(
        files = files[!cat_flags],
        titles = titles[!cat_flags],
        outfile = file.path(
          fig_dir,
          paste0("Appendix_FigureS9_ranking_summary_", safe_name(s), ".png")
        ),
        categorical = FALSE,
        population_mode = FALSE,
        show_legend = TRUE,
        ncol = 2,
        width = 10,
        height = 5.4
      )
    }
  }
})

# 16. Appendix Table S2 and Figure S10 =========================================

run_step("Appendix Table S2 and Figure S10 | Full species population diversity", {
  
  files <- find_files("species_zone_population_long.*\\.csv$")
  f <- pick_file(files, prefer = c("species_zone_population_long"))
  
  if (is.na(f)) {
    stop("No species_zone_population_long.csv found.")
  }
  
  dt <- fread(f)
  
  species_col <- pick_col(dt, c("species", "species_name", "taxon", "tree_species"))
  zone_col <- pick_col(dt, c("zoneID", "zone_id", "zone", "ecotype", "source_zone", "population_zone"))
  
  if (is.na(species_col) || is.na(zone_col)) {
    stop("Cannot identify species and zone columns.")
  }
  
  abundance_col <- pick_col(dt, c("abundance", "n", "count", "cell_count", "pixel_count", "area_km2", "population_size"))
  
  if (!is.na(abundance_col) && is.numeric(dt[[abundance_col]])) {
    pop <- data.table(
      species = as.character(dt[[species_col]]),
      zoneID = extract_zone_id(dt[[zone_col]]),
      abundance = as.numeric(dt[[abundance_col]])
    )
  } else {
    pop <- unique(data.table(
      species = as.character(dt[[species_col]]),
      zoneID = extract_zone_id(dt[[zone_col]])
    ))
    
    pop[, abundance := 1]
  }
  
  pop <- pop[
    !is.na(species) &
      zoneID %in% model_zoneID &
      !is.na(abundance) &
      is.finite(abundance)
  ]
  
  if (nrow(pop) == 0) {
    stop("No valid species-population records.")
  }
  
  pop_agg <- pop[
    ,
    .(abundance = sum(abundance, na.rm = TRUE)),
    by = .(species, zoneID)
  ]
  
  species_summary <- pop_agg[
    ,
    {
      total <- sum(abundance, na.rm = TRUE)
      p <- abundance / total
      H <- -sum(p[p > 0] * log(p[p > 0]))
      
      .(
        n_populations = .N,
        total_abundance = total,
        shannon_H = H
      )
    },
    by = species
  ][order(-n_populations, -shannon_H)]
  
  out_file <- file.path(tab_dir, "Appendix_TableS2_all_species_population_diversity.csv")
  fwrite(species_summary, out_file)
  cat0("[SAVED] ", out_file)
  
  top_species <- species_summary[seq_len(min(25L, .N)), species]
  pdt <- species_summary[species %in% top_species]
  pdt[, species := factor(species, levels = rev(top_species))]
  
  p <- ggplot(
    pdt,
    aes(x = species, y = n_populations)
  ) +
    geom_col(width = 0.7, fill = "#4F6F52") +
    geom_point(aes(y = shannon_H * max(n_populations) / max(shannon_H)), size = 1.8, colour = "#2F5D7C") +
    coord_flip() +
    labs(
      title = "Population diversity across tree species",
      subtitle = "Bars show number of source-zone populations; points show Shannon H rescaled for display.",
      x = NULL,
      y = "Number of populations"
    ) +
    theme_app()
  
  save_gg(
    p,
    file.path(fig_dir, "Appendix_FigureS10_species_population_diversity.png"),
    width = 8,
    height = 7
  )
})

# 17. Appendix Figure S11 =======================================================

run_step("Appendix Figure S11 | Full species x zone abundance heatmap", {
  
  pop_file <- file.path(tab_dir, "Appendix_TableS2_all_species_population_diversity.csv")
  
  files <- find_files("species_zone_population_long.*\\.csv$")
  f <- pick_file(files, prefer = c("species_zone_population_long"))
  
  if (is.na(f)) {
    stop("No species_zone_population_long.csv found.")
  }
  
  dt <- fread(f)
  
  species_col <- pick_col(dt, c("species", "species_name", "taxon", "tree_species"))
  zone_col <- pick_col(dt, c("zoneID", "zone_id", "zone", "ecotype", "source_zone", "population_zone"))
  abundance_col <- pick_col(dt, c("abundance", "n", "count", "cell_count", "pixel_count", "area_km2", "population_size"))
  
  if (is.na(species_col) || is.na(zone_col)) {
    stop("Cannot identify species and zone columns.")
  }
  
  if (!is.na(abundance_col) && is.numeric(dt[[abundance_col]])) {
    pop <- data.table(
      species = as.character(dt[[species_col]]),
      zoneID = extract_zone_id(dt[[zone_col]]),
      abundance = as.numeric(dt[[abundance_col]])
    )
  } else {
    pop <- unique(data.table(
      species = as.character(dt[[species_col]]),
      zoneID = extract_zone_id(dt[[zone_col]])
    ))
    
    pop[, abundance := 1]
  }
  
  pop <- pop[
    !is.na(species) &
      zoneID %in% model_zoneID &
      !is.na(abundance) &
      is.finite(abundance)
  ]
  
  pop <- pop[
    ,
    .(abundance = sum(abundance, na.rm = TRUE)),
    by = .(species, zoneID)
  ]
  
  species_order <- pop[
    ,
    .(n_pop = .N, total = sum(abundance)),
    by = species
  ][order(-n_pop, -total)]$species
  
  pop[, species := factor(species, levels = rev(species_order))]
  pop[, zoneID_chr := factor(as.character(zoneID), levels = as.character(model_zoneID))]
  pop[, log_abundance := log1p(abundance)]
  
  p <- ggplot(
    pop,
    aes(x = zoneID_chr, y = species, fill = log_abundance)
  ) +
    geom_tile(colour = "white", linewidth = 0.12) +
    scale_fill_gradient(
      low = "#F5F5F2",
      high = "#3F4A4A",
      name = "log(1 + abundance)"
    ) +
    labs(
      title = "Full species-by-ecotype population abundance matrix",
      x = "Vegetation zone",
      y = NULL
    ) +
    theme_app(base_size = 8.5) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
      axis.text.y = element_text(size = 5.8),
      panel.grid = element_blank()
    )
  
  save_gg(
    p,
    file.path(fig_dir, "Appendix_FigureS11_full_species_zone_abundance_heatmap.png"),
    width = 11.5,
    height = 10
  )
})

# 18. Appendix Figure S12 =======================================================

run_step("Appendix Figure S12 | Additional population niche maps", {
  
  roots <- c(
    file.path(base_dir, "future tree niche"),
    file.path(base_dir, "future tree niche dual suitability")
  )
  
  roots <- roots[dir.exists(roots)]
  
  if (length(roots) == 0) {
    stop("No future tree niche folders found.")
  }
  
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
    grepl("(pop|population)", files, ignore.case = TRUE) &
      !grepl("species", basename(files), ignore.case = TRUE) &
      !grepl("(rank|summary|lookup|index)", basename(files), ignore.case = TRUE)
  ]
  
  files <- files[file.exists(files)]
  
  if (length(files) == 0) {
    stop("No population niche raster files found.")
  }
  
  score <- rep(0L, length(files))
  score <- score + as.integer(grepl(preferred_method, files, ignore.case = TRUE, fixed = TRUE))
  score <- score + as.integer(grepl("SSP585", files, ignore.case = TRUE, fixed = TRUE))
  score <- score + as.integer(grepl("2071-2100", files, ignore.case = TRUE, fixed = TRUE))
  
  files <- files[order(-score, files)]
  files <- files[seq_len(min(12L, length(files)))]
  
  titles <- tools::file_path_sans_ext(basename(files))
  titles <- gsub("_", " ", titles)
  titles <- substr(titles, 1, 55)
  
  plot_zone_panel_gg(
    files = as.list(files),
    titles = titles,
    outfile = file.path(fig_dir, "Appendix_FigureS12_additional_population_niche_maps.png"),
    categorical = TRUE,
    population_mode = TRUE,
    show_legend = TRUE,
    ncol = 4,
    width = 13.5,
    height = 9
  )
})

# 19. Appendix Figure S13 =======================================================

run_step("Appendix Figure S13 | Additional species-level niche maps", {
  
  roots <- c(
    file.path(base_dir, "future tree niche"),
    file.path(base_dir, "future tree niche dual suitability")
  )
  
  roots <- roots[dir.exists(roots)]
  
  if (length(roots) == 0) {
    stop("No future tree niche folders found.")
  }
  
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
    grepl("species", files, ignore.case = TRUE) &
      !grepl("(population|rank|summary|lookup|index)", basename(files), ignore.case = TRUE)
  ]
  
  files <- files[file.exists(files)]
  
  if (length(files) == 0) {
    stop("No species-level niche raster files found.")
  }
  
  score <- rep(0L, length(files))
  score <- score + as.integer(grepl(preferred_method, files, ignore.case = TRUE, fixed = TRUE))
  score <- score + as.integer(grepl("SSP585", files, ignore.case = TRUE, fixed = TRUE))
  score <- score + as.integer(grepl("2071-2100", files, ignore.case = TRUE, fixed = TRUE))
  
  files <- files[order(-score, files)]
  files <- files[seq_len(min(12L, length(files)))]
  
  titles <- tools::file_path_sans_ext(basename(files))
  titles <- gsub("_", " ", titles)
  titles <- substr(titles, 1, 55)
  
  plot_zone_panel_gg(
    files = as.list(files),
    titles = titles,
    outfile = file.path(fig_dir, "Appendix_FigureS13_additional_species_niche_maps.png"),
    categorical = FALSE,
    population_mode = FALSE,
    show_legend = TRUE,
    ncol = 4,
    width = 13.5,
    height = 9
  )
})

# 20. Save step log =============================================================

cat0("\n============================================================")
cat0("SAVE APPENDIX STEP LOG")
cat0("============================================================")

step_log_dt <- rbindlist(step_log, fill = TRUE)

log_file <- file.path(app_dir, "appendix_step_log.csv")
fwrite(step_log_dt, log_file)

cat0("[SAVED] ", log_file)

cat0("\nCOMPLETE")
cat0("Appendix figure folder: ", fig_dir)
cat0("Appendix table folder: ", tab_dir)
cat0("Appendix step log: ", log_file)