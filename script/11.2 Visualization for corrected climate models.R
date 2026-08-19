# 11.2 Main-text visualization: Plain RF and Plain MF RF
# ==============================================================================
# Final manuscript visualization uses ONLY the two correctly specified binary
# workflows:
#   rf_var : Plain RF
#   mf_var : Plain MF RF
#
# Both climate and soil models use zone-specific variables selected by mcRFop.
# Old incorrectly coded RF results are not used anywhere in this script; they
# remain available only for supplementary responses to reviewers.
#
# Main outputs:
#   Figure 1    Climate/soil independent-test metrics (parallel dot-range plot)
#   Figure 2    Normal-map reconstruction metrics
#   Figure 3    Zone-level climate/soil F1, precision and sensitivity
#   Figure 4    Major normal-map confusion Sankey
#               + zone/category chord PDFs for all normal-map transitions
#   Figure 5a-b Future assigned ecosystem maps, future only (SSP rows)
#   Figure 6    Projected novel ecosystem area
#   Figure 7    Normal-to-future transition shares
#   Figure 8    Assigned-zone species suitable area
#   Figure 9    Continuous dual-suitability species area
#   Figure 10a  Population suitable area
#   Figure 10b  Zone-level macro metrics (balanced accuracy, sensitivity, specificity, precision, F1, TSS), including Multiclass RF
#   Figure 10c  Zone-colored zone-level F1 vs Multiclass RF bubble comparison
#   Figure 11   Reference Top-k agreement and future rank retention
#
# Scientific rules:
#   soil gate applied upstream = 0.2
#   final suitability/novel threshold = 0.4
#   novel only when max dual suitability < 0.4
# ==============================================================================

library(terra)
library(data.table)
library(ggplot2)

rm(list = ls())
gc()


# 0. Paths and settings ==========================================================

base_dir <- "H:/Jing/ecoChina2"
assessment_dir <- file.path(base_dir, "assessment_var")
assigned_result_root <- file.path(base_dir, "future tree niche var")
assigned_table_dir <- file.path(assigned_result_root, "tables")
dual_result_root <- file.path(base_dir, "future tree niche dual suitability var")
dual_table_dir <- file.path(dual_result_root, "tables")
ranking_root <- file.path(base_dir, "dual suit ranking var")
ranking_table_dir <- file.path(ranking_root, "tables")
result_map_root <- file.path(base_dir, "result maps")
reference_file <- file.path(base_dir, "raster/ecosys_ori.tif")
palette_file <- file.path(base_dir, "color_palette_China.csv")
soil_model_root <- file.path(base_dir, "rf_soil")
soil_test_file <- file.path(base_dir, "results", "soil_test_data.csv")
model_zone_metric_file <- file.path(
  assessment_dir,
  "model_test_zone_metrics_var.csv"
)

output_root <- file.path(base_dir, "visualization var threshold0.4")
figure_dir <- file.path(output_root, "figures")
table_dir <- file.path(output_root, "tables")
chord_dir <- file.path(figure_dir, "chord diagrams")
assigned_species_page_dir <- file.path(figure_dir, "assigned species maps")
dual_species_page_dir <- file.path(figure_dir, "dual species maps")
population_detail_dir <- file.path(figure_dir, "population detail")

for (dir_name in c(
  figure_dir,
  table_dir,
  chord_dir,
  assigned_species_page_dir,
  dual_species_page_dir,
  population_detail_dir
)) {
  dir.create(dir_name, recursive = TRUE, showWarnings = FALSE)
}

new_threshold <- 0.4
tie_tol <- 1e-4
novel_value <- 99L
base_seed <- 49L
prob_threshold <- 0.5

method_order <- c(
  "rf_var",
  "mf_var"
)

method_labels <- c(
  rf_var = "Plain RF",
  mf_var = "Plain MF RF",
  multiclass_rf = "Multiclass RF"
)

method_colors <- c(
  rf_var = "#2E5E8C",
  mf_var = "#A36A26",
  multiclass_rf = "#4F4F4F"
)

main_method_colors <- setNames(
  unname(method_colors[method_order]),
  method_labels[method_order]
)

model_zoneID <- c(
  1:7,
  9:50,
  52:55
)

future_order <- c(
  "2011-2040SSP245",
  "2041-2070SSP245",
  "2071-2100SSP245",
  "2011-2040SSP585",
  "2041-2070SSP585",
  "2071-2100SSP585"
)

scenario_order <- c(
  "normal",
  future_order
)

period_levels <- c(
  "2011-2040",
  "2041-2070",
  "2071-2100"
)

ssp_levels <- c(
  "SSP245",
  "SSP585"
)

species_per_page <- 10L
species_page_columns <- 5L
species_page_rows <- 2L

multiclass_reference_map_file <- file.path(
  result_map_root,
  "multiclass_rf",
  "assigned_zone_normal_multiclass_rf.tif"
)


# 1. Helpers ====================================================================

safe_name <- function(x) {
  
  x <- gsub(
    "[^A-Za-z0-9_]+",
    "_",
    x
  )
  
  x <- gsub(
    "^_+|_+$",
    "",
    x
  )
  
  ifelse(
    nchar(x) == 0,
    "unnamed",
    x
  )
}


require_file <- function(file) {
  
  if (!file.exists(file)) {
    stop(
      "Missing required file: ",
      file
    )
  }
  
  file
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
    dpi = 300
  )
}


assigned_map_file <- function(
    method,
    scenario) {
  
  file.path(
    result_map_root,
    method,
    paste0(
      "assigned_zone_",
      scenario,
      "_threshold",
      new_threshold,
      "_tol",
      tie_tol,
      "_novel",
      novel_value,
      "_maskNA8_noNovelNormal.tif"
    )
  )
}


scenario_fields <- function(scenario) {
  
  data.table(
    period = sub(
      "SSP.*$",
      "",
      scenario
    ),
    ssp = sub(
      "^.*(SSP[0-9]+)$",
      "\\1",
      scenario
    )
  )
}


plot_zone_map <- function(
    raster_file,
    title,
    palette) {
  
  if (!file.exists(raster_file)) {
    
    plot.new()
    
    title(
      main = paste0(
        title,
        "\nmissing"
      )
    )
    
    return(invisible(FALSE))
  }
  
  raster <- rast(
    raster_file
  )
  
  indexed <- subst(
    raster,
    from = palette$zoneID,
    to = seq_len(
      nrow(palette)
    ),
    others = NA
  )
  
  plot(
    indexed,
    col = palette$COLOR,
    breaks = seq(
      0.5,
      nrow(palette) + 0.5,
      by = 1
    ),
    legend = FALSE,
    axes = FALSE,
    box = FALSE,
    main = title
  )
  
  invisible(TRUE)
}


plot_binary_species_page <- function(
    species_names,
    raster_files,
    page_title,
    output_file,
    reference_mask,
    fill_color = "#2E8B57") {
  
  png(
    output_file,
    width = 3200,
    height = 1500,
    res = 250
  )
  
  old_par <- par(
    no.readonly = TRUE
  )
  
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  
  par(
    mfrow = c(
      species_page_rows,
      species_page_columns
    ),
    mar = c(
      0.5,
      0.5,
      2,
      0.5
    ),
    oma = c(
      0,
      0,
      2.5,
      0
    )
  )
  
  for (index in seq_along(
    species_names
  )) {
    
    plot(
      reference_mask,
      col = "grey94",
      legend = FALSE,
      axes = FALSE,
      box = FALSE,
      main = species_names[index]
    )
    
    if (file.exists(
      raster_files[index]
    )) {
      
      species_raster <- rast(
        raster_files[index]
      )
      
      plot(
        species_raster,
        add = TRUE,
        col = grDevices::adjustcolor(
          fill_color,
          alpha.f = 0.30
        ),
        legend = FALSE
      )
    }
  }
  
  if (length(species_names) <
      species_per_page) {
    
    for (unused in seq_len(
      species_per_page -
      length(species_names)
    )) {
      plot.new()
    }
  }
  
  mtext(
    page_title,
    outer = TRUE,
    line = 0.5,
    cex = 1.2
  )
  
  invisible(TRUE)
}


global_mean <- function(raster) {
  
  as.numeric(
    global(
      raster,
      fun = "mean",
      na.rm = TRUE
    )[1, 1]
  )
}



div <- function(a, b) {
  
  n <- max(
    length(a),
    length(b)
  )
  
  a <- rep(
    a,
    length.out = n
  )
  
  b <- rep(
    b,
    length.out = n
  )
  
  out <- rep(
    NA_real_,
    n
  )
  
  valid <- (
    is.finite(b) &
      b > 0
  )
  
  out[valid] <- (
    a[valid] /
      b[valid]
  )
  
  out
}


zone_color_vector <- function(zone_values) {
  
  zone_values <- sort(unique(as.integer(zone_values)))
  zone_values <- zone_values[!is.na(zone_values)]
  
  color_map <- setNames(
    rep(NA_character_, length(zone_values)),
    as.character(zone_values)
  )
  
  matched <- match(
    zone_values,
    palette_map$zoneID
  )
  
  good <- !is.na(matched)
  
  color_map[as.character(zone_values[good])] <-
    palette_map$COLOR[matched[good]]
  
  color_map
}


calculate_zone_metrics_from_maps <- function(original_map, predicted_map) {
  
  original_modeled <- subst(
    original_map,
    from = model_zoneID,
    to = model_zoneID,
    others = NA
  )
  
  predicted_modeled <- subst(
    predicted_map,
    from = model_zoneID,
    to = model_zoneID,
    others = NA
  )
  
  names(original_modeled) <- "original_zone"
  names(predicted_modeled) <- "predicted_zone"
  
  confusion <- as.data.table(
    crosstab(
      c(
        original_modeled,
        predicted_modeled
      ),
      long = TRUE,
      useNA = FALSE
    )
  )
  
  setnames(
    confusion,
    names(confusion),
    c(
      "original_zone",
      "predicted_zone",
      "n"
    )
  )
  
  confusion[
    ,
    `:=`(
      original_zone = as.integer(original_zone),
      predicted_zone = as.integer(predicted_zone),
      n = as.numeric(n)
    )
  ]
  
  total_compared <- sum(
    confusion$n,
    na.rm = TRUE
  )
  
  rbindlist(
    lapply(
      model_zoneID,
      function(zone_value) {
        
        tp <- confusion[
          original_zone == zone_value & predicted_zone == zone_value,
          sum(n, na.rm = TRUE)
        ]
        
        fn <- confusion[
          original_zone == zone_value & predicted_zone != zone_value,
          sum(n, na.rm = TRUE)
        ]
        
        fp <- confusion[
          original_zone != zone_value & predicted_zone == zone_value,
          sum(n, na.rm = TRUE)
        ]
        
        tn <- total_compared - tp - fn - fp
        
        recall_value <- div(tp, tp + fn)
        specificity_value <- div(tn, tn + fp)
        precision_value <- div(tp, tp + fp)
        
        data.table(
          zone = zone_value,
          f1 = div(
            2 * precision_value * recall_value,
            precision_value + recall_value
          ),
          balanced_accuracy = div(
            recall_value + specificity_value,
            2
          )
        )
      }
    ),
    fill = TRUE
  )
}


plot_chord <- function(
    flow,
    item_order,
    item_color,
    item_label,
    out_file,
    main,
    label_cex = 0.65,
    label_track_height = 0.13) {
  
  flow <- as.data.table(copy(flow))
  
  flow[
    ,
    `:=`(
      from = as.character(from),
      to = as.character(to),
      n = as.numeric(n)
    )
  ]
  
  flow <- flow[n > 0]
  
  from_id <- item_order[item_order %in% unique(flow$from)]
  to_id <- rev(item_order[item_order %in% unique(flow$to)])
  
  from_sector <- paste0("O_", from_id)
  to_sector <- paste0("P_", to_id)
  sector_order <- c(from_sector, to_sector)
  
  flow[
    ,
    `:=`(
      from_sector = paste0("O_", from),
      to_sector = paste0("P_", to)
    )
  ]
  
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
  
  on.exit({
    circlize::circos.clear()
    dev.off()
  }, add = TRUE)
  
  par(mar = c(0.5, 0.5, 2.8, 0.5), xpd = NA)
  
  circlize::circos.par(
    start.degree = 90,
    clock.wise = FALSE,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.002, 0.002),
    canvas.xlim = c(-1.38, 1.38),
    canvas.ylim = c(-1.28, 1.28),
    points.overflow.warning = FALSE
  )
  
  circlize::chordDiagram(
    x = as.data.frame(flow[, .(from_sector, to_sector, n)]),
    order = sector_order,
    grid.col = sector_color,
    grid.border = NA,
    col = link_color,
    transparency = 0.72,
    directional = 1,
    direction.type = "diffHeight",
    diffHeight = circlize::mm_h(1),
    link.target.prop = FALSE,
    link.sort = "default",
    link.decreasing = TRUE,
    link.largest.ontop = TRUE,
    annotationTrack = "grid",
    annotationTrackHeight = circlize::mm_h(2),
    preAllocateTracks = list(track.height = label_track_height),
    big.gap = 14,
    small.gap = 0.15,
    reduce = -1
  )
  
  circlize::circos.trackPlotRegion(
    track.index = 1,
    bg.border = NA,
    panel.fun = function(x, y) {
      
      sector <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      ylim <- circlize::get.cell.meta.data("ylim")
      
      circlize::circos.text(
        x = mean(xlim),
        y = mean(ylim),
        labels = unname(sector_label[sector]),
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0.5, 0.5),
        cex = label_cex
      )
    }
  )
  
  mtext(main, side = 3, line = 0.5, font = 2, cex = 1.25)
  text(-1.27, 0, labels = "Original", srt = 90, font = 2, cex = 1.15)
  text(1.27, 0, labels = "Assigned", srt = 270, font = 2, cex = 1.15)
}


auc_rank <- function(y, p) {
  
  n1 <- as.numeric(sum(y == 1))
  n0 <- as.numeric(sum(y == 0))
  
  if (n1 == 0 || n0 == 0) {
    return(NA_real_)
  }
  
  (
    sum(rank(p, ties.method = "average")[y == 1]) -
      n1 * (n1 + 1) / 2
  ) / (n1 * n0)
}


load_rf_object <- function(file, object_name) {
  
  if (!file.exists(file)) {
    return(NULL)
  }
  
  model_env <- new.env()
  load(file, envir = model_env)
  
  if (!exists(object_name, envir = model_env)) {
    return(NULL)
  }
  
  get(object_name, envir = model_env)
}


balance_test <- function(test, zone, seed) {
  
  positive <- which(test$zoneID == zone)
  
  negative <- which(
    test$zoneID %in% model_zoneID &
      test$zoneID != zone
  )
  
  if (!length(positive) || !length(negative)) {
    return(integer())
  }
  
  set.seed(seed)
  
  if (length(negative) > length(positive)) {
    negative <- negative[
      sample.int(length(negative), length(positive))
    ]
  }
  
  c(positive, negative)
}


assess_selected_soil_models <- function() {
  
  soil_test <- as.data.frame(
    fread(
      require_file(soil_test_file)
    )
  )
  
  soil_test$zoneID <- as.numeric(
    as.character(soil_test$zoneID)
  )
  
  model_config <- data.table(
    method = method_order,
    file_prefix = c(
      "soil_plain_zone",
      "soil_mf_zone"
    ),
    object_name = c(
      "soil_plain",
      "soil_mf"
    )
  )
  
  test_index <- setNames(
    lapply(
      model_zoneID,
      function(zone_value) {
        balance_test(
          soil_test,
          zone_value,
          base_seed + 1000L + zone_value
        )
      }
    ),
    model_zoneID
  )
  
  results <- list()
  
  for (config_row in seq_len(nrow(model_config))) {
    
    method_name <- model_config$method[[config_row]]
    
    for (zone_value in model_zoneID) {
      
      model_file <- file.path(
        soil_model_root,
        paste0(
          model_config$file_prefix[[config_row]],
          zone_value,
          ".Rdata"
        )
      )
      
      model <- load_rf_object(
        model_file,
        model_config$object_name[[config_row]]
      )
      
      if (is.null(model)) {
        cat("[SKIP SOIL MODEL] ", method_name, " | zone ", zone_value, "\n", sep = "")
        next
      }
      
      # The current soil_plain and soil_mf models store the mcRFop-selected
      # zone-specific predictor set in $varlist.
      varlist <- model$varlist
      
      if (is.null(varlist) || !length(varlist)) {
        stop(
          "Selected-variable list missing from soil model: ",
          model_file
        )
      }
      
      if (!all(varlist %in% names(soil_test))) {
        stop(
          "Soil test data are missing selected predictors for ",
          method_name,
          " zone ",
          zone_value
        )
      }
      
      index <- test_index[[as.character(zone_value)]]
      
      if (!length(index)) {
        next
      }
      
      x <- soil_test[index, varlist, drop = FALSE]
      y <- as.integer(soil_test$zoneID[index] == zone_value)
      
      keep <- complete.cases(x)
      x <- x[keep, , drop = FALSE]
      y <- y[keep]
      
      if (!nrow(x) || length(unique(y)) < 2) {
        next
      }
      
      probability <- predict(
        model,
        x,
        type = "prob"
      )
      
      if (!("1" %in% colnames(probability))) {
        stop(
          "No presence-probability column in soil model: ",
          model_file
        )
      }
      
      probability <- as.numeric(probability[, "1"])
      finite <- is.finite(probability)
      probability <- probability[finite]
      y <- y[finite]
      
      prediction <- as.integer(
        probability >= prob_threshold
      )
      
      TP <- sum(y == 1 & prediction == 1)
      TN <- sum(y == 0 & prediction == 0)
      FP <- sum(y == 0 & prediction == 1)
      FN <- sum(y == 1 & prediction == 0)
      
      recall_value <- div(TP, TP + FN)
      specificity_value <- div(TN, TN + FP)
      precision_value <- div(TP, TP + FP)
      
      results[[length(results) + 1L]] <- data.table(
        niche = "soil",
        method = method_name,
        zone = zone_value,
        n_predictors = length(varlist),
        accuracy = div(TP + TN, length(y)),
        balanced_accuracy = div(
          recall_value + specificity_value,
          2
        ),
        recall = recall_value,
        specificity = specificity_value,
        precision = precision_value,
        f1 = div(
          2 * precision_value * recall_value,
          precision_value + recall_value
        ),
        tss = recall_value + specificity_value - 1,
        auc = auc_rank(y, probability)
      )
    }
  }
  
  if (!length(results)) {
    stop("No selected-variable soil models could be assessed.")
  }
  
  rbindlist(results, fill = TRUE)
}


# 2. Read current main-text source data ==========================================

model_zone_metrics <- fread(
  require_file(model_zone_metric_file)
)

model_zone_metrics <- model_zone_metrics[
  method %in% method_order &
    niche %in% c("climate", "soil") &
    zone %in% model_zoneID
]

expected_assessment_grid <- CJ(
  niche = c("climate", "soil"),
  method = method_order,
  zone = model_zoneID,
  unique = TRUE
)

observed_assessment_grid <- unique(
  model_zone_metrics[, .(niche, method, zone)]
)

missing_assessment_rows <- fsetdiff(
  expected_assessment_grid,
  observed_assessment_grid
)

duplicated_assessment_rows <- model_zone_metrics[
  ,
  .N,
  by = .(niche, method, zone)
][N != 1L]

nonfinite_auc_rows <- model_zone_metrics[!is.finite(auc)]

if (
  nrow(missing_assessment_rows) > 0L ||
  nrow(duplicated_assessment_rows) > 0L ||
  nrow(nonfinite_auc_rows) > 0L
) {
  if (nrow(missing_assessment_rows) > 0L) print(missing_assessment_rows)
  if (nrow(duplicated_assessment_rows) > 0L) print(duplicated_assessment_rows)
  if (nrow(nonfinite_auc_rows) > 0L) {
    print(nonfinite_auc_rows[, .(niche, method, zone, n_test, auc)])
  }
  
  stop(
    "Incomplete final assessment. Re-run script 5.1 after applying the ",
    "floating-point AUC fix; all 53 zones must have one finite AUC for each ",
    "niche and method."
  )
}

fwrite(
  model_zone_metrics,
  file.path(
    table_dir,
    "main_text_climate_soil_zone_metrics.csv"
  )
)

map_summary <- fread(
  require_file(
    file.path(
      assessment_dir,
      "normal_map_overall_metrics_var.csv"
    )
  )
)[
  method %in% method_order
]

ecosystem_area <- fread(
  require_file(
    file.path(
      assigned_table_dir,
      "future_ecosystem_area_var.csv"
    )
  )
)[
  method %in% method_order
]

transition <- fread(
  require_file(
    file.path(
      assigned_table_dir,
      "future_ecosystem_transition_var.csv"
    )
  )
)[
  method %in% method_order
]

assigned_species_area <- fread(
  require_file(
    file.path(
      assigned_table_dir,
      "future_species_niche_area_var.csv"
    )
  )
)[
  method %in% method_order
]

assigned_population_area <- fread(
  require_file(
    file.path(
      assigned_table_dir,
      "future_population_niche_area_var.csv"
    )
  )
)[
  method %in% method_order
]

dual_species_area <- fread(
  require_file(
    file.path(
      dual_table_dir,
      "dual_species_niche_area_var.csv"
    )
  )
)[
  method %in% method_order
]

dual_population_area <- fread(
  require_file(
    file.path(
      dual_table_dir,
      "dual_population_niche_area_var.csv"
    )
  )
)[
  method %in% method_order
]

ranking_index_file <- file.path(
  ranking_table_dir,
  "ranking_output_index_var.csv"
)

ranking_index <- if (file.exists(ranking_index_file)) {
  fread(ranking_index_file)[method %in% method_order]
} else {
  data.table()
}

palette <- fread(
  require_file(palette_file)
)

reference_map <- rast(
  require_file(reference_file)
)

reference_mask <- ifel(
  !is.na(reference_map),
  1,
  NA
)

palette[, zoneID := as.integer(zoneID)]

palette_map <- palette[
  zoneID != 8
][order(zoneID)]

if (!(novel_value %in% palette_map$zoneID)) {
  
  palette_map <- rbind(
    palette_map,
    data.table(
      zoneID = novel_value,
      zone = "Novel",
      COLOR = "#333333"
    ),
    fill = TRUE
  )
}

for (table_name in c(
  "ecosystem_area",
  "transition",
  "assigned_species_area",
  "assigned_population_area",
  "dual_species_area",
  "dual_population_area"
)) {
  
  table_object <- get(table_name)
  
  if ("period" %in% names(table_object)) {
    table_object[
      ,
      period := factor(period, levels = period_levels)
    ]
  }
  
  if ("ssp" %in% names(table_object)) {
    table_object[
      ,
      ssp := factor(ssp, levels = ssp_levels)
    ]
  }
  
  assign(table_name, table_object)
}


# 3. Climate and soil assessment ================================================

metric_columns <- c(
  "balanced_accuracy",
  "f1",
  "auc",
  "precision",
  "recall",
  "specificity",
  "tss"
)

metric_labels <- c(
  balanced_accuracy = "Balanced accuracy",
  f1 = "F1",
  auc = "AUC",
  precision = "Precision",
  recall = "Sensitivity",
  specificity = "Specificity",
  tss = "TSS"
)

performance_long <- melt(
  model_zone_metrics,
  id.vars = c(
    "niche",
    "method",
    "zone"
  ),
  measure.vars = metric_columns,
  variable.name = "metric",
  value.name = "value"
)

performance_long <- performance_long[
  is.finite(value)
]

performance_summary <- performance_long[
  ,
  .(
    n_zones = .N,
    mean = mean(value),
    sd = sd(value),
    se = sd(value) / sqrt(.N)
  ),
  by = .(
    niche,
    method,
    metric
  )
]

performance_summary[
  ,
  `:=`(
    ci95_low = mean - 1.96 * se,
    ci95_high = mean + 1.96 * se,
    metric_label = factor(
      metric_labels[metric],
      levels = unname(metric_labels[metric_columns])
    ),
    niche_label = factor(
      niche,
      levels = c("climate", "soil"),
      labels = c("Climate", "Soil")
    ),
    method_label = factor(
      method_labels[method],
      levels = method_labels[method_order]
    )
  )
]

fwrite(
  performance_summary,
  file.path(
    table_dir,
    "Figure_var_1_climate_soil_metric_summary.csv"
  )
)

performance_y_min <- max(
  0,
  floor(
    min(
      performance_summary$ci95_low,
      na.rm = TRUE
    ) * 20
  ) / 20 - 0.05
)

performance_long[
  ,
  `:=`(
    metric_label = factor(
      metric_labels[metric],
      levels = unname(metric_labels[metric_columns])
    ),
    niche_label = factor(
      niche,
      levels = c("climate", "soil"),
      labels = c("Climate", "Soil")
    ),
    method_label = factor(
      method_labels[method],
      levels = method_labels[method_order]
    ),
    zone_factor = factor(
      as.character(zone),
      levels = as.character(model_zoneID)
    )
  )
]

zone_point_colors <- zone_color_vector(model_zoneID)

# Metrics are parallel categories. Zone-level points keep the palette colors and are not connected.
figure_1 <- ggplot() +
  geom_point(
    data = performance_long,
    aes(
      x = metric_label,
      y = value,
      color = zone_factor,
      shape = method_label
    ),
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.52,
      seed = 1
    ),
    size = 1.55,
    alpha = 0.82,
    stroke = 0.15
  ) +
  geom_errorbar(
    data = performance_summary,
    aes(
      x = metric_label,
      y = mean,
      ymin = ci95_low,
      ymax = ci95_high,
      group = method_label
    ),
    position = position_dodge(width = 0.52),
    width = 0.12,
    linewidth = 0.6,
    color = "black"
  ) +
  geom_point(
    data = performance_summary,
    aes(
      x = metric_label,
      y = mean,
      shape = method_label
    ),
    position = position_dodge(width = 0.52),
    size = 3.2,
    fill = "white",
    color = "black",
    stroke = 0.55
  ) +
  facet_wrap(
    ~ niche_label,
    nrow = 1
  ) +
  scale_color_manual(
    values = zone_point_colors,
    guide = "none"
  ) +
  scale_shape_manual(
    values = c(21, 24)
  ) +
  coord_cartesian(
    ylim = c(performance_y_min, 1)
  ) +
  labs(
    x = NULL,
    y = "Independent-test metric",
    shape = NULL,
    title = "Climate and soil model performance",
    subtitle = "Colored points are vegetation zones; black symbols and error bars show mean +/- 1.96 SE."
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    ),
    panel.spacing = grid::unit(1.2, "lines")
  )

save_plot(
  figure_1,
  file.path(
    figure_dir,
    "Figure_var_1_climate_soil_model_performance.png"
  ),
  11.0,
  5.2
)


# 4. Normal-map assessment =======================================================

normal_map_zone_metrics <- fread(
  require_file(file.path(assessment_dir, "normal_map_zone_metrics_var.csv"))
)

map_metric_columns <- c(
  "f1",
  "precision",
  "recall",
  "tss"
)

map_metric_labels <- c(
  f1 = "F1",
  precision = "Precision",
  recall = "Sensitivity",
  tss = "TSS"
)

map_zone_long <- melt(
  normal_map_zone_metrics[
    method %in% method_order & zone %in% model_zoneID,
    c("method", "zone", map_metric_columns),
    with = FALSE
  ],
  id.vars = c("method", "zone"),
  measure.vars = map_metric_columns,
  variable.name = "metric",
  value.name = "value"
)

map_zone_long[
  ,
  `:=`(
    method_label = factor(
      method_labels[method],
      levels = method_labels[method_order]
    ),
    metric_label = factor(
      map_metric_labels[metric],
      levels = unname(map_metric_labels[map_metric_columns])
    ),
    zone_factor = factor(
      as.character(zone),
      levels = as.character(model_zoneID)
    )
  )
]

map_zone_summary <- map_zone_long[
  ,
  .(
    n_zones = .N,
    mean = mean(value),
    sd = sd(value),
    se = sd(value) / sqrt(.N)
  ),
  by = .(
    method,
    method_label,
    metric,
    metric_label
  )
]

map_zone_summary[
  ,
  `:=`(
    ci95_low = mean - 1.96 * se,
    ci95_high = mean + 1.96 * se
  )
]

figure_2 <- ggplot() +
  geom_point(
    data = map_zone_long,
    aes(
      x = metric_label,
      y = value,
      color = zone_factor,
      shape = method_label
    ),
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.52,
      seed = 2
    ),
    size = 1.65,
    alpha = 0.84,
    stroke = 0.15
  ) +
  geom_errorbar(
    data = map_zone_summary,
    aes(
      x = metric_label,
      y = mean,
      ymin = ci95_low,
      ymax = ci95_high,
      group = method_label
    ),
    position = position_dodge(width = 0.52),
    width = 0.12,
    linewidth = 0.6,
    color = "black"
  ) +
  geom_point(
    data = map_zone_summary,
    aes(
      x = metric_label,
      y = mean,
      shape = method_label
    ),
    position = position_dodge(width = 0.52),
    size = 3.25,
    fill = "white",
    color = "black",
    stroke = 0.55
  ) +
  scale_color_manual(
    values = zone_color_vector(model_zoneID),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c(21, 24)
  ) +
  coord_cartesian(
    ylim = c(0, 1)
  ) +
  labs(
    x = NULL,
    y = "Zone-level reconstruction metric",
    shape = NULL,
    title = "Normal-period map reconstruction",
    subtitle = "Colored points are vegetation zones; black symbols and error bars show mean +/- 1.96 SE."
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

save_plot(
  figure_2,
  file.path(
    figure_dir,
    "Figure_var_2_normal_map_assessment.png"
  ),
  8.2,
  5.0
)


# 5. Zone-level climate and soil metrics ========================================

zone_metric_long <- performance_long[
  metric %in% c(
    "f1",
    "precision",
    "recall",
    "tss"
  )
]

zone_metric_long[
  ,
  `:=`(
    metric_label = factor(
      metric_labels[metric],
      levels = c("F1", "Precision", "Sensitivity", "TSS")
    ),
    niche_label = factor(
      niche,
      levels = c("climate", "soil"),
      labels = c("Climate", "Soil")
    ),
    method_label = factor(
      method_labels[method],
      levels = method_labels[method_order]
    ),
    zone_factor = factor(
      as.character(zone),
      levels = rev(as.character(model_zoneID))
    )
  )
]

figure_3 <- ggplot(
  zone_metric_long,
  aes(
    x = value,
    y = zone_factor,
    color = method_label,
    shape = method_label
  )
) +
  geom_point(
    position = position_dodge(width = 0.55),
    size = 1.45,
    alpha = 0.85
  ) +
  facet_grid(
    niche_label ~ metric_label,
    scales = "free_x"
  ) +
  scale_color_manual(
    values = main_method_colors
  ) +
  scale_shape_manual(
    values = c(16, 17)
  ) +
  labs(
    x = "Independent-test metric",
    y = "Vegetation zone",
    color = NULL,
    shape = NULL,
    title = "Zone-level climate and soil performance"
  ) +
  theme_bw(base_size = 9.5) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 6.0),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = grid::unit(0.9, "lines")
  )

save_plot(
  figure_3,
  file.path(
    figure_dir,
    "Figure_var_3_zone_level_climate_soil_metrics.png"
  ),
  11.5,
  9.0
)


# 6. Normal assigned maps ========================================================

normal_files <- c(
  original = reference_file,
  rf_var = assigned_map_file("rf_var", "normal"),
  mf_var = assigned_map_file("mf_var", "normal")
)

normal_titles <- c(
  original = "Original vegetation map",
  rf_var = "Plain RF",
  mf_var = "Plain MF RF"
)

png(
  file.path(
    figure_dir,
    "Reference_maps_plain_rf_plain_mf.png"
  ),
  width = 3000,
  height = 1050,
  res = 250
)

old_par <- par(no.readonly = TRUE)

par(
  mfrow = c(1, 3),
  mar = c(1, 1, 2.2, 1)
)

for (map_name in names(normal_files)) {
  plot_zone_map(
    normal_files[map_name],
    normal_titles[map_name],
    palette_map
  )
}

par(old_par)
dev.off()


# 7. Future assigned maps ========================================================

future_titles <- c(
  `2011-2040SSP245` = "2011-2040",
  `2041-2070SSP245` = "2041-2070",
  `2071-2100SSP245` = "2071-2100",
  `2011-2040SSP585` = "2011-2040",
  `2041-2070SSP585` = "2041-2070",
  `2071-2100SSP585` = "2071-2100"
)

for (method in method_order) {
  
  png(
    file.path(
      figure_dir,
      paste0("Figure_var_5_future_maps_", method, ".png")
    ),
    width = 2800,
    height = 1900,
    res = 250
  )
  
  old_par <- par(no.readonly = TRUE)
  
  par(
    mfrow = c(2, 3),
    mar = c(1, 1, 2.4, 1),
    oma = c(0, 0, 2.5, 0)
  )
  
  for (scenario in future_order) {
    
    full_title <- paste(
      sub("SSP.*$", "", scenario),
      sub("^.*(SSP[0-9]+)$", "\\1", scenario)
    )
    
    plot_zone_map(
      assigned_map_file(method, scenario),
      full_title,
      palette_map
    )
  }
  
  mtext(
    paste0("Future top-ranked ecotype analogue maps | ", method_labels[method]),
    outer = TRUE,
    line = 0.7,
    cex = 1.2
  )
  
  par(old_par)
  dev.off()
}


# 7b. Normal confusion Sankey and chord diagrams ================================

if (requireNamespace("ggalluvial", quietly = TRUE)) {
  
  normal_map_confusion <- fread(
    require_file(
      file.path(assessment_dir, "normal_map_confusion_long_var.csv")
    )
  )
  
  rf_flow <- normal_map_confusion[
    method == "rf_var" & original_zone != predicted_zone,
    .(count = sum(n)),
    by = .(from = original_zone, to = predicted_zone)
  ]
  
  setorder(rf_flow, -count, from, to)
  rf_flow <- rf_flow[seq_len(min(20L, nrow(rf_flow)))]
  
  if (nrow(rf_flow) > 0) {
    
    rf_flow[
      ,
      `:=`(
        from_chr = as.character(from),
        to_chr = as.character(to)
      )
    ]
    
    figure_4 <- ggplot(
      rf_flow,
      aes(y = count, axis1 = from_chr, axis2 = to_chr)
    ) +
      ggalluvial::geom_alluvium(
        aes(fill = from_chr),
        width = 1 / 12,
        alpha = 0.8,
        show.legend = FALSE
      ) +
      ggalluvial::geom_stratum(
        width = 1 / 8,
        fill = "grey95",
        colour = "grey45"
      ) +
      ggalluvial::stat_stratum(
        geom = "text",
        aes(label = after_stat(stratum)),
        size = 2.5
      ) +
      scale_fill_manual(values = zone_color_vector(rf_flow$from)) +
      scale_x_discrete(
        limits = c("Original", "Assigned"),
        expand = c(0.08, 0.08)
      ) +
      labs(
        x = NULL,
        y = "Pixel count",
        title = "Major normal-map confusion flows | Plain RF"
      ) +
      theme_bw(base_size = 10.5) +
      theme(panel.grid = element_blank())
    
    save_plot(
      figure_4,
      file.path(figure_dir, "Figure_var_4_major_ecotype_confusion_flows.png"),
      8.8,
      6.2
    )
  }
}

if (requireNamespace("circlize", quietly = TRUE)) {
  
  normal_map_confusion <- fread(
    require_file(
      file.path(assessment_dir, "normal_map_confusion_long_var.csv")
    )
  )
  
  category_lookup <- unique(
    palette[
      zoneID %in% model_zoneID,
      .(zoneID, category2 = as.character(category2))
    ],
    by = "zoneID"
  )
  
  zone_to_category <- setNames(
    as.character(category_lookup$category2),
    as.character(category_lookup$zoneID)
  )
  
  zone_order <- as.character(model_zoneID)
  zone_colors <- zone_color_vector(model_zoneID)
  zone_labels <- setNames(as.character(model_zoneID), as.character(model_zoneID))
  
  category_palette <- palette[
    zoneID %in% model_zoneID,
    .(
      first_zone = min(zoneID),
      COLOR = COLOR[which.min(zoneID)]
    ),
    by = category2
  ][order(first_zone)]
  
  category_order <- as.character(category_palette$category2)
  category_colors <- setNames(category_palette$COLOR, category_order)
  category_labels <- setNames(
    gsub("_", " ", category_order, fixed = TRUE),
    category_order
  )
  
  for (method_name in method_order) {
    
    method_dt <- normal_map_confusion[
      method == method_name &
        original_zone %in% model_zoneID &
        predicted_zone %in% model_zoneID &
        n > 0
    ]
    
    zone_flow <- method_dt[
      ,
      .(n = sum(n)),
      by = .(
        from = as.character(original_zone),
        to = as.character(predicted_zone)
      )
    ]
    
    plot_chord(
      zone_flow,
      zone_order,
      zone_colors,
      zone_labels,
      file.path(chord_dir, paste0("normal_map_zone_chord_", method_name, ".pdf")),
      paste0(method_labels[method_name], ": Original zone to assigned zone")
    )
    
    category_flow <- copy(zone_flow)
    category_flow[, from := unname(zone_to_category[from])]
    category_flow[, to := unname(zone_to_category[to])]
    category_flow <- category_flow[
      !is.na(from) & !is.na(to),
      .(n = sum(n)),
      by = .(from, to)
    ]
    
    plot_chord(
      category_flow,
      category_order,
      category_colors,
      category_labels,
      file.path(chord_dir, paste0("normal_map_category_chord_", method_name, ".pdf")),
      paste0(method_labels[method_name], ": Original category to assigned category"),
      label_cex = 0.88,
      label_track_height = 0.18
    )
  }
}


# 8. Novel area and ecosystem transitions =======================================

novel_area <- ecosystem_area[zoneID == novel_value]

novel_area[
  ,
  method_label := factor(
    method_labels[method],
    levels = method_labels[method_order]
  )
]

figure_6 <- ggplot(
  novel_area,
  aes(x = period, y = area_km2, group = method_label, color = method_label)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.1) +
  scale_color_manual(values = main_method_colors) +
  facet_wrap(~ ssp, nrow = 1, scales = "free_y") +
  labs(
    x = "Future period",
    y = expression("Novel area (km"^2*")"),
    color = NULL,
    title = "Projected novel ecosystem area (threshold 0.4)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

save_plot(
  figure_6,
  file.path(figure_dir, "Figure_var_6_novel_area.png"),
  8.5,
  4.8
)

transition_summary <- transition[
  ,
  .(area_km2 = sum(area_km2, na.rm = TRUE)),
  by = .(method, scenario, period, ssp, transition_type)
]

transition_summary[
  ,
  total_area_km2 := sum(area_km2, na.rm = TRUE),
  by = .(method, scenario)
]

transition_summary[, area_share := area_km2 / total_area_km2]

transition_summary[
  ,
  `:=`(
    method_label = factor(
      method_labels[method],
      levels = method_labels[method_order]
    ),
    transition_type = factor(
      transition_type,
      levels = c("stable", "changed", "novel"),
      labels = c("Stable zone", "Changed zone", "Novel")
    )
  )
]

figure_7 <- ggplot(
  transition_summary,
  aes(x = period, y = area_share, fill = transition_type)
) +
  geom_col(width = 0.42, color = "white", linewidth = 0.15) +
  facet_grid(method_label ~ ssp) +
  scale_fill_manual(
    values = c(
      "Stable zone" = "#5B8E7D",
      "Changed zone" = "#D4A373",
      "Novel" = "#333333"
    )
  ) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  labs(
    x = "Future period",
    y = "Share of mapped area",
    fill = NULL,
    title = "Normal-to-future top-ranked analogue transitions"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.spacing = grid::unit(1.1, "lines")
  )

save_plot(
  figure_7,
  file.path(figure_dir, "Figure_var_7_transition_share.png"),
  9.8,
  6.0
)


# 9. Assigned-zone species and population areas =================================

assigned_species_area[
  ,
  method_label := factor(
    method_labels[method],
    levels = method_labels[method_order]
  )
]

figure_8 <- ggplot(
  assigned_species_area,
  aes(
    x = period,
    y = future_area_km2,
    group = method_label,
    shape = method_label,
    linetype = method_label
  )
) +
  geom_line(
    linewidth = 0.78,
    color = "black"
  ) +
  geom_point(
    size = 1.9,
    fill = "white",
    color = "black",
    stroke = 0.5
  ) +
  scale_shape_manual(values = c(21, 24)) +
  scale_linetype_manual(values = c("solid", "22")) +
  facet_grid(Species ~ ssp, scales = "free_y") +
  labs(
    x = "Future period",
    y = expression("Assigned-zone species area (km"^2*")"),
    shape = NULL,
    linetype = NULL,
    title = "Species niches from assigned ecosystem zones"
  ) +
  theme_bw(base_size = 9.2) +
  theme(
    legend.position = "top",
    strip.text.y = element_text(angle = 0, face = "italic"),
    strip.text.x = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "grey94", linewidth = 0.25),
    panel.spacing = grid::unit(0.85, "lines")
  )

save_plot(
  figure_8,
  file.path(figure_dir, "Figure_var_8_assigned_species_area.png"),
  10.4,
  max(8, 1.15 * length(unique(assigned_species_area$Species)))
)


# 10. Continuous dual-suitability niches ========================================

dual_species_long <- melt(
  dual_species_area,
  id.vars = c("Species", "method", "scenario"),
  measure.vars = c("suitable_area_km2", "suitability_weighted_area_km2"),
  variable.name = "area_metric",
  value.name = "area_value"
)

dual_species_long[
  ,
  period := factor(sub("SSP.*$", "", scenario), levels = period_levels)
]

dual_species_long[
  ,
  ssp := factor(sub("^.*(SSP[0-9]+)$", "\\1", scenario), levels = ssp_levels)
]

dual_species_long[
  ,
  `:=`(
    method_label = factor(
      method_labels[method],
      levels = method_labels[method_order]
    ),
    area_metric_label = factor(
      area_metric,
      levels = c("suitable_area_km2", "suitability_weighted_area_km2"),
      labels = c("Area with dual suitability >= 0.4", "Suitability-weighted area")
    ),
    Species = factor(Species, levels = rev(sort(unique(Species))))
  )
]

figure_9 <- ggplot(
  dual_species_long,
  aes(
    x = area_value,
    y = Species,
    shape = method_label
  )
) +
  geom_point(
    position = position_dodge(width = 0.55),
    size = 2.2,
    fill = "white",
    color = "black",
    stroke = 0.55
  ) +
  facet_grid(
    area_metric_label ~ ssp + period,
    scales = "free_x"
  ) +
  scale_shape_manual(values = c(21, 24)) +
  labs(
    x = expression("Species niche area (km"^2*")"),
    y = "Species",
    shape = NULL,
    title = "Species niches from continuous dual suitability",
    subtitle = "Each panel is one future scenario. Circle = Plain RF; triangle = Plain MF RF."
  ) +
  theme_bw(base_size = 9.0) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    strip.text.y = element_text(angle = 0),
    strip.text.x = element_text(face = "bold"),
    panel.spacing = grid::unit(0.8, "lines")
  )

save_plot(
  figure_9,
  file.path(figure_dir, "Figure_var_9_dual_species_area.png"),
  13.0,
  7.8
)

dual_population_area[
  ,
  period := factor(
    sub("SSP.*$", "", scenario),
    levels = period_levels
  )
]

dual_population_area[
  ,
  ssp := factor(
    sub("^.*(SSP[0-9]+)$", "\\1", scenario),
    levels = ssp_levels
  )
]

dual_population_area[
  ,
  method_label := factor(
    method_labels[method],
    levels = method_labels[method_order]
  )
]

# Detailed population visualization: one species per page/file.
# Population labels retain source-zone information because each population
# represents a distinct evolutionary environment.
population_detail <- copy(dual_population_area)

population_detail[
  ,
  population_label := ifelse(
    !is.na(zone_name) & nzchar(zone_name),
    paste0("P", PopulationID, " | Z", source_zone, " ", zone_name),
    paste0("P", PopulationID, " | Z", source_zone)
  )
]

fwrite(
  population_detail,
  file.path(
    table_dir,
    "Figure_var_10a_population_detail_long.csv"
  )
)

population_species <- sort(unique(population_detail$Species))

population_index <- list()

for (species_name in population_species) {
  
  species_dt <- population_detail[
    Species == species_name
  ]
  
  population_order <- species_dt[
    ,
    .(
      source_zone = first(source_zone)
    ),
    by = .(
      PopulationID,
      population_label
    )
  ][
    order(source_zone, PopulationID)
  ]$population_label
  
  species_dt[
    ,
    population_label := factor(
      population_label,
      levels = rev(population_order)
    )
  ]
  
  # 10a-1: absolute suitable area for every population.
  p_area <- ggplot(
    species_dt,
    aes(
      x = suitable_area_km2,
      y = population_label
    )
  ) +
    geom_point(
      shape = 16,
      size = 1.9,
      color = "black"
    ) +
    facet_grid(
      method_label + ssp ~ period,
      scales = "free_x"
    ) +
    labs(
      x = expression("Suitable area (km"^2*")"),
      y = "Population and source zone",
      title = paste0(
        species_name,
        ": population-level future niche area"
      ),
      subtitle = paste(
        "Each point is one source population.",
        "Rows compare Plain RF vs Plain MF RF; within each row, SSP245 and SSP585 are separated."
      )
    ) +
    theme_bw(base_size = 9.6) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 6.8),
      axis.text.x = element_text(size = 8.2),
      panel.spacing = grid::unit(0.75, "lines"),
      legend.position = "none",
      plot.subtitle = element_text(size = 9)
    )
  
  area_file <- file.path(
    population_detail_dir,
    paste0(
      "Figure_var_10a_population_area_",
      safe_name(species_name),
      ".png"
    )
  )
  
  save_plot(
    p_area,
    area_file,
    14.2,
    max(5.6, 2.8 + 0.30 * length(population_order))
  )
  
  # 10a-2: mean dual suitability for the same populations.
  p_suit <- ggplot(
    species_dt,
    aes(
      x = mean_dual_suitability,
      y = population_label
    )
  ) +
    geom_point(
      shape = 16,
      size = 1.9,
      color = "black"
    ) +
    facet_grid(
      method_label + ssp ~ period,
      scales = "fixed"
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2)
    ) +
    labs(
      x = "Mean dual suitability",
      y = "Population and source zone",
      title = paste0(
        species_name,
        ": population-level mean dual suitability"
      ),
      subtitle = paste(
        "Each point is one source population.",
        "This figure emphasizes differences among evolutionary source environments."
      )
    ) +
    theme_bw(base_size = 9.6) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 6.8),
      axis.text.x = element_text(size = 8.2),
      panel.spacing = grid::unit(0.75, "lines"),
      legend.position = "none",
      plot.subtitle = element_text(size = 9)
    )
  
  suit_file <- file.path(
    population_detail_dir,
    paste0(
      "Figure_var_10a_population_mean_suitability_",
      safe_name(species_name),
      ".png"
    )
  )
  
  save_plot(
    p_suit,
    suit_file,
    14.2,
    max(5.6, 2.8 + 0.30 * length(population_order))
  )
  
  population_index[[length(population_index) + 1L]] <- data.table(
    Species = species_name,
    n_populations = length(population_order),
    area_figure = area_file,
    mean_suitability_figure = suit_file
  )
}

population_index <- rbindlist(
  population_index,
  fill = TRUE
)

fwrite(
  population_index,
  file.path(
    table_dir,
    "Figure_var_10a_population_detail_index.csv"
  )
)

# Cross-species population comparison after the within-species detail views.
population_species_summary <- population_detail[
  ,
  .(
    n_populations = uniqueN(PopulationID),
    median_suitable_area_km2 = median(suitable_area_km2, na.rm = TRUE),
    q25_suitable_area_km2 = quantile(suitable_area_km2, 0.25, na.rm = TRUE, names = FALSE),
    q75_suitable_area_km2 = quantile(suitable_area_km2, 0.75, na.rm = TRUE, names = FALSE)
  ),
  by = .(
    Species,
    method_label,
    period,
    ssp
  )
]

figure_10a_summary <- ggplot(
  population_species_summary,
  aes(
    x = median_suitable_area_km2,
    y = reorder(Species, median_suitable_area_km2),
    shape = method_label
  )
) +
  geom_errorbarh(
    aes(
      xmin = q25_suitable_area_km2,
      xmax = q75_suitable_area_km2
    ),
    height = 0.18,
    position = position_dodge(width = 0.5),
    linewidth = 0.5,
    color = "grey45"
  ) +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 2.2,
    fill = "white",
    color = "black",
    stroke = 0.55
  ) +
  facet_grid(
    ssp ~ period,
    scales = "free_x"
  ) +
  scale_shape_manual(values = c(21, 24)) +
  labs(
    x = expression("Median population suitable area (km"^2*")"),
    y = "Species",
    shape = NULL,
    title = "Across-species comparison of population niche area",
    subtitle = "Points are species medians across populations; horizontal bars show the interquartile range."
  ) +
  theme_bw(base_size = 9.2) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.spacing = grid::unit(0.9, "lines")
  )

save_plot(
  figure_10a_summary,
  file.path(
    figure_dir,
    "Figure_var_10a_population_species_summary.png"
  ),
  12.0,
  7.5
)

normal_map_zone_metrics <- fread(
  require_file(file.path(assessment_dir, "normal_map_zone_metrics_var.csv"))
)

# Recalculate all three reference maps on one common valid mask so that the
# overall comparison is directly comparable across Plain RF, Plain MF RF and
# Multiclass RF.
reference_methods <- c(
  method_order,
  "multiclass_rf"
)

reference_map_files <- c(
  setNames(
    vapply(
      method_order,
      function(method_key) {
        assigned_map_file(method_key, "normal")
      },
      character(1)
    ),
    method_order
  ),
  multiclass_rf = multiclass_reference_map_file
)

for (map_file in reference_map_files) {
  require_file(map_file)
}

original_modeled <- subst(
  reference_map,
  from = model_zoneID,
  to = model_zoneID,
  others = NA
)

names(original_modeled) <- "original_zone"

predicted_maps <- lapply(
  reference_methods,
  function(method_key) {
    
    predicted <- rast(reference_map_files[[method_key]])[[1]]
    
    if (!compareGeom(
      original_modeled,
      predicted,
      stopOnError = FALSE
    )) {
      stop(
        "Geometry mismatch for reference-map comparison: ",
        method_key
      )
    }
    
    predicted <- subst(
      predicted,
      from = model_zoneID,
      to = model_zoneID,
      others = NA
    )
    
    names(predicted) <- method_key
    predicted
  }
)

names(predicted_maps) <- reference_methods

valid_original_pixels <- as.numeric(
  global(
    !is.na(original_modeled),
    fun = "sum",
    na.rm = TRUE
  )[1, 1]
)

individual_coverage_dt <- rbindlist(
  lapply(
    reference_methods,
    function(method_key) {
      
      individual_mask <-
        !is.na(original_modeled) &
        !is.na(predicted_maps[[method_key]])
      
      compared_pixels <- as.numeric(
        global(
          individual_mask,
          fun = "sum",
          na.rm = TRUE
        )[1, 1]
      )
      
      data.table(
        method = method_key,
        individual_map_coverage = div(
          compared_pixels,
          valid_original_pixels
        )
      )
    }
  ),
  fill = TRUE
)

common_valid <- !is.na(original_modeled)

for (method_key in reference_methods) {
  common_valid <-
    common_valid &
    !is.na(predicted_maps[[method_key]])
}

common_pixels <- as.numeric(
  global(
    common_valid,
    fun = "sum",
    na.rm = TRUE
  )[1, 1]
)

if (!is.finite(common_pixels) || common_pixels <= 0) {
  stop(
    "No common valid pixels across Plain RF, Plain MF RF and Multiclass RF."
  )
}

original_common <- ifel(
  common_valid,
  original_modeled,
  NA
)

names(original_common) <- "original_zone"

reference_overall_list <- list()
reference_zone_list <- list()

for (method_key in reference_methods) {
  
  predicted_common <- ifel(
    common_valid,
    predicted_maps[[method_key]],
    NA
  )
  
  names(predicted_common) <- "predicted_zone"
  
  confusion <- as.data.table(
    crosstab(
      c(
        original_common,
        predicted_common
      ),
      long = TRUE,
      useNA = FALSE
    )
  )
  
  setnames(
    confusion,
    names(confusion),
    c(
      "original_zone",
      "predicted_zone",
      "n"
    )
  )
  
  confusion[
    ,
    `:=`(
      original_zone = as.integer(original_zone),
      predicted_zone = as.integer(predicted_zone),
      n = as.numeric(n)
    )
  ]
  
  confusion <- confusion[
    original_zone %in% model_zoneID &
      predicted_zone %in% model_zoneID &
      n > 0
  ]
  
  total_compared <- sum(
    confusion$n,
    na.rm = TRUE
  )
  
  zone_metrics <- rbindlist(
    lapply(
      model_zoneID,
      function(zone_value) {
        
        tp <- confusion[
          original_zone == zone_value &
            predicted_zone == zone_value,
          sum(n, na.rm = TRUE)
        ]
        
        fn <- confusion[
          original_zone == zone_value &
            predicted_zone != zone_value,
          sum(n, na.rm = TRUE)
        ]
        
        fp <- confusion[
          original_zone != zone_value &
            predicted_zone == zone_value,
          sum(n, na.rm = TRUE)
        ]
        
        tn <- total_compared - tp - fn - fp
        
        recall_value <- div(
          tp,
          tp + fn
        )
        
        specificity_value <- div(
          tn,
          tn + fp
        )
        
        precision_value <- div(
          tp,
          tp + fp
        )
        
        data.table(
          method = method_key,
          zone = zone_value,
          recall = recall_value,
          specificity = specificity_value,
          precision = precision_value,
          balanced_accuracy = div(
            recall_value + specificity_value,
            2
          ),
          f1 = div(
            2 * precision_value * recall_value,
            precision_value + recall_value
          ),
          tss = recall_value + specificity_value - 1
        )
      }
    ),
    fill = TRUE
  )
  
  exact_accuracy <- div(
    confusion[
      original_zone == predicted_zone,
      sum(n, na.rm = TRUE)
    ],
    total_compared
  )
  
  broad_accuracy <- NA_real_
  
  category_lookup_10b <- unique(
    palette[
      zoneID %in% model_zoneID,
      .(
        zoneID,
        category2 = as.character(category2)
      )
    ],
    by = "zoneID"
  )
  
  if (nrow(category_lookup_10b) > 0) {
    
    category_map <- setNames(
      category_lookup_10b$category2,
      as.character(category_lookup_10b$zoneID)
    )
    
    confusion[
      ,
      `:=`(
        original_category = unname(
          category_map[as.character(original_zone)]
        ),
        predicted_category = unname(
          category_map[as.character(predicted_zone)]
        )
      )
    ]
    
    if (
      !anyNA(confusion$original_category) &&
      !anyNA(confusion$predicted_category)
    ) {
      broad_accuracy <- div(
        confusion[
          original_category == predicted_category,
          sum(n, na.rm = TRUE)
        ],
        total_compared
      )
    }
  }
  
  reference_zone_list[[method_key]] <- zone_metrics
  
  reference_overall_list[[method_key]] <- data.table(
    method = method_key,
    individual_map_coverage = individual_coverage_dt[
      method == method_key,
      individual_map_coverage
    ],
    common_comparison_coverage = div(
      common_pixels,
      valid_original_pixels
    ),
    exact_zone_accuracy = exact_accuracy,
    broad_category_accuracy = broad_accuracy,
    macro_balanced_accuracy = mean(
      zone_metrics$balanced_accuracy,
      na.rm = TRUE
    ),
    macro_recall = mean(
      zone_metrics$recall,
      na.rm = TRUE
    ),
    macro_specificity = mean(
      zone_metrics$specificity,
      na.rm = TRUE
    ),
    macro_precision = mean(
      zone_metrics$precision,
      na.rm = TRUE
    ),
    macro_f1 = mean(
      zone_metrics$f1,
      na.rm = TRUE
    ),
    macro_tss = mean(
      zone_metrics$tss,
      na.rm = TRUE
    )
  )
}

reference_overall_metrics <- rbindlist(
  reference_overall_list,
  fill = TRUE
)

reference_zone_metrics_common <- rbindlist(
  reference_zone_list,
  fill = TRUE
)

reference_overall_metrics[
  ,
  method_label := factor(
    method_labels[method],
    levels = rev(method_labels[reference_methods])
  )
]

macro_se_dt <- reference_zone_metrics_common[
  ,
  .(
    macro_balanced_accuracy_se = sd(balanced_accuracy, na.rm = TRUE) / sqrt(sum(is.finite(balanced_accuracy))),
    macro_recall_se = sd(recall, na.rm = TRUE) / sqrt(sum(is.finite(recall))),
    macro_specificity_se = sd(specificity, na.rm = TRUE) / sqrt(sum(is.finite(specificity))),
    macro_precision_se = sd(precision, na.rm = TRUE) / sqrt(sum(is.finite(precision))),
    macro_f1_se = sd(f1, na.rm = TRUE) / sqrt(sum(is.finite(f1))),
    macro_tss_se = sd(tss, na.rm = TRUE) / sqrt(sum(is.finite(tss)))
  ),
  by = method
]

reference_overall_metrics <- merge(
  reference_overall_metrics,
  macro_se_dt,
  by = "method",
  all.x = TRUE,
  sort = FALSE
)

fwrite(
  reference_overall_metrics,
  file.path(
    table_dir,
    "Figure_var_10b_reference_map_overall_metrics.csv"
  )
)

fwrite(
  reference_zone_metrics_common,
  file.path(
    table_dir,
    "Figure_var_10b_reference_map_zone_metrics_common_mask.csv"
  )
)

metric_labels_10b <- c(
  individual_map_coverage = "Individual map coverage",
  common_comparison_coverage = "Common comparison coverage",
  exact_zone_accuracy = "Exact-zone accuracy",
  broad_category_accuracy = "Broad-category accuracy",
  macro_balanced_accuracy = "Macro balanced accuracy",
  macro_recall = "Macro sensitivity",
  macro_specificity = "Macro specificity",
  macro_precision = "Macro precision",
  macro_f1 = "Macro F1",
  macro_tss = "Macro TSS"
)

metric_order_10b <- names(metric_labels_10b)

error_lookup_10b <- c(
  individual_map_coverage = NA_character_,
  common_comparison_coverage = NA_character_,
  exact_zone_accuracy = NA_character_,
  broad_category_accuracy = NA_character_,
  macro_balanced_accuracy = "macro_balanced_accuracy_se",
  macro_recall = "macro_recall_se",
  macro_specificity = "macro_specificity_se",
  macro_precision = "macro_precision_se",
  macro_f1 = "macro_f1_se",
  macro_tss = "macro_tss_se"
)

overall_long_10b <- rbindlist(
  lapply(metric_order_10b, function(metric_name) {
    se_col <- error_lookup_10b[[metric_name]]
    value_vec <- reference_overall_metrics[[metric_name]]
    se_vec <- if (is.na(se_col)) {
      rep(NA_real_, nrow(reference_overall_metrics))
    } else {
      reference_overall_metrics[[se_col]]
    }
    data.table(
      method = reference_overall_metrics$method,
      method_label = reference_overall_metrics$method_label,
      metric = metric_name,
      metric_label = metric_labels_10b[[metric_name]],
      value = value_vec,
      se = se_vec
    )
  }),
  fill = TRUE
)

overall_long_10b[
  ,
  `:=`(
    metric_label = factor(
      metric_label,
      levels = unname(metric_labels_10b[metric_order_10b])
    ),
    method_label = factor(
      method_label,
      levels = rev(method_labels[reference_methods])
    ),
    xmin = pmax(0, value - 1.96 * se),
    xmax = pmin(1, value + 1.96 * se)
  )
]

method_shapes_10b <- c(
  "Plain RF" = 16,
  "Plain MF RF" = 17,
  "Multiclass RF" = 15
)

figure_10b <- ggplot(
  overall_long_10b,
  aes(
    x = value,
    y = method_label,
    shape = method_label
  )
) +
  geom_errorbarh(
    data = overall_long_10b[is.finite(xmin) & is.finite(xmax)],
    aes(xmin = xmin, xmax = xmax),
    height = 0.16,
    linewidth = 0.5,
    color = "grey45"
  ) +
  geom_point(
    size = 2.9,
    color = "black"
  ) +
  facet_wrap(
    ~ metric_label,
    ncol = 2
  ) +
  scale_shape_manual(
    values = method_shapes_10b
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  labs(
    title = "Overall reference-map reconstruction performance",
    subtitle = "Compared on the same original raster; reconstruction agreement, not independent validation.",
    x = "Metric value",
    y = NULL,
    shape = "Model"
  ) +
  theme_bw(base_size = 10.2) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9.5)
  )

save_plot(
  figure_10b,
  file.path(
    figure_dir,
    "Figure_var_10b_reference_map_overall_metrics.png"
  ),
  10.4,
  9.6
)

# Keep the zone-level F1 bubble plot as a separate diagnostic figure.
bubble_dt <- copy(
  reference_zone_metrics_common[
    ,
    .(
      method,
      zone,
      f1
    )
  ]
)

ref_area <- cellSize(
  reference_map,
  unit = "km"
)

ref_modeled <- subst(
  reference_map,
  from = model_zoneID,
  to = model_zoneID,
  others = NA
)

area_dt <- as.data.table(
  zonal(
    ref_area,
    ref_modeled,
    fun = "sum",
    na.rm = TRUE
  )
)

setnames(
  area_dt,
  names(area_dt)[1:2],
  c(
    "zone",
    "area_km2"
  )
)

area_dt[
  ,
  `:=`(
    zone = as.integer(zone),
    area_km2 = as.numeric(area_km2)
  )
]

bubble_dt <- merge(
  bubble_dt,
  area_dt,
  by = "zone",
  all.x = TRUE,
  sort = FALSE
)

bubble_dt[
  ,
  `:=`(
    method_label = factor(
      method_labels[method],
      levels = method_labels[reference_methods]
    ),
    zone_factor = factor(
      as.character(zone),
      levels = rev(
        as.character(
          sort(unique(zone))
        )
      )
    )
  )
]

bubble_fill_values <- c(
  setNames(
    unname(method_colors[method_order]),
    method_labels[method_order]
  ),
  "Multiclass RF" = method_colors[["multiclass_rf"]]
)

figure_10c <- ggplot(
  bubble_dt,
  aes(
    x = f1,
    y = zone_factor,
    size = area_km2,
    fill = method_label
  )
) +
  geom_point(
    shape = 21,
    color = "grey20",
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = bubble_fill_values
  ) +
  facet_grid(
    . ~ method_label
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  labs(
    x = "Zone-level F1",
    y = "Vegetation zone",
    size = expression("Original area (km"^2*")"),
    fill = NULL,
    title = "Reference-map F1 by vegetation zone"
  ) +
  theme_bw(base_size = 9.5) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.spacing = grid::unit(2.6, "lines"),
    axis.text.y = element_text(size = 6.4),
    axis.text.x = element_text(size = 8.0),
    strip.text.x = element_text(face = "bold", size = 10)
  )

save_plot(
  figure_10c,
  file.path(
    figure_dir,
    "Figure_var_10c_reference_map_F1_bubble.png"
  ),
  15.5,
  7.2
)

# Keep the earlier zone-colored binary-vs-multiclass plot with 45-degree line.
multiclass_compare <- merge(
  bubble_dt[
    method %in% method_order,
    .(
      method,
      method_label,
      zone,
      area_km2,
      binary_f1 = f1
    )
  ],
  bubble_dt[
    method == "multiclass_rf",
    .(
      zone,
      multiclass_f1 = f1
    )
  ],
  by = "zone",
  all.x = TRUE,
  sort = FALSE
)

multiclass_compare[
  ,
  zone_factor := factor(
    as.character(zone),
    levels = as.character(model_zoneID)
  )
]

figure_10d <- ggplot(
  multiclass_compare,
  aes(
    x = multiclass_f1,
    y = binary_f1,
    size = area_km2,
    fill = zone_factor
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.55,
    color = "grey45"
  ) +
  geom_point(
    shape = 21,
    color = "grey20",
    alpha = 0.88,
    stroke = 0.22
  ) +
  facet_grid(
    . ~ method_label
  ) +
  scale_fill_manual(
    values = zone_color_vector(model_zoneID),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  labs(
    x = "Multiclass RF zone-level F1",
    y = "Binary-model zone-level F1",
    size = expression("Original area (km"^2*")"),
    title = "Zone-level F1: binary workflow vs Multiclass RF",
    subtitle = "Dashed line indicates equal F1."
  ) +
  theme_bw(base_size = 9.7) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.spacing = grid::unit(2.0, "lines"),
    strip.text.x = element_text(face = "bold", size = 10)
  )

save_plot(
  figure_10d,
  file.path(
    figure_dir,
    "Figure_var_10d_reference_map_F1_vs_multiclass_bubble.png"
  ),
  11.8,
  6.2
)


# 11. Reference Top-k agreement and future rank retention =====================

topk_analysis_dir <- file.path(
  assessment_dir,
  "future_topk_analysis"
)

reference_topk <- fread(
  require_file(
    file.path(
      topk_analysis_dir,
      "reference_map_topk_agreement.csv"
    )
  )
)[
  method %in% method_order
]

future_topk <- fread(
  require_file(
    file.path(
      topk_analysis_dir,
      "future_topk_retention_overall.csv"
    )
  )
)[
  method %in% method_order
]

required_rank_cutoffs <- c(
  1L,
  2L,
  3L,
  5L
)

expected_reference_topk <- CJ(
  method = method_order,
  rank_cutoff = required_rank_cutoffs,
  unique = TRUE
)

observed_reference_topk <- unique(
  reference_topk[
    ,
    .(
      method,
      rank_cutoff
    )
  ]
)

missing_reference_topk <- fsetdiff(
  expected_reference_topk,
  observed_reference_topk
)

if (nrow(missing_reference_topk) > 0) {
  print(missing_reference_topk)
  stop(
    "Incomplete reference Top-k results. Re-run script 8.2."
  )
}

expected_future_topk <- CJ(
  method = method_order,
  scenario = future_order,
  unique = TRUE
)

observed_future_topk <- unique(
  future_topk[
    ,
    .(
      method,
      scenario
    )
  ]
)

missing_future_topk <- fsetdiff(
  expected_future_topk,
  observed_future_topk
)

if (nrow(missing_future_topk) > 0) {
  print(missing_future_topk)
  stop(
    "Incomplete future Top-k results. Re-run script 8.2."
  )
}

# 11a. Reference-period Top-k agreement ----------------------------------------

reference_topk[
  ,
  `:=`(
    method_label = factor(
      method_labels[method],
      levels = method_labels[method_order]
    ),
    rank_label = factor(
      paste0(
        "Top-",
        rank_cutoff
      ),
      levels = paste0(
        "Top-",
        required_rank_cutoffs
      )
    )
  )
]

fwrite(
  reference_topk,
  file.path(
    table_dir,
    "Figure_var_11a_reference_topk_agreement.csv"
  )
)

figure_11a <- ggplot(
  reference_topk,
  aes(
    x = rank_cutoff,
    y = area_share,
    group = method_label,
    shape = method_label,
    linetype = method_label
  )
) +
  geom_line(
    linewidth = 0.8,
    color = "black"
  ) +
  geom_point(
    size = 2.5,
    fill = "white",
    color = "black",
    stroke = 0.55
  ) +
  scale_shape_manual(
    values = c(21, 24)
  ) +
  scale_linetype_manual(
    values = c("solid", "22")
  ) +
  scale_x_continuous(
    breaks = required_rank_cutoffs,
    labels = paste0(
      "Top-",
      required_rank_cutoffs
    )
  ) +
  scale_y_continuous(
    limits = c(0.5, 1),
    breaks = seq(
      0.5,
      1,
      by = 0.1
    ),
    labels = function(x) {
      paste0(
        round(
          100 * x
        ),
        "%"
      )
    }
  ) +
  labs(
    x = NULL,
    y = "Share of evaluated mapped area",
    shape = NULL,
    linetype = NULL,
    title = "Reference-period ecotype agreement by rank",
    subtitle = "Top-1 is exact assigned-zone agreement; higher ranks show whether the observed ecotype remains among the leading environmental analogues."
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

save_plot(
  figure_11a,
  file.path(
    figure_dir,
    "Figure_var_11a_reference_topk_agreement.png"
  ),
  8.5,
  5.0
)


# 11b. Future displacement severity --------------------------------------------

future_topk[
  ,
  `:=`(
    period = factor(
      period,
      levels = period_levels
    ),
    ssp = factor(
      ssp,
      levels = ssp_levels
    ),
    method_label = factor(
      method_labels[method],
      levels = method_labels[method_order]
    ),
    stable_top1 =
      stable_top1_area_share,
    changed_former_top3 =
      changed_existing_area_share *
      former_top3_among_changed_area,
    changed_former_rank4_5 =
      changed_existing_area_share *
      (
        former_top5_among_changed_area -
        former_top3_among_changed_area
      ),
    changed_former_below_top5 =
      changed_existing_area_share *
      (
        1 -
        former_top5_among_changed_area
      ),
    novel = novel_area_share
  )
]

future_topk_long <- melt(
  future_topk,
  id.vars = c(
    "method",
    "method_label",
    "scenario",
    "period",
    "ssp"
  ),
  measure.vars = c(
    "stable_top1",
    "changed_former_top3",
    "changed_former_rank4_5",
    "changed_former_below_top5",
    "novel"
  ),
  variable.name = "retention_class",
  value.name = "area_share"
)

future_topk_long[
  ,
  retention_label := factor(
    retention_class,
    levels = c(
      "stable_top1",
      "changed_former_top3",
      "changed_former_rank4_5",
      "changed_former_below_top5",
      "novel"
    ),
    labels = c(
      "Stable Top-1",
      "Changed; former in Top-3",
      "Changed; former rank 4-5",
      "Changed; former below Top-5",
      "Novel"
    )
  )
]

share_check <- future_topk_long[
  ,
  .(
    total_share = sum(
      area_share,
      na.rm = TRUE
    )
  ),
  by = .(
    method,
    scenario
  )
]

if (any(
  abs(
    share_check$total_share -
      1
  ) > 1e-6
)) {
  print(share_check)
  stop(
    "Future Top-k retention classes do not sum to one."
  )
}

fwrite(
  future_topk,
  file.path(
    table_dir,
    "Figure_var_11b_future_topk_retention_overall.csv"
  )
)

fwrite(
  future_topk_long,
  file.path(
    table_dir,
    "Figure_var_11b_future_topk_retention_long.csv"
  )
)

figure_11b <- ggplot(
  future_topk_long,
  aes(
    x = period,
    y = area_share,
    fill = retention_label
  )
) +
  geom_col(
    width = 0.58,
    color = "white",
    linewidth = 0.15
  ) +
  facet_grid(
    method_label ~ ssp
  ) +
  scale_fill_manual(
    values = c(
      "Stable Top-1" = "#5B8E7D",
      "Changed; former in Top-3" = "#D9C27E",
      "Changed; former rank 4-5" = "#D4A373",
      "Changed; former below Top-5" = "#B5654D",
      "Novel" = "#333333"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(
      0,
      1,
      by = 0.2
    ),
    labels = function(x) {
      paste0(
        round(
          100 * x
        ),
        "%"
      )
    }
  ) +
  labs(
    x = "Future period",
    y = "Share of mapped area",
    fill = NULL,
    title = "Future displacement of normal-period top-ranked analogues",
    subtitle = "Changed cells are separated by whether the former top-ranked ecotype remains within the future Top-3 or Top-5."
  ) +
  theme_bw(
    base_size = 10.5
  ) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = grid::unit(
      1.0,
      "lines"
    )
  )

save_plot(
  figure_11b,
  file.path(
    figure_dir,
    "Figure_var_11b_future_topk_retention.png"
  ),
  10.2,
  6.2
)


# 12. Multi-page species maps ====================================================

species_names <- sort(
  unique(
    assigned_species_area$Species
  )
)

number_of_pages <- ceiling(
  length(species_names) /
    species_per_page
)

for (method in method_order) {
  for (scenario in future_order) {
    for (page_index in seq_len(
      number_of_pages
    )) {
      
      start_index <-
        (page_index - 1L) *
        species_per_page +
        1L
      
      end_index <- min(
        page_index *
          species_per_page,
        length(species_names)
      )
      
      page_species <- species_names[
        start_index:end_index
      ]
      
      assigned_files <- file.path(
        assigned_result_root,
        method,
        scenario,
        "species niche",
        paste0(
          page_species,
          "_species_niche.tif"
        )
      )
      
      dual_files <- file.path(
        dual_result_root,
        method,
        scenario,
        "species niche binary",
        paste0(
          safe_name(page_species),
          "_species_dual_binary_threshold",
          new_threshold,
          ".tif"
        )
      )
      
      plot_binary_species_page(
        page_species,
        assigned_files,
        paste0(
          "Assigned-zone species niches | ",
          method_labels[method],
          " | ",
          scenario,
          " | page ",
          page_index
        ),
        file.path(
          assigned_species_page_dir,
          paste0(
            "assigned_species_",
            method,
            "_",
            scenario,
            "_page",
            page_index,
            ".png"
          )
        ),
        reference_mask,
        fill_color = "#2E8B57"
      )
      
      plot_binary_species_page(
        page_species,
        dual_files,
        paste0(
          "Dual-suitability species niches >= 0.4 | ",
          method_labels[method],
          " | ",
          scenario,
          " | page ",
          page_index
        ),
        file.path(
          dual_species_page_dir,
          paste0(
            "dual_species_",
            method,
            "_",
            scenario,
            "_page",
            page_index,
            ".png"
          )
        ),
        reference_mask,
        fill_color = "#2E5E8C"
      )
    }
  }
}


# 13. Save figure-source tables ==================================================

fwrite(
  performance_long,
  file.path(
    table_dir,
    "figure_climate_soil_zone_metrics_long.csv"
  )
)

fwrite(
  map_zone_long,
  file.path(
    table_dir,
    "figure_normal_map_assessment_data_var.csv"
  )
)

fwrite(
  novel_area,
  file.path(
    table_dir,
    "figure_novel_area_data_var.csv"
  )
)

fwrite(
  transition_summary,
  file.path(
    table_dir,
    "figure_transition_share_data_var.csv"
  )
)

fwrite(
  assigned_species_area,
  file.path(
    table_dir,
    "figure_assigned_species_area_data_var.csv"
  )
)

fwrite(
  dual_species_long,
  file.path(
    table_dir,
    "figure_dual_species_area_data_var.csv"
  )
)

fwrite(
  dual_population_area,
  file.path(
    table_dir,
    "figure_dual_population_area_data_var.csv"
  )
)

if (exists("ranking_margin_long")) {
  fwrite(
    ranking_margin_long,
    file.path(
      table_dir,
      "figure_ranking_uncertainty_data_var.csv"
    )
  )
}

cat(
  "\nCOMPLETE\n",
  "Threshold: ",
  new_threshold,
  "\n",
  "Main-text models: Plain RF, Plain MF RF\n",
  "Assessment source: final 5.1 selected-variable climate + soil Plain RF / Plain MF RF models\n",
  "Climate assessment source: rf_var / mf_var independent-test results\n",
  "Figures: ",
  figure_dir,
  "\n",
  "Tables: ",
  table_dir,
  "\n",
  sep = ""
)
