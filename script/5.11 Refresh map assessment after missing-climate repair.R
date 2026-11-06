# Refresh ONLY section 3 (assigned-map assessment) of script 5.1.
# Reads corrected reference-period assigned maps; no fitted model is loaded,
# no held-out prediction/test is repeated, and no raster is written.
# Pixel-based metrics and output schemas deliberately match the existing 5.1.
local({
  refresh_map_assessment <- function() {
    for (pkg in c("terra", "data.table")) {
      if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing R package: ", pkg)
    }
    library(terra)
    library(data.table)
    root <- normalizePath(Sys.getenv("ECOCHINA2_DIR", unset = getwd()),
                          winslash = "/", mustWork = TRUE)
    while (!file.exists(file.path(root, "script", "1. data.R"))) {
      parent <- dirname(root)
      if (parent == root) stop("Run inside ecoChina2 or set ECOCHINA2_DIR.")
      root <- parent
    }
    assessment_dir <- Sys.getenv("ECOCHINA2_MAP_ASSESSMENT_OUT",
      unset = file.path(root, "assessment_var"))
    reference_file <- file.path(root, "raster", "ecosys_ori.tif")
    reference <- rast(reference_file)
    zoneID <- c(1:7, 9:50, 52:55)
    method_order <- c("rf_var", "mf_var")
    method_labels <- c(rf_var = "Plain RF", mf_var = "Plain MF RF")
    prediction_files <- file.path(root, "result maps", method_order,
      "assigned_zone_normal_threshold0.4_tol1e-04_novel99_maskNA8_noNovelNormal.tif")
    if (!all(file.exists(prediction_files))) {
      stop("Missing corrected normal map(s):\n",
        paste(prediction_files[!file.exists(prediction_files)], collapse = "\n"))
    }
    for (file in prediction_files) {
      if (!isTRUE(compareGeom(reference, rast(file), stopOnError = FALSE))) {
        stop("Geometry mismatch for normal map: ", file)
      }
    }
    div <- function(a, b) ifelse(is.finite(b) & b > 0, a / b, NA_real_)
    mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
    nz <- length(zoneID)
    counts <- lapply(method_order, function(method) matrix(0, nz, nz))
    names(counts) <- method_order
    compared_pixels <- setNames(numeric(length(method_order)), method_order)
    valid_original_pixels <- 0
    input <- c(reference, rast(prediction_files))
    rows_per_block <- max(1L, min(64L, floor(2e6 / ncol(reference))))
    starts <- seq.int(1L, nrow(reference), by = rows_per_block)
    readStart(input)
    input_open <- TRUE
    on.exit(if (input_open) readStop(input), add = TRUE)
    for (b in seq_along(starts)) {
      start <- starts[b]
      nr <- min(rows_per_block, nrow(reference) - start + 1L)
      values <- readValues(input, row = start, nrows = nr, mat = TRUE)
      original_index <- match(values[, 1], zoneID, nomatch = 0L)
      valid_original_pixels <- valid_original_pixels + sum(original_index > 0L)
      for (m in seq_along(method_order)) {
        predicted_index <- match(values[, m + 1L], zoneID, nomatch = 0L)
        good <- original_index > 0L & predicted_index > 0L
        compared_pixels[m] <- compared_pixels[m] + sum(good)
        key <- original_index[good] + nz * (predicted_index[good] - 1L)
        counts[[m]] <- counts[[m]] + matrix(tabulate(key, nbins = nz * nz), nz, nz)
      }
      if (b %% 10L == 0L || b == length(starts)) {
        cat(sprintf("[MAP ASSESSMENT] %.0f%%\n", 100 * b / length(starts)))
        flush.console()
      }
    }
    readStop(input); input_open <- FALSE
    if (any(compared_pixels <= 0)) stop("A normal map has no valid modeled-zone comparison cells.")
    dir.create(assessment_dir, recursive = TRUE, showWarnings = FALSE)
    map_confusion_results <- map_zone_results <- map_overall_results <- list()
    for (m in seq_along(method_order)) {
      method_key <- method_order[m]
      matrix_counts <- counts[[m]]
      present <- which(matrix_counts > 0, arr.ind = TRUE)
      confusion <- data.table(original_zone = as.integer(zoneID[present[, 1]]),
        predicted_zone = as.integer(zoneID[present[, 2]]),
        n = as.numeric(matrix_counts[present]), method = method_key,
        method_label = method_labels[[method_key]])
      setorder(confusion, original_zone, predicted_zone)
      total <- sum(confusion$n)
      if (total != compared_pixels[m]) stop("Confusion counts do not reconcile: ", method_key)
      zone_metrics <- rbindlist(lapply(zoneID, function(zone_value) {
        TP <- confusion[original_zone == zone_value & predicted_zone == zone_value, sum(n)]
        FN <- confusion[original_zone == zone_value & predicted_zone != zone_value, sum(n)]
        FP <- confusion[original_zone != zone_value & predicted_zone == zone_value, sum(n)]
        TN <- total - TP - FN - FP
        recall <- div(TP, TP + FN)
        specificity <- div(TN, TN + FP)
        precision <- div(TP, TP + FP)
        data.table(method = method_key, method_label = method_labels[[method_key]],
          zone = zone_value, original_pixels = TP + FN, predicted_pixels = TP + FP,
          TP = TP, TN = TN, FP = FP, FN = FN, accuracy = div(TP + TN, total),
          balanced_accuracy = div(recall + specificity, 2), recall = recall,
          specificity = specificity, precision = precision,
          f1 = div(2 * precision * recall, precision + recall),
          tss = recall + specificity - 1)
      }))
      overall <- data.table(method = method_key, method_label = method_labels[[method_key]],
        valid_original_pixels = valid_original_pixels, compared_pixels = compared_pixels[m],
        missing_predictions = valid_original_pixels - compared_pixels[m],
        coverage = div(compared_pixels[m], valid_original_pixels),
        exact_accuracy = confusion[original_zone == predicted_zone, sum(n)] / total,
        macro_balanced_accuracy = mean_na(zone_metrics$balanced_accuracy),
        macro_recall = mean_na(zone_metrics$recall),
        macro_specificity = mean_na(zone_metrics$specificity),
        macro_precision = mean_na(zone_metrics$precision), macro_f1 = mean_na(zone_metrics$f1),
        macro_tss = mean_na(zone_metrics$tss), source_map = prediction_files[m])
      map_confusion_results[[method_key]] <- confusion
      map_zone_results[[method_key]] <- zone_metrics
      map_overall_results[[method_key]] <- overall
      fwrite(dcast(confusion, original_zone ~ predicted_zone, value.var = "n", fill = 0),
        file.path(assessment_dir, paste0("normal_map_confusion_matrix_", method_key, ".csv")))
    }
    normal_map_confusion <- rbindlist(map_confusion_results, fill = TRUE)
    normal_map_zone_metrics <- rbindlist(map_zone_results, fill = TRUE)
    normal_map_overall_metrics <- rbindlist(map_overall_results, fill = TRUE)
    setorder(normal_map_confusion, method, original_zone, predicted_zone)
    setorder(normal_map_zone_metrics, method, zone)
    setorder(normal_map_overall_metrics, method)
    fwrite(normal_map_confusion, file.path(assessment_dir, "normal_map_confusion_long_var.csv"))
    fwrite(normal_map_zone_metrics, file.path(assessment_dir, "normal_map_zone_metrics_var.csv"))
    fwrite(normal_map_overall_metrics, file.path(assessment_dir, "normal_map_overall_metrics_var.csv"))
    errors <- normal_map_confusion[original_zone != predicted_zone]
    errors_from <- errors[, .(error_pixels = sum(n)),
      by = .(method, method_label, original_zone)][order(method, -error_pixels)]
    errors_into <- errors[, .(error_pixels = sum(n)),
      by = .(method, method_label, predicted_zone)][order(method, -error_pixels)]
    fwrite(errors_from, file.path(assessment_dir, "normal_map_errors_from_original_zone_var.csv"))
    fwrite(errors_into, file.path(assessment_dir, "normal_map_errors_into_assigned_zone_var.csv"))
    cat("\n[NORMAL MAP ASSESSMENT COMPLETE]\n")
    print(normal_map_overall_metrics)
    cat("Held-out test/model files and project rasters were not modified.\n")
  }
  refresh_map_assessment()
})
