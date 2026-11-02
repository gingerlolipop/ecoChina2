# Diagnose missing reference-period climate suitability at a few saved audit cells.
# No training, national prediction, gate-sensitivity rerun, or result overwrite.
# Copy to script/ and source("script/trace_climate_na_points.R") from ecoChina2.

local({
  run_trace <- function() {
    for (pkg in c("terra", "data.table", "randomForest")) {
      if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing R package: ", pkg)
    }
    library(terra)
    library(data.table)
    find_root <- function() {
      p <- Sys.getenv("ECOCHINA2_DIR", unset = getwd())
      p <- normalizePath(p, winslash = "/", mustWork = TRUE)
      repeat {
        if (file.exists(file.path(p, "script", "1. data.R"))) return(p)
        q <- dirname(p)
        if (q == p) stop("Run inside ecoChina2 or set ECOCHINA2_DIR.")
        p <- q
      }
    }
    root <- find_root()
    audit_dir <- Sys.getenv("ECOCHINA2_TRACE_AUDIT_DIR", unset = file.path(
      root, "assessment", "soil_domain_audit_20261101_135508"))
    mask_file <- file.path(audit_dir, "domain_normal.tif")
    reference_file <- file.path(root, "raster", "ecosys_ori.tif")
    if (!all(file.exists(c(mask_file, reference_file)))) {
      stop("Need the existing local audit domain_normal.tif and raster/ecosys_ori.tif.\n",
           "Set ECOCHINA2_TRACE_AUDIT_DIR if the audit directory differs.")
    }
    out <- Sys.getenv("ECOCHINA2_TRACE_OUT", unset = file.path(root, "assessment",
      paste0("climate_na_trace_", format(Sys.time(), "%Y%m%d_%H%M%S"))))
    if (dir.exists(out)) stop("Choose a new ECOCHINA2_TRACE_OUT: ", out)
    dir.create(out, recursive = TRUE)
    tmp <- file.path(out, "temporary")
    dir.create(tmp)
    old_options <- terraOptions(print = FALSE)
    terraOptions(tempdir = tmp)
    on.exit({
      restore_tmp <- old_options$tempdir
      if (is.null(restore_tmp) || !dir.exists(restore_tmp)) restore_tmp <- tempdir()
      terraOptions(tempdir = restore_tmp)
      unlink(tmp, recursive = TRUE)
    }, add = TRUE)

    domain <- rast(mask_file)
    reference <- rast(reference_file)
    if (!compareGeom(domain, reference, stopOnError = FALSE)) {
      stop("Audit mask and reference geometry differ; no automatic repair attempted.")
    }
    cat("[TRACE 1/3] Read one saved integer mask, not the suitability stacks.\n")
    readStart(domain)
    ids <- list()
    starts <- seq.int(1L, nrow(domain), by = 128L)
    for (j in seq_along(starts)) {
      nr <- min(128L, nrow(domain) - starts[j] + 1L)
      v <- readValues(domain, row = starts[j], nrows = nr, mat = FALSE)
      hit <- which(v == 3L)
      if (length(hit)) ids[[length(ids) + 1L]] <- (starts[j] - 1) * ncol(domain) + hit
      if (j %% 10L == 0L || j == length(starts)) {
        cat(sprintf("  Mask scan: %.0f%%\n", 100 * j / length(starts)))
        flush.console()
      }
    }
    readStop(domain)
    ids <- unlist(ids, use.names = FALSE)
    if (!length(ids)) stop("No code-3 cells found in the existing mask.")
    xy <- xyFromCell(reference, ids)
    ll <- if (is.lonlat(reference)) xy else terra::project(
      xy, from = crs(reference), to = "EPSG:4326")
    candidates <- data.table(cell = ids, x = xy[, 1], y = xy[, 2],
      longitude = ll[, 1], latitude = ll[, 2])
    candidates[, grid := paste(floor(longitude), floor(latitude), sep = "_")]
    bins <- candidates[, .N, by = grid][order(-N)]
    n_problem <- as.integer(Sys.getenv("ECOCHINA2_TRACE_POINTS", "8"))
    if (is.na(n_problem) || n_problem < 1L || n_problem > 32L) {
      stop("ECOCHINA2_TRACE_POINTS must be between 1 and 32.")
    }
    chosen <- rbindlist(lapply(head(bins$grid, n_problem), function(key) {
      z <- candidates[grid == key]
      z[ceiling(nrow(z) / 2)]
    }))
    chosen[, point_type := "gate020_only_problem"]
    # Nearby controls are selected from a sparse grid of offsets, not asserted
    # to be the mathematically nearest valid cells.
    offsets <- c(-128, -64, -32, -16, -8, -4, -2, -1, 0, 1, 2, 4, 8, 16, 32, 64, 128)
    delta <- as.matrix(expand.grid(dx = offsets, dy = offsets))
    controls <- list()
    for (j in seq_len(nrow(chosen))) {
      test_xy <- cbind(chosen$x[j] + delta[, 1] * res(reference)[1],
                       chosen$y[j] + delta[, 2] * res(reference)[2])
      codes <- terra::extract(domain, test_xy)[[1]]
      good <- which(!is.na(codes) & codes == 1L)
      if (length(good)) {
        k <- good[which.min(rowSums(delta[good, , drop = FALSE]^2))]
        cc <- cellFromXY(reference, test_xy[k, , drop = FALSE])
        cxy <- xyFromCell(reference, cc)
        cll <- if (is.lonlat(reference)) cxy else terra::project(
          cxy, from = crs(reference), to = "EPSG:4326")
        controls[[length(controls) + 1L]] <- data.table(cell = cc,
          x = cxy[, 1], y = cxy[, 2], longitude = cll[, 1], latitude = cll[, 2],
          grid = NA_character_, point_type = "nearby_common_valid_control")
      }
    }
    points <- unique(rbindlist(c(list(chosen), controls)), by = "cell")
    points[, sample_id := seq_len(.N)]
    pv <- vect(as.data.frame(points[, .(x, y)]), geom = c("x", "y"), crs = crs(reference))
    points[, reference_zone := terra::extract(reference, pv, ID = FALSE)[[1]]]
    fwrite(points, file.path(out, "sample_points.csv"))
    writeLines(c(paste("All code-3 cells in existing mask:", length(ids)),
      paste("Problem cells sampled:", nrow(chosen)),
      "A sample is used for diagnosis, not for estimating national frequencies.",
      "Raw-directory provenance is not inferred from file existence.",
      capture.output(sessionInfo())), file.path(out, "run_info.txt"))

    cat("[TRACE 2/3] Inspect fitted variable lists and small raw-data windows.\n")
    model_zones <- as.integer(strsplit(Sys.getenv("ECOCHINA2_TRACE_MODEL_ZONES", "15,30"), ",")[[1]])
    if (anyNA(model_zones) || length(model_zones) > 4L) stop("Choose at most four model zone IDs.")
    climate_suit_root <- Sys.getenv("ECOCHINA2_AUDIT_CLIMATE_SUITABILITY_DIR",
      unset = file.path(root, "clim suitability", "mf_var"))
    soil_suit_root <- Sys.getenv("ECOCHINA2_AUDIT_SOIL_SUITABILITY_DIR",
      unset = file.path(root, "soil suitability", "plain_mf", "normal"))
    add_metadata <- function(files, kind, model_zone, source_id = NA_integer_) {
      info <- file.info(files)
      z <- data.table(model_zone = model_zone, source_id = source_id, kind = kind,
        path = files, exists = file.exists(files), size_bytes = info$size,
        modified = as.character(info$mtime))
      dest <- file.path(out, "input_files.csv")
      fwrite(z, dest, append = file.exists(dest), col.names = !file.exists(dest))
    }
    record_geometry <- function(x, label, model_zone = NA_integer_, source_id = NA_integer_) {
      z <- data.table(dataset = label, model_zone = model_zone, source_id = source_id,
        rows = nrow(x), columns = ncol(x), xmin = xmin(x), xmax = xmax(x),
        ymin = ymin(x), ymax = ymax(x), x_resolution = res(x)[1], y_resolution = res(x)[2],
        x_origin = origin(x)[1], y_origin = origin(x)[2], crs = crs(x))
      dest <- file.path(out, "raster_geometry.csv")
      fwrite(z, dest, append = file.exists(dest), col.names = !file.exists(dest))
    }
    record_geometry(reference, "reference_ecotype")
    get_point_values <- function(file) {
      if (!file.exists(file)) return(rep(NA_real_, nrow(points)))
      z <- rast(file)
      terra::extract(z, terra::project(pv, crs(z)), ID = FALSE)[[1]]
    }
    results_file <- file.path(out, "point_pipeline_check.csv")
    dirs_checked <- list()
    for (zone in model_zones) {
      model_file <- file.path(root, "rf", paste0("clm_mfVar_zone", zone, ".Rdata"))
      if (!file.exists(model_file)) stop("Missing fitted model: ", model_file)
      env <- new.env(parent = baseenv())
      load(model_file, envir = env)
      if (!exists("clm_mfVar", envir = env, inherits = FALSE)) stop("Missing clm_mfVar object.")
      model <- get("clm_mfVar", envir = env)
      vars <- model$varlist
      if (is.null(vars) || !length(vars)) vars <- rownames(model$importance)
      if (is.null(vars) || !length(vars)) stop("No fitted variable names in model: ", zone)
      vars <- as.character(vars)
      fwrite(data.table(variable = vars), file.path(out, paste0("model_zone", zone, "_variables.csv")))
      cat("  Model zone", zone, "uses", length(vars), "variables.\n")
      climate_file <- file.path(climate_suit_root, "normal", paste0("clim_suit_zone", zone, ".tif"))
      soil_file <- file.path(soil_suit_root, paste0("soil_suit_zone", zone, ".tif"))
      dual_file <- file.path(root, "dual suit", "mf_var", "normal",
                             paste0("dual_suitability_zone", zone, ".tif"))
      if (!file.exists(climate_file)) stop("Missing stored climate suitability: ", climate_file)
      record_geometry(rast(climate_file), "stored_climate_suitability", zone)
      saved_climate <- get_point_values(climate_file)
      saved_soil <- get_point_values(soil_file)
      saved_dual <- get_point_values(dual_file)
      add_metadata(c(model_file, climate_file, soil_file, dual_file), "model_or_saved_surface", zone)
      env_raw <- Sys.getenv("ECOCHINA2_CLIMATE_DIR", unset = "")
      paths <- unique(c(if (nzchar(env_raw)) c(env_raw, file.path(env_raw, "Normal_1961_1990")),
        file.path(dirname(root), "ecoChina", "play", "China", "ClimateData", "CN", "800m", "Normal_1961_1990"),
        file.path(root, "data", "rasters", "climate", "Normal_1961_1990")))
      full <- vapply(paths, function(p) all(file.exists(file.path(p, paste0(vars, ".tif")))), logical(1))
      directory_check <- data.table(model_zone = zone, directory = paths,
        exists = dir.exists(paths), contains_all_selected_variables = full)
      dirs_checked[[length(dirs_checked) + 1L]] <- directory_check
      fwrite(rbindlist(dirs_checked), file.path(out, "raw_directory_candidates.csv"))
      if (!any(full)) stop("No raw-climate directory contains all fitted variables.\n",
        "Set ECOCHINA2_CLIMATE_DIR to the original raw-climate directory, then rerun.\n",
        "The directories checked are recorded in: ", out)
      raw_dirs <- unique(normalizePath(paths[full], winslash = "/", mustWork = TRUE))
      for (source_id in seq_along(raw_dirs)) {
        raw_files <- file.path(raw_dirs[source_id], paste0(vars, ".tif"))
        add_metadata(raw_files, "raw_climate", zone, source_id)
        raw <- rast(raw_files)
        names(raw) <- vars
        record_geometry(raw, "raw_climate", zone, source_id)
        p_raw <- terra::project(pv, crs(raw))
        raw_xy <- crds(p_raw)
        raw_values <- as.data.frame(terra::extract(raw, p_raw, ID = FALSE))
        fwrite(cbind(data.table(sample_id = points$sample_id), as.data.table(raw_values)),
          file.path(out, paste0("raw_values_zone", zone, "_source", source_id, ".csv")))
        metadata <- data.table(variable = vars, file = raw_files,
          nonfinite_at_problem_points = vapply(raw_values, function(z) sum(
            !is.finite(z[points$point_type == "gate020_only_problem"])), integer(1)),
          nonfinite_at_control_points = vapply(raw_values, function(z) sum(
            !is.finite(z[points$point_type == "nearby_common_valid_control"])), integer(1)))
        metadata[, nodata_flag := vapply(seq_along(vars), function(k) {
          paste(terra::NAflag(raw[[k]]), collapse = ";")
        }, character(1))]
        fwrite(metadata, file.path(out, paste0("raw_variable_missingness_zone", zone,
                                              "_source", source_id, ".csv")))
        same_crs <- compareGeom(raw[[1]], reference, ext = FALSE, rowcol = FALSE,
                               res = FALSE, stopOnError = FALSE)
        same_grid <- compareGeom(raw[[1]], reference, stopOnError = FALSE)
        rxy <- res(raw)
        for (j in seq_len(nrow(points))) {
          cat("  Zone", zone, "raw source", source_id, "point", j, "/", nrow(points), "\n")
          flush.console()
          native_cell <- cellFromXY(raw, raw_xy[j, , drop = FALSE])
          aligned_value <- native_prediction <- NA_real_
          window_complete <- window_total <- 0L
          error_text <- ""
          if (is.finite(native_cell) && same_crs) {
            tryCatch({
              rr <- rowFromCell(raw, native_cell)
              cc <- colFromCell(raw, native_cell)
              rmin <- max(1L, rr - 3L); rmax <- min(nrow(raw), rr + 3L)
              cmin <- max(1L, cc - 3L); cmax <- min(ncol(raw), cc + 3L)
              box <- ext(xmin(raw) + (cmin - 1) * rxy[1], xmin(raw) + cmax * rxy[1],
                         ymax(raw) - rmax * rxy[2], ymax(raw) - (rmin - 1) * rxy[2])
              win <- crop(raw, box)
              w <- values(win, mat = TRUE)
              complete <- rowSums(is.finite(w)) == length(vars)
              window_total <- nrow(w); window_complete <- sum(complete)
              pred <- rep(NA_real_, nrow(w))
              if (any(complete)) {
                input <- as.data.frame(w[complete, , drop = FALSE]); names(input) <- vars
                pp <- stats::predict(model, input, type = "prob")
                if (!"1" %in% colnames(pp)) stop("Fitted model has no class '1'.")
                pred[complete] <- pp[, "1"]
              }
              local_pred <- rast(win[[1]]); values(local_pred) <- pred
              native_prediction <- terra::extract(local_pred, raw_xy[j, , drop = FALSE])[[1]]
              if (same_grid) {
                aligned_value <- native_prediction
              } else {
                rc <- rowFromCell(reference, points$cell[j])
                cl <- colFromCell(reference, points$cell[j])
                target <- rast(nrows = 1L, ncols = 1L,
                  xmin = xmin(reference) + (cl - 1) * res(reference)[1],
                  xmax = xmin(reference) + cl * res(reference)[1],
                  ymin = ymax(reference) - rc * res(reference)[2],
                  ymax = ymax(reference) - (rc - 1) * res(reference)[2], crs = crs(reference))
                aligned_value <- values(resample(local_pred, target, method = "bilinear"))[1]
              }
            }, error = function(e) error_text <<- conditionMessage(e))
          }
          missing_vars <- vars[!is.finite(as.numeric(unlist(
            raw_values[j, , drop = FALSE], use.names = FALSE)))]
          correct_dual <- if (is.finite(saved_climate[j]) && is.finite(saved_soil[j])) {
            if (saved_soil[j] > 0.20) saved_climate[j] else 0
          } else NA_real_
          diagnosis <- if (!same_crs) "raw_and_reference_CRS_differ" else if (nzchar(error_text)) {
            "local_recompute_error"
          } else if (!is.finite(saved_climate[j]) && is.finite(aligned_value)) {
            "stored_surface_disagrees_with_local_recompute"
          } else if (is.finite(saved_climate[j]) && !is.finite(aligned_value)) {
            "stored_climate_finite_but_local_recompute_missing"
          } else if (is.finite(saved_climate[j]) && is.finite(aligned_value) &&
                     abs(saved_climate[j] - aligned_value) > 1e-5) {
            "stored_and_local_finite_predictions_differ"
          } else if (length(missing_vars)) "raw_selected_variable_missing_at_native_cell" else if (
            !is.finite(saved_climate[j])) "alignment_mask_or_cache_requires_review" else "saved_climate_finite"
          result <- data.table(model_zone = zone, source_id = source_id, raw_directory = raw_dirs[source_id],
            sample_id = points$sample_id[j], point_type = points$point_type[j],
            cell = points$cell[j], longitude = points$longitude[j], latitude = points$latitude[j],
            reference_zone = points$reference_zone[j], n_selected_variables = length(vars),
            raw_missing_variables = paste(missing_vars, collapse = ";"),
            raw_and_reference_same_crs = same_crs, raw_and_reference_same_grid = same_grid,
            complete_cells_in_local_window = window_complete, total_cells_in_local_window = window_total,
            stored_climate_suitability = saved_climate[j], stored_soil_suitability = saved_soil[j],
            stored_dual_suitability = saved_dual[j], strict_dual_from_stored_inputs = correct_dual,
            local_native_prediction = native_prediction, local_aligned_prediction = aligned_value,
            absolute_stored_vs_local_difference = if (is.finite(saved_climate[j]) &&
              is.finite(aligned_value)) abs(saved_climate[j] - aligned_value) else NA_real_,
            diagnosis = diagnosis, error = error_text)
          fwrite(result, results_file, append = file.exists(results_file),
                 col.names = !file.exists(results_file))
        }
        rm(raw)
        gc()
      }
      rm(model, env)
      gc()
    }
    cat("[TRACE 3/3 SAVED]", out, "\n")
    cat("This is a small reference-period diagnostic; no project prediction or model was overwritten.\n")
  }
  run_trace()
})
