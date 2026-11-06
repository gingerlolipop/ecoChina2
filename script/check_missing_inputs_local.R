# EcoChina2: local, read-only missing-input check.
# Save this file in H:/Jing/ecoChina2/script/, then run:
#   setwd("H:/Jing/ecoChina2")
#   source("script/check_missing_inputs_local.R")
# Requires the already installed terra and data.table packages.
# No fitting, national prediction, interpolation/filling, or input replacement.
# Writes small reports only under a NEW assessment/local_input_check_* folder.
# Optional settings (set BEFORE source):
#   Sys.setenv(ECOCHINA2_CLIMATE_DIR=".../climate")
#   Sys.setenv(ECOCHINA2_SOIL_RASTER_DIR=".../soil TIFF folder")
#   Sys.setenv(ECOCHINA2_CHECK_STATUS_FILE=".../availability_mf_var_normal.tif")
#   Sys.setenv(ECOCHINA2_CHECK_PER_ZONE="3")  # 1..20, for each missing type
#   Sys.setenv(ECOCHINA2_CHECK_DOWNLOAD_REFERENCE=".../original_download_reference.tif")

local({
  run_check <- function() {
    for (p in c("terra", "data.table")) {
      if (!requireNamespace(p, quietly = TRUE)) stop("Missing installed package: ", p)
    }
    library(terra)
    library(data.table)

    first_file <- function(paths) {
      paths <- unique(paths[!is.na(paths) & nzchar(paths)])
      hit <- paths[file.exists(paths) & !dir.exists(paths)]
      if (length(hit)) hit[1] else NA_character_
    }
    root <- normalizePath(Sys.getenv("ECOCHINA2_DIR", unset = getwd()),
                          winslash = "/", mustWork = TRUE)
    repeat {
      if (file.exists(file.path(root, "script", "1. data.R"))) break
      up <- dirname(root)
      if (up == root) stop("Run from ecoChina2 or set ECOCHINA2_DIR.")
      root <- up
    }
    n_sample <- suppressWarnings(as.integer(Sys.getenv("ECOCHINA2_CHECK_PER_ZONE", "3")))
    if (length(n_sample) != 1L || is.na(n_sample) || n_sample < 1L || n_sample > 20L) {
      stop("ECOCHINA2_CHECK_PER_ZONE must be an integer from 1 through 20.")
    }
    modeled_zones <- c(1:7, 9:50, 52:55)
    reference_file <- file.path(root, "raster", "ecosys_ori.tif")
    status_file <- first_file(c(Sys.getenv("ECOCHINA2_CHECK_STATUS_FILE", ""),
      file.path(root, "assessment", "missing_climate_repair_v1",
                "availability_mf_var_normal.tif")))
    if (!file.exists(reference_file)) stop("Missing reference: ", reference_file)
    if (is.na(status_file)) {
      stop("The existing repaired availability map was not found.\n",
        "Do NOT rerun prediction or training for this check. Locate the existing\n",
        "availability_mf_var_normal.tif and set ECOCHINA2_CHECK_STATUS_FILE to its path.")
    }
    reference <- rast(reference_file)
    status <- rast(status_file)[[1]]
    if (!compareGeom(reference, status, stopOnError = FALSE)) {
      stop("Availability/reference geometry mismatch; no automatic resampling attempted.")
    }
    out_base <- file.path(root, "assessment", paste0("local_input_check_",
                                                   format(Sys.time(), "%Y%m%d_%H%M%S")))
    out <- out_base
    suffix <- 0L
    while (dir.exists(out)) {
      suffix <- suffix + 1L
      out <- paste0(out_base, "_", suffix)
    }
    dir.create(out, recursive = TRUE)
    log_file <- file.path(out, "00_summary.txt")
    log <- function(...) {
      line <- paste0(...)
      cat(line, "\n")
      cat(line, "\n", file = log_file, append = TRUE)
      flush.console()
    }
    log("READ-ONLY LOCAL INPUT CHECK")
    log("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
    log("No fitted models or existing outputs are modified.")
    log("Reference: ", reference_file)
    log("Status map: ", status_file)
    writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))

    errors <- list()
    record_error <- function(stage, path, msg) {
      errors[[length(errors) + 1L]] <<- data.table(stage = stage, path = path, error = msg)
    }
    geometry <- list()
    manifests <- list()
    record_file <- function(path, kind) {
      z <- file.info(path)
      manifests[[length(manifests) + 1L]] <<- data.table(kind = kind, path = path,
        size_bytes = z$size, modified = as.character(z$mtime))
    }
    same_crs <- function(a, b) compareGeom(a[[1]], b[[1]], crs = TRUE,
      ext = FALSE, rowcol = FALSE, res = FALSE, stopOnError = FALSE)
    record_geometry <- function(r, path, kind) {
      geometry[[length(geometry) + 1L]] <<- data.table(kind = kind, path = path,
        rows = nrow(r), columns = ncol(r), xmin = xmin(r), xmax = xmax(r),
        ymin = ymin(r), ymax = ymax(r), x_resolution = res(r)[1], y_resolution = res(r)[2],
        x_origin = origin(r)[1], y_origin = origin(r)[2],
        same_crs_as_vegetation = same_crs(r, reference),
        same_geometry_as_vegetation = compareGeom(r[[1]], reference[[1]], stopOnError = FALSE),
        nodata_flag = paste(NAflag(r), collapse = ";"), crs = crs(r))
      record_file(path, kind)
    }
    record_geometry(reference, reference_file, "vegetation")
    record_geometry(status, status_file, "availability")

    # One bounded-memory scan of two EXISTING layers, not 53 suitability layers.
    # Uniform top-priority sampling within each modeled zone and availability type.
    # Kept priorities are generated independently for each cell; only the best k persist.
    counts <- list()
    reservoir <- new.env(parent = emptyenv())
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) saved_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
      if (had_seed) assign(".Random.seed", saved_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(list = ".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(20261103)
    input <- c(reference, status)
    readStart(input)
    reading <- TRUE
    on.exit(if (reading) try(readStop(input), silent = TRUE), add = TRUE)
    starts <- seq.int(1L, nrow(input), by = 128L)
    log("[1/3] Scan existing reference/status layers and sample by zone.")
    for (j in seq_along(starts)) {
      nr <- min(128L, nrow(input) - starts[j] + 1L)
      v <- readValues(input, row = starts[j], nrows = nr, mat = TRUE)
      keep <- which(is.finite(v[, 1]))
      if (length(keep)) {
        a <- data.table(reference_zone = as.integer(v[keep, 1]),
          availability_status = as.integer(v[keep, 2]))
        counts[[length(counts) + 1L]] <- a[, .(pixels = .N),
          by = .(reference_zone, availability_status)]
        eligible <- which(a$reference_zone %in% modeled_zones &
                          a$availability_status %in% c(1L, 2L, 3L, 4L))
        if (length(eligible)) {
          candidates <- a[eligible]
          candidates[, cell := (starts[j] - 1) * ncol(input) + keep[eligible]]
          candidates[, priority := runif(.N)]
          keys <- unique(candidates[, .(reference_zone, availability_status)])
          for (k in seq_len(nrow(keys))) {
            z <- keys$reference_zone[k]
            s <- keys$availability_status[k]
            key <- paste(z, s, sep = "_")
            cap <- if (s == 1L) 1L else n_sample
            add <- candidates[reference_zone == z & availability_status == s]
            setorder(add, priority)
            add <- head(add, cap)
            if (exists(key, envir = reservoir, inherits = FALSE)) {
              add <- rbind(get(key, envir = reservoir, inherits = FALSE), add)
              setorder(add, priority)
              add <- head(add, cap)
            }
            assign(key, add, envir = reservoir)
          }
        }
      }
      if (j %% 8L == 0L || j == length(starts)) {
        log("  Status scan: ", round(100 * j / length(starts)), "%")
      }
    }
    readStop(input)
    reading <- FALSE
    keys <- ls(reservoir, all.names = TRUE)
    if (!length(keys)) stop("No modeled-zone sample points available; report: ", out)
    summary_counts <- rbindlist(counts)[, .(pixels = sum(pixels)),
      by = .(reference_zone, availability_status)]
    summary_counts[, modeled_zone := reference_zone %in% modeled_zones]
    setorder(summary_counts, reference_zone, availability_status)
    fwrite(summary_counts, file.path(out, "01_availability_by_zone.csv"))
    points <- rbindlist(lapply(keys, get, envir = reservoir, inherits = FALSE))
    points[, priority := NULL]
    setorder(points, availability_status, reference_zone, cell)
    xy <- xyFromCell(reference, points$cell)
    points[, `:=`(sample_id = seq_len(.N), x = xy[, 1], y = xy[, 2])]
    pts <- vect(as.data.frame(points[, .(x, y)]), geom = c("x", "y"), crs = crs(reference))
    ll <- crds(project(pts, "EPSG:4326"))
    points[, `:=`(longitude = ll[, 1], latitude = ll[, 2])]
    fwrite(points, file.path(out, "02_sample_points.csv"))
    log("  Sample points: ", nrow(points), "; exact counts are in 01_availability_by_zone.csv.")

    climate_env <- Sys.getenv("ECOCHINA2_CLIMATE_DIR", "")
    climate_options <- if (nzchar(climate_env)) c(
      file.path(climate_env, "Normal_1961_1990"), climate_env) else c(
        file.path(dirname(root), "ecoChina", "play", "China", "ClimateData", "CN", "800m", "Normal_1961_1990"),
        file.path(root, "data", "rasters", "climate", "Normal_1961_1990"))
    soil_env <- Sys.getenv("ECOCHINA2_SOIL_RASTER_DIR", "")
    soil_options <- if (nzchar(soil_env)) soil_env else c(
      file.path(dirname(root), "soil rasters", "tif2"),
      file.path(root, "data", "rasters", "soil"))
    choose_directory <- function(options, kind) {
      files <- lapply(options, function(p) if (dir.exists(p)) {
        list.files(p, pattern = "[.]tif(f)?$", full.names = TRUE, ignore.case = TRUE)
      } else character())
      fwrite(data.table(kind = kind, directory = options,
        exists = dir.exists(options), tiff_files = lengths(files)),
        file.path(out, paste0(kind, "_directory_candidates.csv")))
      hit <- which(lengths(files) > 0L)
      if (length(hit)) files[[hit[1]]] else character()
    }
    climate_files <- choose_directory(climate_options, "climate")
    soil_all <- choose_directory(soil_options, "soil")
    soil_vars <- c("T_GRAVEL", "T_SAND", "T_SILT", "T_CLAY", "T_REF_BULK_DENSITY",
      "T_OC", "T_PH_H2O", "T_CEC_CLAY", "T_CEC_SOIL", "T_BS", "T_TEB",
      "T_CACO3", "T_CASO4", "T_ESP", "T_ECE")
    soil_names <- tools::file_path_sans_ext(basename(soil_all))
    soil_files <- soil_all[toupper(soil_names) %in% soil_vars]
    missing_soil_vars <- setdiff(soil_vars, toupper(soil_names))
    log("[2/3] Extract original raster values at the sampled locations.")
    log("  Climate TIFFs found: ", length(climate_files), "; topsoil TIFFs found: ", length(soil_files))
    if (!length(climate_files)) record_error("discovery", "climate",
      "No climate TIFFs found. Set ECOCHINA2_CLIMATE_DIR and rerun this check only.")
    if (length(missing_soil_vars)) record_error("discovery", "soil",
      paste("Missing training topsoil layers:", paste(missing_soil_vars, collapse = ", ")))
    values_out <- list()
    inspect <- function(path, kind) {
      r <- rast(path)
      if (nlyr(r) != 1L) stop("Expected a single-band variable TIFF; found ", nlyr(r), " bands.")
      record_geometry(r, path, kind)
      p <- if (same_crs(r, reference)) pts else project(pts, crs(r))
      xy_native <- crds(p)
      native_cell <- cellFromXY(r, xy_native)
      nearest <- extract(r, p, ID = FALSE, method = "simple")[[1]]
      linear <- extract(r, p, ID = FALSE, method = "bilinear")[[1]]
      data.table(sample_id = points$sample_id, kind = kind,
        variable = tools::file_path_sans_ext(basename(path)),
        native_cell = native_cell, inside_extent = is.finite(native_cell),
        native_value = nearest, bilinear_value = linear,
        native_finite = is.finite(nearest) & nearest != -9999,
        bilinear_finite = is.finite(linear) & linear != -9999, path = path)
    }
    all_files <- c(climate_files, soil_files)
    kinds <- c(rep("raw_climate", length(climate_files)), rep("raw_topsoil", length(soil_files)))
    for (i in seq_along(all_files)) {
      value <- tryCatch(inspect(all_files[i], kinds[i]), error = function(e) {
        record_error("point_extraction", all_files[i], conditionMessage(e)); NULL
      })
      if (!is.null(value)) values_out[[length(values_out) + 1L]] <- value
      if (i %% 10L == 0L || i == length(all_files)) log("  Input layer ", i, "/", length(all_files))
    }
    download_ref <- Sys.getenv("ECOCHINA2_CHECK_DOWNLOAD_REFERENCE", "")
    if (nzchar(download_ref)) {
      if (file.exists(download_ref)) {
        value <- tryCatch(inspect(download_ref, "optional_download_reference"), error = function(e) {
          record_error("download_reference", download_ref, conditionMessage(e)); NULL
        })
        if (!is.null(value)) values_out[[length(values_out) + 1L]] <- value
      } else record_error("download_reference", download_ref, "Configured file does not exist.")
    }
    if (length(values_out)) {
      val <- rbindlist(values_out)
      fwrite(val, file.path(out, "03_sample_raw_values.csv"))
      point_summary <- val[, .(layers_checked = .N, layers_inside_extent = sum(inside_extent),
        native_finite_layers = sum(native_finite), bilinear_finite_layers = sum(bilinear_finite),
        all_native_layers_finite = all(native_finite),
        all_bilinear_layers_finite = all(bilinear_finite),
        native_missing_variables = paste(variable[!native_finite], collapse = ";"),
        bilinear_missing_variables = paste(variable[!bilinear_finite], collapse = ";")),
        by = .(sample_id, kind)]
      point_summary <- merge(points, point_summary, by = "sample_id", all.y = TRUE, sort = FALSE)
      fwrite(point_summary, file.path(out, "04_point_input_summary.csv"))
      layer_summary <- merge(val, points[, .(sample_id, reference_zone, availability_status)],
                             by = "sample_id", all.x = TRUE, sort = FALSE)[,
        .(sampled_points = .N, native_missing_points = sum(!native_finite),
          bilinear_missing_points = sum(!bilinear_finite)),
        by = .(kind, variable, availability_status)]
      fwrite(layer_summary, file.path(out, "05_layer_missingness.csv"))
    }
    fwrite(rbindlist(geometry, fill = TRUE), file.path(out, "06_raster_geometry.csv"))
    fwrite(rbindlist(manifests, fill = TRUE), file.path(out, "07_input_files.csv"))
    err <- if (length(errors)) rbindlist(errors) else data.table(stage = character(),
      path = character(), error = character())
    fwrite(err, file.path(out, "08_errors.csv"))
    log("[3/3] Reports saved.")
    log("Exact modeled-zone missing-climate cells: ", summary_counts[
      modeled_zone & availability_status %in% c(2L, 3L), sum(pixels)])
    log("Errors or unavailable inputs: ", nrow(err), "; see 08_errors.csv.")
    log("INTERPRETATION: status 1=all candidates complete; 2=no climate candidate;")
    log("3=some climate candidates missing; 4=climate complete but soil incomplete; 0=excluded footprint.")
    log("Raw checks cover all climate TIFFs FOUND in the selected reference-period directory,")
    log("and the 15 named topsoil training variables, not every fitted model's selected-variable set.")
    log("Simple extraction reproduces native-cell lookup. Bilinear values are diagnostic only;")
    log("they are not saved predictions, filled input rasters, or evidence that a cache is wrong.")
    log("The legacy -9999 sentinel is counted as missing as in climate training; zero is retained.")
    log("Point-level findings are samples; only 01_availability_by_zone.csv has national counts.")
    log("A broad raster extent does not imply finite data. Matching grids do not imply matching masks.")
    log("Without the actual 2024 download reference, this report cannot establish its provenance.")
    log("READY TO SHARE: only the small report folder, or the ZIP below; no TIFFs or model files.")
    report_files <- list.files(out, full.names = FALSE)
    oldwd <- getwd()
    archive <- paste0(out, ".zip")
    zip_ok <- tryCatch({
      setwd(out)
      utils::zip(archive, files = report_files, flags = "-q")
      file.exists(archive) && file.info(archive)$size > 0
    }, error = function(e) FALSE, finally = setwd(oldwd))
    cat("\nReport folder: ", normalizePath(out, winslash = "/"), "\n", sep = "")
    if (zip_ok) cat("Share ZIP: ", normalizePath(archive, winslash = "/"), "\n", sep = "")
    else cat("ZIP was unavailable; manually compress the small report folder. No rerun needed.\n")
    invisible(out)
  }
  run_check()
})
