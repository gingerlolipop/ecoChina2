# EcoChina2: small-window prediction check using EXISTING fitted models.
# Put this file in H:/Jing/ecoChina2/script/ and run:
#   setwd("H:/Jing/ecoChina2")
#   source("script/check_window_predictions_local.R")
# Automatically uses the newest assessment/local_input_check_* report.
# Default: at most 20 existing sample points; climate and soil models 15 and 30.
# Models are loaded, never fitted or saved. Existing project files are read only.
# Only a NEW assessment/window_prediction_check_* folder receives output.
# Small temporary TIFFs reproduce the production FLT4S/GDAL resampling path;
# they are deleted from that new folder after each successful job.
# Optional settings, before source():
#   Sys.setenv(ECOCHINA2_WINDOW_INPUT_REPORT=".../local_input_check_20261103_130905")
#   Sys.setenv(ECOCHINA2_WINDOW_POINTS_PER_GROUP="4")   # 1..20
#   Sys.setenv(ECOCHINA2_WINDOW_MODEL_ZONES="15,30")   # 1..4 modeled zone IDs
#   Sys.setenv(ECOCHINA2_CLIMATE_DIR=".../climate")
#   Sys.setenv(ECOCHINA2_SOIL_RASTER_DIR=".../soil TIFF folder")

local({
  # Pure helpers are also exposed through TEST_MODE for isolated verification.
  float32 <- function(x) {
    readBin(writeBin(as.numeric(x), raw(), size = 4L), what = "numeric",
            n = length(x), size = 4L)
  }
  same_value <- function(a, b, tolerance = 1e-5) {
    if (!is.finite(a) && !is.finite(b)) return(TRUE)
    is.finite(a) && is.finite(b) && abs(a - b) <= tolerance
  }
  classify <- function(stored, small, large, error = "", tolerance = 1e-5) {
    if (nzchar(error)) return("job_error")
    if (!same_value(small, large, tolerance)) return("window_context_not_stable")
    if (!is.finite(stored) && !is.finite(large)) return("stored_and_recomputed_missing")
    if (!is.finite(stored) && is.finite(large)) return("stored_missing_local_pipeline_finite")
    if (is.finite(stored) && !is.finite(large)) return("stored_finite_local_pipeline_missing")
    if (abs(stored - large) > tolerance) return("stored_and_recomputed_values_differ")
    "stored_and_recomputed_agree"
  }
  select_spread <- function(d, cap) {
    if (nrow(d) <= cap) return(d)
    # Even spacing through zone/latitude/longitude order; diagnostic selection,
    # not a probability sample and not a national area/frequency estimator.
    d[unique(as.integer(round(seq(1, nrow(d), length.out = cap)))), ]
  }

  main <- function() {
    for (pkg in c("terra", "data.table", "randomForest")) {
      if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing installed package: ", pkg)
    }
    library(terra)
    library(data.table)
    # Loading the package registers randomForest's predict method.
    loadNamespace("randomForest")
    root <- normalizePath(Sys.getenv("ECOCHINA2_DIR", unset = getwd()),
                          winslash = "/", mustWork = TRUE)
    repeat {
      if (file.exists(file.path(root, "script", "1. data.R"))) break
      up <- dirname(root)
      if (identical(up, root)) stop("Run from ecoChina2 or set ECOCHINA2_DIR.")
      root <- up
    }
    modeled <- c(1:7, 9:50, 52:55)
    model_zones <- suppressWarnings(as.integer(strsplit(
      Sys.getenv("ECOCHINA2_WINDOW_MODEL_ZONES", "15,30"), ",", fixed = TRUE)[[1]]))
    if (!length(model_zones) || length(model_zones) > 4L || anyNA(model_zones) ||
        anyDuplicated(model_zones) || any(!model_zones %in% modeled)) {
      stop("WINDOW_MODEL_ZONES must contain 1..4 distinct modeled zone IDs.")
    }
    cap <- suppressWarnings(as.integer(Sys.getenv("ECOCHINA2_WINDOW_POINTS_PER_GROUP", "4")))
    if (length(cap) != 1L || is.na(cap) || cap < 1L || cap > 20L) {
      stop("WINDOW_POINTS_PER_GROUP must be an integer from 1 through 20.")
    }
    report <- Sys.getenv("ECOCHINA2_WINDOW_INPUT_REPORT", "")
    if (!nzchar(report)) {
      candidates <- list.dirs(file.path(root, "assessment"), recursive = FALSE, full.names = TRUE)
      candidates <- candidates[grepl("^local_input_check_", basename(candidates))]
      candidates <- candidates[vapply(candidates, function(p) all(file.exists(file.path(p,
        c("02_sample_points.csv", "04_point_input_summary.csv", "07_input_files.csv")))), logical(1))]
      if (!length(candidates)) stop("Existing local_input_check_* report not found.\n",
        "Set ECOCHINA2_WINDOW_INPUT_REPORT to that report folder; do not rerun training.")
      report <- tail(sort(candidates), 1L)
    }
    report <- normalizePath(report, winslash = "/", mustWork = TRUE)
    reference_file <- file.path(root, "raster", "ecosys_ori.tif")
    if (!file.exists(reference_file)) stop("Missing vegetation reference: ", reference_file)
    reference <- rast(reference_file)
    if (nlyr(reference) != 1L || !nzchar(crs(reference))) stop("Invalid vegetation reference.")
    points <- fread(file.path(report, "02_sample_points.csv"))
    summary <- fread(file.path(report, "04_point_input_summary.csv"))
    prior <- fread(file.path(report, "07_input_files.csv"))
    needed_points <- c("sample_id", "cell", "reference_zone", "availability_status", "x", "y")
    needed_summary <- c("sample_id", "kind", "layers_checked", "native_finite_layers",
                        "bilinear_finite_layers", "availability_status")
    if (!all(needed_points %in% names(points)) || !all(needed_summary %in% names(summary)) ||
        !all(c("kind", "path", "size_bytes", "modified") %in% names(prior))) {
      stop("The selected earlier report has an unexpected schema.")
    }
    if (anyDuplicated(points$sample_id) || anyDuplicated(summary[, .(sample_id, kind)])) {
      stop("Duplicate sample IDs or duplicate sample/kind rows in the earlier report.")
    }
    points <- points[reference_zone %in% modeled]
    if (!nrow(points)) stop("No modeled-zone points in the earlier report.")
    xy_expected <- xyFromCell(reference, points$cell)
    xy_scale <- max(res(reference)) * 1e-6
    if (any(!is.finite(points$x) | !is.finite(points$y)) ||
        any(abs(points$x - xy_expected[, 1]) > xy_scale) ||
        any(abs(points$y - xy_expected[, 2]) > xy_scale)) {
      stop("Sample cell IDs/coordinates do not match the current reference.")
    }
    observed <- extract(reference, as.matrix(points[, .(x, y)]))[[1]]
    if (any(!is.finite(observed)) || any(observed != points$reference_zone)) {
      stop("Sample vegetation labels changed since the earlier report.")
    }
    cl <- summary[kind == "raw_climate"]
    so <- summary[kind == "raw_topsoil"]
    pools <- list(
      valid_control = points[availability_status == 1L],
      climate_missing_bilinear_support = points[sample_id %in% cl[
        availability_status == 2L & native_finite_layers == 0L & layers_checked > 0L &
        bilinear_finite_layers == layers_checked, sample_id]],
      climate_missing_no_bilinear_support = points[sample_id %in% cl[
        availability_status == 2L & native_finite_layers == 0L &
        bilinear_finite_layers == 0L, sample_id]],
      soil_missing_bilinear_support = points[sample_id %in% so[
        availability_status == 4L & native_finite_layers == 0L & layers_checked > 0L &
        bilinear_finite_layers == layers_checked, sample_id]],
      soil_missing_no_bilinear_support = points[sample_id %in% so[
        availability_status == 4L & native_finite_layers == 0L &
        bilinear_finite_layers == 0L, sample_id]])
    group_counts <- rbindlist(lapply(names(pools), function(group) data.table(
      group = group, available_sample_points = nrow(pools[[group]]),
      selected_points = min(cap, nrow(pools[[group]])))))
    if (!nrow(pools$valid_control)) stop("At least one valid control is required.")
    selected <- rbindlist(lapply(names(pools), function(group) {
      d <- copy(pools[[group]])
      setorder(d, reference_zone, y, x)
      d <- select_spread(d, cap)
      d[, group := rep(group, .N)]
      d
    }), fill = TRUE)
    if (anyDuplicated(selected$sample_id)) stop("A point belongs to more than one diagnostic group.")
    out_base <- file.path(root, "assessment", paste0("window_prediction_check_",
                                                   format(Sys.time(), "%Y%m%d_%H%M%S")))
    out <- out_base
    suffix <- 0L
    while (dir.exists(out)) { suffix <- suffix + 1L; out <- paste0(out_base, "_", suffix) }
    dir.create(out, recursive = TRUE)
    temporary <- file.path(out, "temporary")
    dir.create(temporary)
    log <- function(...) {
      line <- paste0(...)
      cat(line, "\n")
      cat(line, "\n", file = file.path(out, "00_summary.txt"), append = TRUE)
      flush.console()
    }
    log("SMALL-WINDOW PREDICTION CHECK: reference period; main mf_var workflow")
    log("Earlier report: ", report)
    log("Existing fitted models are loaded only. No training, filling or project-output replacement.")
    log("Selected points: ", nrow(selected), "; fitted-model zone IDs: ", paste(model_zones, collapse = ","))
    log("Each point uses native radius 4 / target radius 1 and native radius 8 / target radius 3.")
    fwrite(group_counts, file.path(out, "01_selection_groups.csv"))
    fwrite(selected, file.path(out, "02_selected_points.csv"))
    writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
    errors <- list()
    metadata <- list()
    signatures <- list()
    variable_records <- list()
    rows <- list()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) saved_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
      if (had_seed) assign(".Random.seed", saved_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(list = ".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    add_error <- function(stage, kind, zone, sample_id, msg) {
      errors[[length(errors) + 1L]] <<- data.table(stage = stage, kind = kind,
        model_zone = zone, sample_id = sample_id, error = msg)
    }
    register_files <- function(paths, kind) {
      paths <- unique(paths)
      z <- file.info(paths)
      signatures[[length(signatures) + 1L]] <<- data.table(kind = kind, path = paths,
        exists_before = file.exists(paths), size_before = z$size,
        mtime_before = as.numeric(z$mtime), modified_before = as.character(z$mtime))
    }
    register_files(c(reference_file, file.path(report,
      c("02_sample_points.csv", "04_point_input_summary.csv", "07_input_files.csv"))), "reference_or_prior_report")
    same_crs <- function(a, b) compareGeom(a[[1]], b[[1]], crs = TRUE,
      ext = FALSE, rowcol = FALSE, res = FALSE, stopOnError = FALSE)
    # Geometry from exact source row/column edges; never reconstruct from rounded lon/lat.
    local_window <- function(r, cell, radius) {
      rr <- rowFromCell(r, cell); cc <- colFromCell(r, cell)
      r1 <- max(1L, rr - radius); r2 <- min(nrow(r), rr + radius)
      c1 <- max(1L, cc - radius); c2 <- min(ncol(r), cc + radius)
      box <- ext(xmin(r) + (c1 - 1) * res(r)[1], xmin(r) + c2 * res(r)[1],
                 ymax(r) - r2 * res(r)[2], ymax(r) - (r1 - 1) * res(r)[2])
      crop(r, box, snap = "near")
    }
    get_varlist <- function(model) {
      vars <- model$varlist
      if (is.null(vars) || !length(vars)) vars <- rownames(model$importance)
      vars <- as.character(vars)
      if (!length(vars) || anyNA(vars) || any(!nzchar(vars)) || anyDuplicated(vars)) {
        stop("The fitted model has no valid unique predictor list.")
      }
      vars
    }
    raw_paths <- function(kind, vars) {
      env <- Sys.getenv(if (kind == "climate") "ECOCHINA2_CLIMATE_DIR" else "ECOCHINA2_SOIL_RASTER_DIR", "")
      prior_kind <- if (kind == "climate") "raw_climate" else "raw_topsoil"
      known <- unique(prior[kind == prior_kind, path])
      directories <- if (nzchar(env)) {
        if (kind == "climate") c(file.path(env, "Normal_1961_1990"), env) else env
      } else unique(dirname(known))
      for (d in directories) {
        if (!dir.exists(d)) next
        files <- list.files(d, "[.]tif(f)?$", full.names = TRUE, ignore.case = TRUE)
        stems <- tools::file_path_sans_ext(basename(files))
        if (anyDuplicated(tolower(stems))) next
        pos <- match(tolower(vars), tolower(stems))
        if (all(!is.na(pos))) return(files[pos])
      }
      stop("No single ", kind, " directory contains all fitted variables. Checked: ",
           paste(directories, collapse = "; "))
    }
    # Write only explicitly generated small rasters within this new temporary folder.
    temporary_path <- function(prefix) {
      f <- tempfile(prefix, tmpdir = temporary, fileext = ".tif")
      if (!identical(normalizePath(dirname(f), winslash = "/", mustWork = TRUE),
                     normalizePath(temporary, winslash = "/", mustWork = TRUE))) stop("Invalid temporary path.")
      f
    }
    temporary_write <- function(r, prefix) {
      f <- temporary_path(prefix)
      writeRaster(r, f, overwrite = FALSE, wopt = list(datatype = "FLT4S", gdal = "COMPRESS=NONE"))
      f
    }
    one_window <- function(raw, model, vars, point, native_radius, target_radius, raw_same_grid) {
      xy <- matrix(c(point$x, point$y), nrow = 1L)
      native_cell <- cellFromXY(raw, xy)
      if (!is.finite(native_cell)) stop("Point is outside original predictor extent.")
      win <- local_window(raw, native_cell, native_radius)
      w <- values(win, mat = TRUE)
      if (any(w == -9999, na.rm = TRUE)) stop("Undeclared -9999 sentinel in raw window; no automatic reinterpretation attempted.")
      if (any(is.infinite(w))) stop("Infinite raw predictor value found.")
      complete <- rowSums(is.finite(w)) == length(vars)
      prediction <- rep(NA_real_, nrow(w))
      if (any(complete)) {
        newdata <- as.data.frame(w[complete, , drop = FALSE]); names(newdata) <- vars
        pp <- stats::predict(model, newdata, type = "prob")
        if (!is.matrix(pp) || !"1" %in% colnames(pp)) stop("Fitted model has no class-1 probability.")
        prediction[complete] <- pp[, "1"]
        if (any(!is.finite(prediction[complete])) || any(prediction[complete] < 0 | prediction[complete] > 1)) {
          stop("Invalid fitted-model probability.")
        }
      }
      local_prediction <- rast(win[[1]])
      values(local_prediction) <- prediction
      native_file <- temporary_write(local_prediction, "native_")
      cleanup <- native_file
      on.exit({
        for (f in cleanup) {
          if (identical(dirname(f), temporary)) unlink(c(f, paste0(f, ".aux.xml")), force = TRUE)
        }
      }, add = TRUE)
      native_prediction <- rast(native_file)
      p_native <- as.numeric(extract(native_prediction, xy, method = "simple")[[1]])
      p_extract_bilinear <- as.numeric(extract(native_prediction, xy, method = "bilinear")[[1]])
      if (raw_same_grid) {
        p_resampled <- p_native
      } else {
        target <- local_window(reference, point$cell, target_radius)
        aligned_file <- temporary_path("aligned_")
        cleanup <- c(cleanup, aligned_file)
        aligned <- resample(native_prediction, target, method = "bilinear",
          filename = aligned_file, overwrite = FALSE,
          wopt = list(datatype = "FLT4S", gdal = "COMPRESS=NONE"))
        p_resampled <- as.numeric(extract(rast(aligned_file), xy, method = "simple")[[1]])
      }
      raw_native <- as.numeric(unlist(extract(win, xy, method = "simple")[1, ], use.names = FALSE))
      raw_linear <- as.numeric(unlist(extract(win, xy, method = "bilinear")[1, ], use.names = FALSE))
      answer <- list(prediction = p_resampled, native_prediction = p_native,
        probability_extract_bilinear = p_extract_bilinear,
        native_complete_cells = sum(complete), native_window_cells = nrow(w),
        raw_native_finite_variables = sum(is.finite(raw_native)),
        raw_bilinear_finite_variables = sum(is.finite(raw_linear)))
      # Release file-backed raster handles before on.exit cleanup on Windows.
      rm(native_prediction, local_prediction, win)
      if (exists("aligned")) rm(aligned)
      answer
    }
    for (kind_name in c("climate", "soil")) {
      for (zone in model_zones) {
        log("[MODEL] ", kind_name, " zone ", zone)
        model_file <- file.path(root, if (kind_name == "climate") "rf" else "rf_soil",
          paste0(if (kind_name == "climate") "clm_mfVar_zone" else "soil_mf_zone", zone, ".Rdata"))
        object_name <- if (kind_name == "climate") "clm_mfVar" else "soil_mf"
        saved_file <- file.path(root, if (kind_name == "climate") "clim suitability" else "soil suitability",
          if (kind_name == "climate") "mf_var" else "plain_mf", "normal",
          paste0(if (kind_name == "climate") "clim_suit_zone" else "soil_suit_zone", zone, ".tif"))
        model_error <- ""
        prepared <- tryCatch({
          if (!all(file.exists(c(model_file, saved_file)))) stop("Missing existing model or suitability file: ",
            paste(c(model_file, saved_file)[!file.exists(c(model_file, saved_file))], collapse = "; "))
          env <- new.env(parent = baseenv())
          load(model_file, envir = env)
          if (!exists(object_name, envir = env, inherits = FALSE)) stop("Expected object absent: ", object_name)
          model <- get(object_name, envir = env, inherits = FALSE)
          vars <- get_varlist(model)
          files <- raw_paths(kind_name, vars)
          raw <- rast(files); names(raw) <- vars
          saved <- rast(saved_file)
          if (nlyr(saved) != 1L || !compareGeom(saved, reference, stopOnError = FALSE)) stop("Saved suitability geometry differs from vegetation.")
          if (!same_crs(raw, reference)) stop("Predictor CRS differs; production uses resample, so no automatic CRS conversion attempted.")
          raw_same_grid <- compareGeom(raw[[1]], reference, stopOnError = FALSE)
          register_files(c(model_file, saved_file, files), paste0(kind_name, "_model_", zone))
          metadata[[length(metadata) + 1L]] <- data.table(kind = kind_name, model_zone = zone,
            model_file = model_file, saved_suitability = saved_file,
            model_class = paste(class(model), collapse = ";"), n_variables = length(vars),
            variables = paste(vars, collapse = ";"), raw_same_grid_as_vegetation = raw_same_grid)
          variable_records[[length(variable_records) + 1L]] <- data.table(kind = kind_name,
            model_zone = zone, variable = vars, raw_file = files)
          list(model = model, vars = vars, raw = raw, saved = saved, raw_same_grid = raw_same_grid)
        }, error = function(e) { model_error <<- conditionMessage(e); NULL })
        if (is.null(prepared)) {
          add_error("model_preparation", kind_name, zone, NA_integer_, model_error)
          log("  [ERROR] ", model_error)
          next
        }
        for (j in seq_len(nrow(selected))) {
          point <- as.list(selected[j])
          result <- data.table(sample_id = point$sample_id, group = point$group,
            reference_zone = point$reference_zone, availability_status = point$availability_status,
            cell = point$cell, x = point$x, y = point$y, kind = kind_name, model_zone = zone,
            n_variables = length(prepared$vars), stored_suitability = NA_real_,
            local_small = NA_real_, local_large = NA_real_, local_native = NA_real_,
            local_probability_extract_bilinear = NA_real_,
            native_complete_cells_large = NA_integer_, native_window_cells_large = NA_integer_,
            raw_native_finite_variables = NA_integer_, raw_bilinear_finite_variables = NA_integer_,
            error = "")
          result <- tryCatch({
            xy <- matrix(c(point$x, point$y), nrow = 1L)
            result[, stored_suitability := as.numeric(extract(prepared$saved, xy)[[1]])]
            small <- one_window(prepared$raw, prepared$model, prepared$vars, point, 4L, 1L, prepared$raw_same_grid)
            large <- one_window(prepared$raw, prepared$model, prepared$vars, point, 8L, 3L, prepared$raw_same_grid)
            result[, `:=`(local_small = small$prediction, local_large = large$prediction,
              local_native = large$native_prediction,
              local_probability_extract_bilinear = large$probability_extract_bilinear,
              native_complete_cells_large = large$native_complete_cells,
              native_window_cells_large = large$native_window_cells,
              raw_native_finite_variables = large$raw_native_finite_variables,
              raw_bilinear_finite_variables = large$raw_bilinear_finite_variables)]
            result
          }, error = function(e) {
            msg <- conditionMessage(e)
            result[, error := msg]
            add_error("point_window", kind_name, zone, point$sample_id, msg)
            result
          })
          result[, diagnosis := classify(stored_suitability, local_small, local_large, error)]
          result[, window_context_stable := !nzchar(error) && same_value(local_small, local_large)]
          result[, absolute_difference := if (is.finite(stored_suitability) && is.finite(local_large)) {
            abs(stored_suitability - local_large)
          } else NA_real_]
          rows[[length(rows) + 1L]] <- result
          if (j %% 4L == 0L || j == nrow(selected)) log("  Checked points ", j, "/", nrow(selected))
        }
        # Save current progress before moving to the next model.
        fwrite(rbindlist(rows, fill = TRUE), file.path(out, "03_window_comparison.csv"))
        rm(prepared)
        # 'model', 'raw' and 'saved' created in tryCatch are local to main; release them.
        for (nm in intersect(c("model", "raw", "saved", "env"), ls())) rm(list = nm)
        gc()
      }
    }
    comparison <- if (length(rows)) rbindlist(rows, fill = TRUE) else data.table()
    expected_models <- CJ(kind = c("climate", "soil"), model_zone = model_zones)
    controls <- if (nrow(comparison)) comparison[group == "valid_control", .(
      control_points_checked = .N,
      control_points_passed = sum(diagnosis == "stored_and_recomputed_agree"),
      max_control_absolute_difference = if (any(is.finite(absolute_difference))) max(absolute_difference, na.rm = TRUE) else NA_real_),
      by = .(kind, model_zone)] else data.table(kind = character(), model_zone = integer(),
        control_points_checked = integer(), control_points_passed = integer(), max_control_absolute_difference = numeric())
    controls <- merge(expected_models, controls, by = c("kind", "model_zone"), all.x = TRUE)
    n_control <- nrow(selected[group == "valid_control"])
    controls[, controls_pass := !is.na(control_points_checked) & control_points_checked == n_control &
               control_points_passed == n_control]
    controls[is.na(controls_pass), controls_pass := FALSE]
    fwrite(controls, file.path(out, "04_control_validation.csv"))
    if (nrow(comparison)) {
      comparison <- merge(comparison, controls[, .(kind, model_zone, controls_pass)],
                          by = c("kind", "model_zone"), all.x = TRUE, sort = FALSE)
      fwrite(comparison, file.path(out, "03_window_comparison.csv"))
      diagnosis_summary <- comparison[, .(checked_point_model_pairs = .N),
        by = .(kind, model_zone, group, controls_pass, diagnosis)]
      fwrite(diagnosis_summary, file.path(out, "05_diagnosis_counts.csv"))
    }
    if (length(metadata)) fwrite(rbindlist(metadata), file.path(out, "06_models.csv"))
    if (length(variable_records)) fwrite(rbindlist(variable_records), file.path(out, "07_fitted_variables.csv"))
    signature <- unique(rbindlist(signatures), by = c("kind", "path"))
    end_info <- file.info(signature$path)
    signature[, `:=`(exists_after = file.exists(path), size_after = end_info$size,
                     mtime_after = as.numeric(end_info$mtime))]
    signature[, unchanged_during_check := exists_before & exists_after &
      size_before == size_after & mtime_before == mtime_after]
    fwrite(signature, file.path(out, "08_input_file_checks.csv"))
    prior_check <- merge(prior, signature[, .SD[1], by = path][,
      .(path, size_before, modified_before)], by = "path", all.x = TRUE, sort = FALSE)
    prior_check[, `:=`(used_in_window_check = !is.na(size_before),
      size_matches_previous_report = size_bytes == size_before,
      displayed_mtime_matches_previous_report = modified == modified_before)]
    fwrite(prior_check, file.path(out, "09_previous_input_checks.csv"))
    error_table <- if (length(errors)) rbindlist(errors) else data.table(stage = character(),
      kind = character(), model_zone = integer(), sample_id = integer(), error = character())
    fwrite(error_table, file.path(out, "10_errors.csv"))
    log("CHECK FINISHED. Successful control model sets: ", sum(controls$controls_pass), "/", nrow(controls))
    log("Errors: ", nrow(error_table), "; input files changed during check: ",
        sum(!signature$unchanged_during_check | is.na(signature$unchanged_during_check)))
    log("INTERPRETATION:")
    log("A model's missing-point comparison is interpretable only alongside its control validation,")
    log("stable window contexts, input-version records, and any errors.")
    log("stored_and_recomputed_missing: the current native-prediction/resampling route retains NA.")
    log("stored_missing_local_pipeline_finite: inspect source/model versions and global-vs-local warp behavior;")
    log("this is a diagnostic flag, not authorization or proof that an old national cell should be replaced.")
    log("Bilinear extraction of raw variables or probabilities can differ from raster resampling.")
    log("No prediction on interpolated predictor values is substituted for the original workflow.")
    log("Only selected points and models were checked; no national recovery fraction or retraining decision is proved.")
    log("Share the small CSV/TXT report ZIP. Any remaining small temporary TIFFs are excluded from it.")
    archive <- paste0(out, ".zip")
    original_wd <- getwd()
    zip_ok <- tryCatch({
      setwd(out)
      report_files <- list.files(out, pattern = "[.](csv|txt)$", full.names = FALSE, ignore.case = TRUE)
      utils::zip(archive, files = report_files, flags = "-q")
      file.exists(archive) && file.info(archive)$size > 0
    }, error = function(e) FALSE, finally = setwd(original_wd))
    cat("\nReport folder: ", normalizePath(out, winslash = "/"), "\n", sep = "")
    if (zip_ok) cat("Share ZIP: ", normalizePath(archive, winslash = "/"), "\n", sep = "")
    else cat("ZIP unavailable; manually compress only this report folder's CSV/TXT files. No rerun needed.\n")
    invisible(out)
  }
  if (identical(Sys.getenv("ECOCHINA2_WINDOW_TEST_MODE"), "true")) {
    list(float32 = float32, same_value = same_value, classify = classify,
         select_spread = select_spread, main = main)
  } else main()
})
