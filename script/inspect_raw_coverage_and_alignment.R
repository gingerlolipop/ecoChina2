# Read-only follow-up to trace_climate_na_points.R. No models are loaded/fitted.
# No project raster is changed, filled, resampled nationally, or overwritten.
# Copy this ONE file to script/ and source it from the ecoChina2 project root.
# Requires terra and data.table; uses the previous trace's saved sample points.
local({
  run_check <- function() {
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
    assessment <- file.path(root, "assessment")
    trace_dir <- Sys.getenv("ECOCHINA2_COVERAGE_TRACE_DIR", unset = "")
    if (!nzchar(trace_dir)) {
      candidates <- list.dirs(assessment, recursive = FALSE, full.names = TRUE)
      candidates <- candidates[grepl("^climate_na_trace_", basename(candidates)) &
        file.exists(file.path(candidates, "sample_points.csv"))]
      if (!length(candidates)) stop("No previous climate_na_trace_*/sample_points.csv found.\n",
        "Set ECOCHINA2_COVERAGE_TRACE_DIR to the folder from the previous diagnostic.")
      trace_dir <- tail(sort(candidates), 1L)
    }
    trace_dir <- normalizePath(trace_dir, winslash = "/", mustWork = TRUE)
    radius <- as.integer(Sys.getenv("ECOCHINA2_COVERAGE_WINDOW_RADIUS", "16"))
    if (is.na(radius) || radius < 3L || radius > 64L) {
      stop("ECOCHINA2_COVERAGE_WINDOW_RADIUS must be between 3 and 64 cells.")
    }
    out <- Sys.getenv("ECOCHINA2_COVERAGE_OUT", unset = file.path(assessment,
      paste0("raw_coverage_alignment_", format(Sys.time(), "%Y%m%d_%H%M%S"))))
    if (file.exists(out)) stop("Output already exists; choose a NEW ECOCHINA2_COVERAGE_OUT: ", out)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    temp <- file.path(out, "temporary"); dir.create(temp)
    old_options <- terraOptions(print = FALSE)
    terraOptions(tempdir = temp, memfrac = 0.15)
    on.exit({
      old_temp <- old_options$tempdir
      if (is.null(old_temp) || !dir.exists(old_temp)) old_temp <- tempdir()
      terraOptions(tempdir = old_temp, memfrac = old_options$memfrac)
      unlink(temp, recursive = TRUE)
    }, add = TRUE)
    reference_file <- file.path(root, "raster", "ecosys_ori.tif")
    reference <- rast(reference_file)
    if (!nzchar(crs(reference))) stop("Reference CRS is missing; no CRS is assigned by this script.")
    points_file <- file.path(trace_dir, "sample_points.csv")
    points <- fread(points_file)
    if (!all(c("sample_id", "x", "y", "point_type") %in% names(points)) ||
        anyDuplicated(points$sample_id) || !nrow(points) || nrow(points) > 128L) {
      stop("Invalid saved sample_points.csv (expect 1-128 uniquely numbered points).")
    }
    if (any(!is.finite(points$x) | !is.finite(points$y))) stop("Non-finite sample coordinates.")
    point_vector <- vect(as.data.frame(points[, .(x, y)]), geom = c("x", "y"), crs = crs(reference))
    ll <- crds(project(point_vector, "EPSG:4326"))
    points[, `:=`(longitude = ll[, 1], latitude = ll[, 2])]
    points[, reference_zone := extract(reference, point_vector, ID = FALSE)[[1]]]
    points[, cell := cellFromXY(reference, as.matrix(points[, .(x, y)]))]
    if (any(!is.finite(points$cell))) stop("Saved points fall outside the current reference raster.")
    fwrite(points, file.path(out, "sample_points.csv"))
    variable_files <- list.files(trace_dir, "^model_zone[0-9]+_variables[.]csv$", full.names = TRUE)
    if (!length(variable_files)) stop("Missing saved model_zone*_variables.csv in ", trace_dir)
    model_vars <- lapply(variable_files, function(file) {
      a <- fread(file)
      if (!"variable" %in% names(a) || anyNA(a$variable) || !nrow(a) ||
          any(!nzchar(a$variable)) || anyDuplicated(a$variable)) stop("Invalid variable list: ", file)
      as.character(a$variable)
    })
    names(model_vars) <- sub("^model_zone([0-9]+)_variables[.]csv$", "\\1", basename(variable_files))
    union_vars <- sort(unique(unlist(model_vars, use.names = FALSE)))
    fwrite(rbindlist(lapply(names(model_vars), function(z) data.table(
      model_zone = as.integer(z), variable = model_vars[[z]]))), file.path(out, "climate_variables_by_model.csv"))

    # Inspect every valid source copy rather than silently asserting cache provenance.
    env_climate <- Sys.getenv("ECOCHINA2_CLIMATE_DIR", unset = "")
    climate_candidates <- if (nzchar(env_climate)) {
      c(env_climate, file.path(env_climate, "Normal_1961_1990"))
    } else {
      previous_candidates <- file.path(trace_dir, "raw_directory_candidates.csv")
      old_dirs <- if (file.exists(previous_candidates)) fread(previous_candidates)$directory else character()
      previous_points <- file.path(trace_dir, "point_pipeline_check.csv")
      checked_dirs <- if (file.exists(previous_points)) fread(previous_points)$raw_directory else character()
      c(checked_dirs, old_dirs, file.path(dirname(root), "ecoChina", "play", "China",
        "ClimateData", "CN", "800m", "Normal_1961_1990"),
        file.path(root, "data", "rasters", "climate", "Normal_1961_1990"))
    }
    climate_candidates <- unique(climate_candidates[!is.na(climate_candidates) & nzchar(climate_candidates)])
    full <- vapply(climate_candidates, function(p) all(file.exists(file.path(p,
      paste0(union_vars, ".tif")))), logical(1))
    fwrite(data.table(directory = climate_candidates, exists = dir.exists(climate_candidates),
      contains_all_inspected_climate_variables = full), file.path(out, "climate_directory_candidates.csv"))
    if (!any(full)) stop("No climate directory contains the saved variable lists. Set ECOCHINA2_CLIMATE_DIR.")
    climate_dirs <- unique(normalizePath(climate_candidates[full], winslash = "/", mustWork = TRUE))
    soil_env <- Sys.getenv("ECOCHINA2_SOIL_RASTER_DIR", unset = "")
    soil_candidates <- if (nzchar(soil_env)) soil_env else c("H:/Jing/soil rasters/tif2",
      file.path(root, "data", "rasters", "soil"))
    soil_candidates <- unique(soil_candidates)
    soil_listing <- lapply(soil_candidates, function(p) list.files(p, "[.]tif$", full.names = TRUE, ignore.case = TRUE))
    usable <- lengths(soil_listing) > 0L
    fwrite(data.table(directory = soil_candidates, exists = dir.exists(soil_candidates),
      tif_files = lengths(soil_listing), chosen = usable & seq_along(usable) == match(TRUE, usable)),
      file.path(out, "soil_directory_candidates.csv"))
    if (!any(usable)) stop("No soil TIFFs found. Set ECOCHINA2_SOIL_RASTER_DIR.")
    soil_files <- sort(soil_listing[[which(usable)[1]]])
    climate_files <- lapply(climate_dirs, function(p) file.path(p, paste0(union_vars, ".tif")))
    all_inputs <- c(reference_file, points_file, variable_files, soil_files,
                    unlist(climate_files, use.names = FALSE))
    append_csv <- function(x, name) {
      destination <- file.path(out, name)
      fwrite(x, destination, append = file.exists(destination), col.names = !file.exists(destination))
    }
    same_crs <- function(a, b) isTRUE(compareGeom(a, b, ext = FALSE, rowcol = FALSE,
                                                res = FALSE, stopOnError = FALSE))
    same_lattice <- function(a, b) {
      if (!same_crs(a, b) || any(abs(res(a) - res(b)) > 1e-10)) return(FALSE)
      d <- (origin(a) - origin(b)) / res(b)
      all(abs(d - round(d)) < 1e-6)
    }
    soil_template <- rast(soil_files[1])
    climate_template <- rast(climate_files[[1]][1])
    layer_rows <- list(); point_rows <- list()
    record_layer <- function(file, kind, variable, source_id) {
      info <- file.info(file)
      common <- list(kind = kind, source_id = source_id, variable = variable,
        file = normalizePath(file, winslash = "/", mustWork = TRUE),
        size_bytes = as.numeric(info$size), modified = as.character(info$mtime))
      error <- ""
      z <- tryCatch({
        x <- rast(file)
        if (nlyr(x) != 1L) stop("Expected a single-layer input TIFF.")
        geometry <- cbind(as.data.table(common), data.table(rows = nrow(x), columns = ncol(x),
          xmin = xmin(x), xmax = xmax(x), ymin = ymin(x), ymax = ymax(x),
          x_resolution = res(x)[1], y_resolution = res(x)[2],
          x_origin = origin(x)[1], y_origin = origin(x)[2],
          crs_present = nzchar(crs(x)), crs = gsub("[\r\n]+", " ", crs(x)),
          same_crs_as_reference = same_crs(x, reference),
          same_lattice_as_reference = same_lattice(x, reference),
          same_geometry_as_reference = isTRUE(compareGeom(x, reference, stopOnError = FALSE)),
          same_crs_as_first_soil = same_crs(x, soil_template),
          same_lattice_as_first_soil = same_lattice(x, soil_template),
          same_geometry_as_first_soil = isTRUE(compareGeom(x, soil_template, stopOnError = FALSE)),
          same_crs_as_first_climate = same_crs(x, climate_template),
          same_lattice_as_first_climate = same_lattice(x, climate_template),
          same_geometry_as_first_climate = isTRUE(compareGeom(x, climate_template, stopOnError = FALSE)), error = ""))
        pp <- if (nzchar(crs(x))) project(point_vector, crs(x)) else NULL
        if (is.null(pp)) stop("Input CRS is missing; not guessed.")
        xy <- crds(pp); cells <- cellFromXY(x, xy)
        value <- extract(x, pp, ID = FALSE)[[1]]
        ptable <- data.table(kind = kind, source_id = source_id, variable = variable,
          file = common$file, sample_id = points$sample_id, point_type = points$point_type,
          source_x = xy[, 1], source_y = xy[, 2], source_cell = cells,
          inside_source_extent = is.finite(cells), value = value,
          finite = is.finite(value), error = "")
        list(geometry = geometry, points = ptable)
      }, error = function(e) { error <<- conditionMessage(e); NULL })
      if (is.null(z)) {
        z <- list(geometry = cbind(as.data.table(common), data.table(error = error)), points = data.table(
          kind = kind, source_id = source_id, variable = variable, file = common$file,
          sample_id = points$sample_id, point_type = points$point_type, source_x = NA_real_,
          source_y = NA_real_, source_cell = NA_real_, inside_source_extent = NA,
          value = NA_real_, finite = FALSE, error = error))
      }
      layer_rows[[length(layer_rows) + 1L]] <<- z$geometry
      point_rows[[length(point_rows) + 1L]] <<- z$points
      invisible(NULL)
    }
    cat("[CHECK 1/3] Raw soil headers and values at saved points; no national scan.\n")
    for (i in seq_along(soil_files)) {
      record_layer(soil_files[i], "raw_soil", tools::file_path_sans_ext(basename(soil_files[i])), 1L)
      cat(sprintf("  Soil layer %d/%d: %s\n", i, length(soil_files), basename(soil_files[i])))
      flush.console()
    }
    for (source_id in seq_along(climate_files)) {
      for (i in seq_along(union_vars)) record_layer(climate_files[[source_id]][i],
        "raw_climate", union_vars[i], source_id)
    }
    headers <- rbindlist(layer_rows, fill = TRUE)
    pvalues <- rbindlist(point_rows, fill = TRUE)
    fwrite(headers, file.path(out, "raw_layer_geometry.csv"))
    fwrite(headers[kind == "raw_soil"], file.path(out, "soil_layer_geometry.csv"))
    fwrite(pvalues, file.path(out, "raw_layer_point_values.csv"))
    soil_point <- pvalues[kind == "raw_soil", .(soil_layers_inspected = .N,
      soil_layers_read_errors = sum(nzchar(error)), soil_layers_finite = sum(finite)), by = sample_id]
    soil_point <- soil_point[match(points$sample_id, sample_id)]
    # Soil counts are descriptive across listed layers, NOT model applicability.
    fwrite(soil_point, file.path(out, "soil_point_coverage_descriptive.csv"))
    saved_root <- Sys.getenv("ECOCHINA2_AUDIT_CLIMATE_SUITABILITY_DIR", unset = file.path(root, "clim suitability", "mf_var"))
    soil_saved_root <- Sys.getenv("ECOCHINA2_AUDIT_SOIL_SUITABILITY_DIR", unset = file.path(root, "soil suitability", "plain_mf", "normal"))
    get_saved <- function(file) {
      if (!file.exists(file)) return(list(values = rep(NA_real_, nrow(points)), error = "saved_file_missing"))
      tryCatch({
        x <- rast(file)
        if (nlyr(x) != 1L || !nzchar(crs(x))) stop("Invalid saved-surface geometry/CRS.")
        list(values = extract(x, project(point_vector, crs(x)), ID = FALSE)[[1]], error = "")
      }, error = function(e) list(values = rep(NA_real_, nrow(points)), error = conditionMessage(e)))
    }
    target_cell <- function(cell) {
      rr <- rowFromCell(reference, cell); cc <- colFromCell(reference, cell)
      rast(nrows = 1L, ncols = 1L,
        xmin = xmin(reference) + (cc - 1) * res(reference)[1],
        xmax = xmin(reference) + cc * res(reference)[1],
        ymin = ymax(reference) - rr * res(reference)[2],
        ymax = ymax(reference) - (rr - 1) * res(reference)[2], crs = crs(reference))
    }
    to_ll <- function(xy, source_crs) project(xy, from = source_crs, to = "EPSG:4326")
    great_circle_km <- function(xy, lon, lat) {
      a <- xy[, 2] * pi / 180; b <- lat * pi / 180
      dlon <- (xy[, 1] - lon) * pi / 180
      h <- sin((a - b) / 2)^2 + cos(a) * cos(b) * sin(dlon / 2)^2
      6371.0088 * 2 * asin(sqrt(pmin(1, pmax(0, h))))
    }
    cat("[CHECK 2/3] Small native climate windows, separately for each model variable list.\n")
    total_jobs <- length(climate_dirs) * length(model_vars) * nrow(points)
    job <- 0L; summary_rows <- list()
    for (source_id in seq_along(climate_dirs)) {
      for (zone in names(model_vars)) {
        vars <- model_vars[[zone]]
        files <- file.path(climate_dirs[source_id], paste0(vars, ".tif"))
        current_source <- source_id
        h <- headers[kind == "raw_climate" & source_id == current_source & variable %in% vars]
        stack_error <- if (any(nzchar(h$error))) "raw_layer_read_error" else ""
        raw <- NULL
        if (!nzchar(stack_error)) {
          raw <- tryCatch({
            template <- rast(files[1])
            if (!all(vapply(files, function(f) isTRUE(compareGeom(rast(f), template,
                stopOnError = FALSE)), logical(1)))) stop("inconsistent_model_input_grids")
            x <- rast(files); names(x) <- vars; x
          }, error = function(e) { stack_error <<- conditionMessage(e); NULL })
        }
        saved_file <- file.path(saved_root, "normal", paste0("clim_suit_zone", zone, ".tif"))
        soil_saved_file <- file.path(soil_saved_root, paste0("soil_suit_zone", zone, ".tif"))
        all_inputs <- c(all_inputs, c(saved_file, soil_saved_file)[file.exists(c(saved_file, soil_saved_file))])
        saved <- get_saved(saved_file); saved_soil <- get_saved(soil_saved_file)
        raw_points <- if (!is.null(raw)) project(point_vector, crs(raw)) else NULL
        native_xy <- if (!is.null(raw_points)) crds(raw_points) else matrix(NA_real_, nrow(points), 2L)
        native_ids <- if (!is.null(raw)) cellFromXY(raw, native_xy) else rep(NA_real_, nrow(points))
        for (j in seq_len(nrow(points))) {
          job <- job + 1L
          cat(sprintf("  Window %d/%d: source %d, model zone %s, sample %s\n",
                      job, total_jobs, source_id, zone, points$sample_id[j])); flush.console()
          row <- data.table(source_id = source_id, raw_directory = climate_dirs[source_id],
            model_zone = as.integer(zone), sample_id = points$sample_id[j], point_type = points$point_type[j],
            longitude = points$longitude[j], latitude = points$latitude[j], reference_zone = points$reference_zone[j],
            n_selected_climate_variables = length(vars),
            soil_layers_inspected = soil_point$soil_layers_inspected[j],
            soil_layers_finite_descriptive = soil_point$soil_layers_finite[j],
            native_cell = native_ids[j], inside_raw_extent = is.finite(native_ids[j]),
            n_finite_at_native_cell = NA_integer_, all_selected_finite_at_native_cell = NA,
            window_radius_cells = radius, window_total_cells = NA_integer_, window_complete_cells = NA_integer_,
            window_reference_present_cells = NA_integer_, window_reference_present_climate_missing_cells = NA_integer_,
            nearest_complete_in_window_km_approx = NA_real_, nearest_complete_longitude = NA_real_,
            nearest_complete_latitude = NA_real_, nearest_complete_reference_zone = NA_real_,
            same_crs_as_reference = if (!is.null(raw)) same_crs(raw, reference) else NA,
            nearest_mask_value_at_reference_cell = NA_real_, bilinear_complete_mask_value_at_reference_cell = NA_real_,
            bilinear_dummy_available_at_reference_cell = NA,
            stored_climate_suitability = saved$values[j], stored_climate_available = is.finite(saved$values[j]),
            stored_soil_suitability = saved_soil$values[j], saved_climate_error = saved$error,
            saved_soil_error = saved_soil$error, support_comparison = "not_evaluated", error = stack_error)
          if (!nzchar(stack_error) && is.finite(native_ids[j])) {
            tryCatch({
              rr <- rowFromCell(raw, native_ids[j]); cc <- colFromCell(raw, native_ids[j])
              r1 <- max(1L, rr - radius); r2 <- min(nrow(raw), rr + radius)
              c1 <- max(1L, cc - radius); c2 <- min(ncol(raw), cc + radius)
              box <- ext(xmin(raw) + (c1 - 1) * res(raw)[1], xmin(raw) + c2 * res(raw)[1],
                         ymax(raw) - r2 * res(raw)[2], ymax(raw) - (r1 - 1) * res(raw)[2])
              win <- crop(raw, box, snap = "near")
              v <- values(win, mat = TRUE)
              finite_count <- rowSums(is.finite(v)); complete <- finite_count == length(vars)
              center <- cellFromXY(win, native_xy[j, , drop = FALSE])
              if (!is.finite(center)) stop("Local window does not contain target native cell.")
              row[, `:=`(n_finite_at_native_cell = finite_count[center],
                all_selected_finite_at_native_cell = complete[center],
                window_total_cells = nrow(v), window_complete_cells = sum(complete))]
              xy <- xyFromCell(win, seq_len(ncell(win)))
              ll_native <- to_ll(xy, crs(win))
              vp <- vect(as.data.frame(xy), geom = c("x", "y"), crs = crs(win))
              rz <- extract(reference, project(vp, crs(reference)), ID = FALSE)[[1]]
              row[, `:=`(window_reference_present_cells = sum(is.finite(rz)),
                window_reference_present_climate_missing_cells = sum(is.finite(rz) & !complete))]
              if (any(complete)) {
                good <- which(complete)
                d <- great_circle_km(ll_native[good, , drop = FALSE], points$longitude[j], points$latitude[j])
                nearest <- good[which.min(d)]
                row[, `:=`(nearest_complete_in_window_km_approx = min(d),
                  nearest_complete_longitude = ll_native[nearest, 1], nearest_complete_latitude = ll_native[nearest, 2],
                  nearest_complete_reference_zone = as.numeric(rz[nearest]))]
              }
              append_csv(data.table(source_id = source_id, model_zone = as.integer(zone),
                sample_id = points$sample_id[j], native_cell = cellFromXY(raw, xy),
                longitude = ll_native[, 1], latitude = ll_native[, 2],
                n_finite_selected_variables = finite_count, n_selected_variables = length(vars),
                complete_climate_inputs = complete, reference_zone_at_native_center = rz),
                "neighbourhood_coverage.csv")
              if (isTRUE(row$same_crs_as_reference)) {
                target <- target_cell(points$cell[j])
                binary <- rast(win[[1]]); values(binary) <- as.integer(complete)
                dummy <- binary; values(dummy) <- ifelse(complete, 1, NA_real_)
                near_value <- as.numeric(values(resample(binary, target, method = "near")))[1]
                fraction <- as.numeric(values(resample(binary, target, method = "bilinear")))[1]
                support <- is.finite(as.numeric(values(resample(dummy, target, method = "bilinear")))[1])
                row[, `:=`(nearest_mask_value_at_reference_cell = near_value,
                  bilinear_complete_mask_value_at_reference_cell = fraction,
                  bilinear_dummy_available_at_reference_cell = support)]
                diagnosis <- if (nzchar(saved$error)) "saved_surface_unavailable_for_comparison" else if (
                  support && !is.finite(saved$values[j])) "local_support_but_saved_missing_review_provenance_and_masks" else if (
                  !support && is.finite(saved$values[j])) "saved_finite_but_current_local_support_missing_review_provenance" else if (
                  !support) "raw_window_has_no_bilinear_complete_input_support" else if (
                  !complete[center]) "native_containing_cell_incomplete_but_bilinear_support_exists" else "raw_inputs_and_saved_surface_available"
                row[, support_comparison := diagnosis]
              } else row[, support_comparison := "CRS_differs_original_resample_path_requires_review"]
            }, error = function(e) {
              row[, `:=`(error = conditionMessage(e), support_comparison = "local_check_error")]
            })
          } else if (nzchar(stack_error)) row[, support_comparison := "raw_model_input_grid_or_read_error"] else {
            row[, support_comparison := "point_outside_raw_climate_extent"]
          }
          summary_rows[[length(summary_rows) + 1L]] <- row
          fwrite(rbindlist(summary_rows, fill = TRUE), file.path(out, "point_coverage_summary.csv"))
        }
      }
    }
    result <- rbindlist(summary_rows, fill = TRUE)
    fwrite(result[, .N, by = .(source_id, model_zone, point_type, support_comparison)],
      file.path(out, "diagnosis_counts.csv"))
    all_inputs <- unique(all_inputs)
    info <- file.info(all_inputs)
    fwrite(data.table(path = normalizePath(all_inputs, winslash = "/", mustWork = TRUE),
      size_bytes = as.numeric(info$size), modified = as.character(info$mtime)), file.path(out, "input_files.csv"))
    writeLines(c(
      "READ-ONLY RAW COVERAGE AND ALIGNMENT CHECK", paste("Run:", Sys.time()),
      paste("Previous trace:", trace_dir), paste("Saved sample points:", nrow(points)),
      paste("Native window radius:", radius, "cells; maximum width:", 2L * radius + 1L),
      paste("Soil input layers:", length(soil_files)), paste("Climate source copies inspected:", length(climate_dirs)),
      "No model was loaded, trained, or used to predict suitability.",
      "No raw data, cached suitability, assigned map, existing CSV, or manuscript was modified.",
      "Only this NEW report folder contains new output.", "",
      "INTERPRETATION:",
      "Climate completeness uses each model's own selected variables. The union is only an I/O list.",
      "Soil finite-layer counts describe every listed soil TIFF, not a model-specific input requirement.",
      "same_lattice compares CRS, resolution and origin; extents may differ without a lattice mismatch.",
      "The nearest complete native cell is searched ONLY within the stated window. Distance is approximate great-circle km.",
      "Native-cell reference zones are sampled at cell centers, not exact overlap areas.",
      "A bilinear 0/1 mask is a local coverage indicator, NOT suitability or a model prediction.",
      "The finite/NA dummy tests local interpolation support, NOT whether a cached prediction is definitively wrong.",
      "Mismatch flags require checking source version, native masks and cache provenance.",
      "No missing value is filled; nearest cells are diagnostic only.",
      "These sampled windows do not estimate national frequencies, areas, or the cause of all 19,506 legacy cells.",
      "No water/coastal explanation is asserted without independent source masks.", "",
      capture.output(sessionInfo())), file.path(out, "run_info.txt"))
    cat("[CHECK 3/3 SAVED]", normalizePath(out, winslash = "/"), "\n")
    cat("No model, project raster, or existing analysis result was overwritten.\n")
    cat("For review, share this folder's CSV files and run_info.txt; no TIFF upload is needed.\n")
  }
  run_check()
})
