# Climate, soil and dual-suitability prediction for the final Multi-Forest.
# Established raster paths are retained and valid outputs are reused by default.

library(data.table)
library(randomForest)
library(terra)

rm(list = ls())
gc()

find_project_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "script", "1. data.R"))) return(path)
    parent <- dirname(path)
    if (parent == path) stop("Run inside the repository or set ECOCHINA2_DIR.")
    path <- parent
  }
}

env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = if (default) "true" else "false")
  tolower(trimws(value)) %in% c("1", "true", "yes", "y")
}

env_root <- Sys.getenv("ECOCHINA2_DIR", unset = "")
base_dir <- if (nzchar(env_root)) {
  normalizePath(env_root, winslash = "/", mustWork = TRUE)
} else {
  find_project_root()
}

reference_file <- file.path(base_dir, "raster", "ecosys_ori.tif")
if (!file.exists(reference_file)) stop("Run script/1. data.R first.")
reference_map <- rast(reference_file)

env_climate <- Sys.getenv("ECOCHINA2_CLIMATE_DIR", unset = "")
climate_root <- if (nzchar(env_climate)) {
  env_climate
} else {
  file.path(base_dir, "data", "rasters", "climate")
}

env_soil <- Sys.getenv("ECOCHINA2_SOIL_RASTER_DIR", unset = "")
soil_raster_dir <- if (nzchar(env_soil)) {
  env_soil
} else {
  file.path(base_dir, "data", "rasters", "soil")
}

climate_model_dir <- file.path(base_dir, "rf")
soil_model_dir <- file.path(base_dir, "rf_soil")
climate_root_out <- file.path(base_dir, "clim suitability", "mf_var")
soil_out <- file.path(base_dir, "soil suitability", "plain_mf", "normal")
dual_root_out <- file.path(base_dir, "dual suit", "mf_var")
sensitivity_root <- file.path(base_dir, "dual suit sensitivity")
assessment_dir <- file.path(base_dir, "assessment")

dir.create(climate_root_out, recursive = TRUE, showWarnings = FALSE)
dir.create(soil_out, recursive = TRUE, showWarnings = FALSE)
dir.create(dual_root_out, recursive = TRUE, showWarnings = FALSE)
dir.create(assessment_dir, recursive = TRUE, showWarnings = FALSE)

tmp_dir <- file.path(base_dir, "tmp_prediction_var")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
terraOptions(tempdir = tmp_dir, memfrac = 0.15)

zoneID <- c(1:7, 9:50, 52:55)
soil_gate <- 0.20
soil_gates <- c(0, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50)
novel_threshold <- 0.40

detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(detected_cores)) detected_cores <- 1L
predict_cores <- min(4L, max(1L, detected_cores - 1L))

reuse_suitability_requested <- env_flag("ECOCHINA2_REUSE_SUITABILITY", TRUE)
force_suitability <- env_flag("ECOCHINA2_FORCE_SUITABILITY", FALSE) ||
  !reuse_suitability_requested
force_dual <- env_flag("ECOCHINA2_FORCE_DUAL", FALSE) || force_suitability
reuse_suitability <- reuse_suitability_requested && !force_suitability
reuse_dual <- env_flag("ECOCHINA2_REUSE_DUAL", TRUE) && !force_dual
force_sensitivity <- env_flag("ECOCHINA2_FORCE_SENSITIVITY", FALSE)
write_gate_rasters <- env_flag("ECOCHINA2_WRITE_GATE_RASTERS", FALSE)

scenario_table <- data.table(
  scenario = c(
    "normal",
    "2011-2040SSP245", "2041-2070SSP245", "2071-2100SSP245",
    "2011-2040SSP585", "2041-2070SSP585", "2071-2100SSP585"
  ),
  climate_folder = c(
    "Normal_1961_1990",
    "8GCMs_ensemble_ssp245_2011-2040",
    "8GCMs_ensemble_ssp245_2041-2070",
    "8GCMs_ensemble_ssp245_2071-2100",
    "8GCMs_ensemble_ssp585_2011-2040",
    "8GCMs_ensemble_ssp585_2041-2070",
    "8GCMs_ensemble_ssp585_2071-2100"
  )
)

valid_raster <- function(file, template = reference_map) {
  if (!file.exists(file) || is.na(file.info(file)$size) || file.info(file)$size <= 0) {
    return(FALSE)
  }
  tryCatch({
    x <- rast(file)
    nlyr(x) == 1L && compareGeom(x, template, stopOnError = FALSE)
  }, error = function(e) FALSE)
}

remove_raster_files <- function(file) {
  stem <- tools::file_path_sans_ext(file)
  targets <- unique(c(
    file, paste0(file, c(".aux.xml", ".ovr", ".msk")),
    paste0(stem, c(".aux.xml", ".ovr", ".tfw"))
  ))
  for (attempt in seq_len(5L)) {
    suppressWarnings(unlink(targets[file.exists(targets)], force = TRUE))
    if (!file.exists(file)) return(invisible(TRUE))
    gc()
    Sys.sleep(0.2 * attempt)
  }
  if (file.exists(file)) stop("Cannot replace locked raster: ", file)
  invisible(TRUE)
}

load_model <- function(file, object_name) {
  if (!file.exists(file)) stop("Missing model: ", file)
  e <- new.env()
  load(file, envir = e)
  if (!exists(object_name, envir = e, inherits = FALSE)) {
    stop("Object ", object_name, " not found in: ", file)
  }
  get(object_name, envir = e, inherits = FALSE)
}

get_varlist <- function(model, model_file) {
  vars <- model$varlist
  if (is.null(vars) || !length(vars)) vars <- rownames(model$importance)
  if (is.null(vars) || !length(vars)) stop("No predictors found in: ", model_file)
  if (length(vars) != 9L) {
    stop("Expected 9 selected predictors in ", model_file, "; found ", length(vars), ".")
  }
  as.character(vars)
}

get_stack <- function(varlist, raster_dir) {
  files <- file.path(raster_dir, paste0(varlist, ".tif"))
  missing <- files[!file.exists(files)]
  if (length(missing)) {
    stop("Missing predictor raster(s):\n", paste(missing, collapse = "\n"))
  }
  x <- rast(files)
  names(x) <- varlist
  x
}

probability_one <- function(model, data) {
  probability <- predict(model, as.data.frame(data), type = "prob")
  if (!"1" %in% colnames(probability)) stop("Model has no class '1'.")
  as.numeric(probability[, "1"])
}

predict_suitability <- function(model_file, object_name, raster_dir,
                                output_file, label) {
  if (reuse_suitability && valid_raster(output_file)) {
    cat("[USE SUITABILITY]", output_file, "\n")
    return(rast(output_file))
  }
  if (!dir.exists(raster_dir)) stop("Missing predictor directory: ", raster_dir)

  model <- load_model(model_file, object_name)
  varlist <- get_varlist(model, model_file)
  predictors <- get_stack(varlist, raster_dir)
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(output_file)) remove_raster_files(output_file)

  if (compareGeom(predictors[[1]], reference_map, stopOnError = FALSE)) {
    prediction <- terra::predict(
      predictors, model, fun = probability_one,
      filename = output_file, overwrite = TRUE, cores = predict_cores,
      cpkgs = "randomForest", na.rm = TRUE,
      wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW")
    )
  } else {
    temporary_file <- tempfile("prediction_", tmpdir = tmp_dir, fileext = ".tif")
    prediction_native <- terra::predict(
      predictors, model, fun = probability_one,
      filename = temporary_file, overwrite = TRUE, cores = predict_cores,
      cpkgs = "randomForest", na.rm = TRUE,
      wopt = list(datatype = "FLT4S")
    )
    prediction <- resample(
      prediction_native, reference_map, method = "bilinear",
      filename = output_file, overwrite = TRUE,
      wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW")
    )
    rm(prediction_native)
    remove_raster_files(temporary_file)
  }
  names(prediction) <- label
  rm(model, predictors)
  gc()
  rast(output_file)
}

dual_suitability <- function(climate, soil, gate) {
  if (gate == 0) return(climate)
  ifel(is.na(soil), NA, ifel(soil > gate, climate, 0))
}

sensitivity_file <- file.path(assessment_dir, "soil_gate_sensitivity.csv")
sensitivity_complete <- function(file) {
  if (!file.exists(file)) return(FALSE)
  tryCatch({
    x <- fread(file)
    required <- c("scenario", "soil_gate", "valid_pixels", "mean_max_dual")
    if (length(setdiff(required, names(x)))) return(FALSE)
    observed <- unique(x[, .(scenario, soil_gate = round(soil_gate, 2))])
    expected <- CJ(
      scenario = scenario_table$scenario,
      soil_gate = round(soil_gates, 2),
      unique = TRUE
    )
    nrow(observed) == nrow(expected) &&
      nrow(merge(expected, observed, by = c("scenario", "soil_gate"), all = FALSE)) ==
        nrow(expected)
  }, error = function(e) FALSE)
}

primary_suitability_files <- c(
  file.path(soil_out, paste0("soil_suit_zone", zoneID, ".tif")),
  unlist(lapply(seq_len(nrow(scenario_table)), function(job) {
    file.path(
      climate_root_out,
      scenario_table$scenario[job],
      paste0("clim_suit_zone", zoneID, ".tif")
    )
  }))
)
alternative_gates <- soil_gates[soil_gates != soil_gate]
gate_raster_jobs <- rbindlist(lapply(alternative_gates, function(gate) {
  rbindlist(lapply(scenario_table$scenario, function(scenario) {
    data.table(
      output = file.path(
        sensitivity_root,
        paste0("gate_", gsub("\\.", "p", sprintf("%.2f", gate))),
        scenario,
        paste0("dual_suitability_zone", zoneID, ".tif")
      ),
      climate = file.path(
        climate_root_out, scenario, paste0("clim_suit_zone", zoneID, ".tif")
      ),
      soil = file.path(soil_out, paste0("soil_suit_zone", zoneID, ".tif"))
    )
  }))
}))
gate_raster_files <- gate_raster_jobs$output
gate_rasters_current <- FALSE
if (all(file.exists(unlist(gate_raster_jobs)))) {
  gate_info <- file.info(gate_raster_jobs$output)
  climate_info <- file.info(gate_raster_jobs$climate)
  soil_info <- file.info(gate_raster_jobs$soil)
  gate_rasters_current <-
    all(is.finite(gate_info$size) & gate_info$size > 0) &&
      all(
        as.numeric(gate_info$mtime) >= pmax(
          as.numeric(climate_info$mtime), as.numeric(soil_info$mtime)
        )
      )
}

sensitivity_current <- sensitivity_complete(sensitivity_file)
if (sensitivity_current && all(file.exists(primary_suitability_files))) {
  sensitivity_current <- file.info(sensitivity_file)$mtime >=
    max(file.info(primary_suitability_files)$mtime)
}

run_sensitivity <- force_suitability || force_sensitivity ||
  !isTRUE(sensitivity_current) ||
  (write_gate_rasters && !gate_rasters_current)

# Ensure static soil suitability only when a missing primary dual raster or the
# sensitivity analysis needs it.
dual_files_all <- unlist(lapply(scenario_table$scenario, function(scenario) {
  file.path(
    dual_root_out,
    scenario,
    paste0("dual_suitability_zone", zoneID, ".tif")
  )
}))
need_primary_dual <- !reuse_dual || !all(vapply(dual_files_all, valid_raster, logical(1)))
need_soil <- need_primary_dual || run_sensitivity

if (need_soil) {
  for (zone in zoneID) {
    predict_suitability(
      file.path(soil_model_dir, paste0("soil_mf_zone", zone, ".Rdata")),
      "soil_mf",
      soil_raster_dir,
      file.path(soil_out, paste0("soil_suit_zone", zone, ".tif")),
      "soil_suitability"
    )
  }
} else {
  cat("[USE EXISTING] Primary dual rasters and sensitivity table are complete.\n")
}

# Climate and primary dual suitability. Each valid dual raster is skipped before
# its climate model or static-soil raster is loaded unless sensitivity needs it.
for (job in seq_len(nrow(scenario_table))) {
  scenario <- scenario_table$scenario[job]
  climate_dir <- file.path(climate_root, scenario_table$climate_folder[job])
  climate_dir_out <- file.path(climate_root_out, scenario)
  dual_dir_out <- file.path(dual_root_out, scenario)
  dir.create(climate_dir_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(dual_dir_out, recursive = TRUE, showWarnings = FALSE)

  for (zone in zoneID) {
    climate_file <- file.path(climate_dir_out, paste0("clim_suit_zone", zone, ".tif"))
    dual_file <- file.path(dual_dir_out, paste0("dual_suitability_zone", zone, ".tif"))
    dual_ok <- reuse_dual && valid_raster(dual_file)

    if (dual_ok && !run_sensitivity) {
      cat("[USE DUAL]", scenario, "| zone", zone, "\n")
      next
    }

    climate <- predict_suitability(
      file.path(climate_model_dir, paste0("clm_mfVar_zone", zone, ".Rdata")),
      "clm_mfVar",
      climate_dir,
      climate_file,
      if (scenario == "normal") "climate_suitability" else "future_climate_suitability"
    )

    if (dual_ok) {
      cat("[USE DUAL]", scenario, "| zone", zone, "\n")
      rm(climate)
      next
    }

    soil_file <- file.path(soil_out, paste0("soil_suit_zone", zone, ".tif"))
    if (!valid_raster(soil_file)) stop("Missing valid soil suitability: ", soil_file)
    soil <- rast(soil_file)
    if (!compareGeom(climate, soil, stopOnError = FALSE)) {
      soil <- resample(soil, climate, method = "bilinear")
    }

    dual <- dual_suitability(climate, soil, soil_gate)
    names(dual) <- if (scenario == "normal") {
      "dual_suitability"
    } else {
      "future_dual_suitability"
    }
    if (file.exists(dual_file)) remove_raster_files(dual_file)
    writeRaster(
      dual, dual_file, overwrite = TRUE,
      wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW")
    )
    cat("[SAVED DUAL]", scenario, "| zone", zone, "\n")
    rm(climate, soil, dual)
    gc()
  }
}

# Soil-gate sensitivity: one blockwise pass per scenario, only when the complete
# 7-scenario x 8-gate table is absent or explicitly forced.
if (!run_sensitivity) {
  cat("[USE EXISTING]", sensitivity_file, "\n")
} else {
  soil_files <- file.path(soil_out, paste0("soil_suit_zone", zoneID, ".tif"))
  if (!all(vapply(soil_files, valid_raster, logical(1)))) {
    stop("Static soil suitability stack is incomplete.")
  }
  soil_stack <- rast(soil_files)
  names(soil_stack) <- as.character(zoneID)

  analysis_mask <- subst(
    reference_map,
    from = zoneID,
    to = rep(1, length(zoneID)),
    others = NA
  )
  cell_area <- mask(cellSize(reference_map, unit = "km"), analysis_mask)
  gate_results <- vector("list", nrow(scenario_table))

  for (job in seq_len(nrow(scenario_table))) {
    scenario <- scenario_table$scenario[job]
    climate_files <- file.path(
      climate_root_out,
      scenario,
      paste0("clim_suit_zone", zoneID, ".tif")
    )
    if (!all(vapply(climate_files, valid_raster, logical(1)))) {
      stop("Climate suitability stack is incomplete for: ", scenario)
    }
    climate_stack <- rast(climate_files)
    names(climate_stack) <- as.character(zoneID)
    input <- c(climate_stack, soil_stack, reference_map, cell_area)
    nz <- length(zoneID)

    result <- data.table(
      soil_gate = soil_gates,
      valid_pixels = 0, analogue_pixels = 0, below_0.4_pixels = 0,
      valid_area_km2 = 0, analogue_area_km2 = 0, below_0.4_area_km2 = 0,
      sum_max_dual = 0
    )
    rows_per_block <- max(1L, floor(1e7 / (ncol(input) * nlyr(input))))
    block_starts <- seq.int(1L, nrow(input), by = rows_per_block)
    readStart(input)

    for (start_row in block_starts) {
      n_rows <- min(rows_per_block, nrow(input) - start_row + 1L)
      values <- readValues(input, row = start_row, nrows = n_rows, mat = TRUE)
      climate_values <- values[, seq_len(nz), drop = FALSE]
      soil_values <- values[, nz + seq_len(nz), drop = FALSE]
      original <- values[, 2L * nz + 1L]
      area <- values[, 2L * nz + 2L]
      modeled <- original %in% zoneID & is.finite(area)

      for (gate_index in seq_along(soil_gates)) {
        gate <- soil_gates[gate_index]
        candidates <- climate_values
        if (gate > 0) {
          candidates[!is.na(soil_values) & soil_values <= gate] <- 0
          candidates[is.na(soil_values)] <- -Inf
        }
        candidates[!is.finite(candidates)] <- -Inf
        max_col <- max.col(candidates, ties.method = "first")
        max_dual <- candidates[cbind(seq_len(nrow(candidates)), max_col)]
        valid <- modeled & is.finite(max_dual)
        analogue <- valid & max_dual >= novel_threshold
        below <- valid & max_dual < novel_threshold

        result[gate_index, `:=`(
          valid_pixels = valid_pixels + sum(valid),
          analogue_pixels = analogue_pixels + sum(analogue),
          below_0.4_pixels = below_0.4_pixels + sum(below),
          valid_area_km2 = valid_area_km2 + sum(area[valid]),
          analogue_area_km2 = analogue_area_km2 + sum(area[analogue]),
          below_0.4_area_km2 = below_0.4_area_km2 + sum(area[below]),
          sum_max_dual = sum_max_dual + sum(max_dual[valid])
        )]
      }
    }
    readStop(input)

    result[, `:=`(
      scenario = scenario,
      climate_only = soil_gate == 0,
      novel_rule_applied = scenario != "normal",
      analogue_share = analogue_pixels / valid_pixels,
      below_0.4_share = below_0.4_pixels / valid_pixels,
      mean_max_dual = sum_max_dual / valid_pixels
    )]
    result[, sum_max_dual := NULL]
    setcolorder(result, c(
      "scenario", "soil_gate", "climate_only", "novel_rule_applied",
      "valid_pixels", "analogue_pixels", "below_0.4_pixels",
      "valid_area_km2", "analogue_area_km2", "below_0.4_area_km2",
      "analogue_share", "below_0.4_share", "mean_max_dual"
    ))
    gate_results[[job]] <- result

    if (write_gate_rasters) {
      for (gate in soil_gates[soil_gates != soil_gate]) {
        gate_dir <- file.path(
          sensitivity_root,
          paste0("gate_", gsub("\\.", "p", sprintf("%.2f", gate))),
          scenario
        )
        dir.create(gate_dir, recursive = TRUE, showWarnings = FALSE)
        gate_outputs <- file.path(
          gate_dir, paste0("dual_suitability_zone", zoneID, ".tif")
        )
        gate_outputs_current <- vapply(seq_along(zoneID), function(z) {
          valid_raster(gate_outputs[z]) &&
            isTRUE(file.info(gate_outputs[z])$mtime >= max(file.info(c(
              climate_files[z], soil_files[z]
            ))$mtime))
        }, logical(1))
        if (reuse_dual && all(gate_outputs_current)) {
          cat("[USE GATE RASTERS]", scenario, "| gate", gate, "\n")
          next
        }

        gate_stack <- dual_suitability(climate_stack, soil_stack, gate)
        for (z in seq_along(zoneID)) {
          output <- gate_outputs[z]
          if (!reuse_dual || !gate_outputs_current[z]) {
            writeRaster(
              gate_stack[[z]], output, overwrite = TRUE,
              wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW")
            )
          }
        }
        rm(gate_stack)
      }
    }
    cat("[SOIL-GATE SUMMARY]", scenario, "\n")
    rm(climate_stack, input, result)
    gc()
  }

  gate_summary <- rbindlist(gate_results)
  gate_summary[, scenario_order := match(scenario, scenario_table$scenario)]
  setorder(gate_summary, scenario_order, soil_gate)
  gate_summary[, scenario_order := NULL]
  fwrite(gate_summary, sensitivity_file)
  cat("[SAVED]", sensitivity_file, "\n")
}

cat(
  "\nCOMPLETE\n",
  "Method: mf_var (selected-variable Multi-Forest)\n",
  "Static soil source: plain_mf (historical folder name)\n",
  "Primary soil gate: ", soil_gate, "\n",
  "Climate outputs: ", climate_root_out, "\n",
  "Dual outputs: ", dual_root_out, "\n",
  sep = ""
)
