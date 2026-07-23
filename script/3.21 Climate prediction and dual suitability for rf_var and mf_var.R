# 3.21 Climate prediction and dual suitability for rf_var and mf_var
# ==============================================================================
# This supplemental script runs only the two climate models created by script 2.16:
#   rf_var: clm_rfVar_zone*.Rdata
#   mf_var: clm_mfVar_zone*.Rdata
#
# The soil models are unchanged:
#   rf_var uses the existing plain_rf soil suitability rasters.
#   mf_var uses the existing plain_mf soil suitability rasters.
#
# Outputs are written only to the new method folders:
#   clim suitability/rf_var/
#   clim suitability/mf_var/
#   dual suit/rf_var/
#   dual suit/mf_var/
#
# Existing plain_rf, plain_mf, optimized_rf and optimized_mf outputs are untouched.
# ==============================================================================

library(terra)
library(randomForest)
library(data.table)

rm(list = ls())
gc()

base_dir <- "H:/Jing/ecoChina2"

reference_file <- file.path(base_dir, "raster/ecosys_ori.tif")
r <- rast(reference_file)

climdir_normal <- "H:/Jing/ecoChina/play/China/ClimateData/CN/800m/Normal_1961_1990"
clm_mod_dir <- file.path(base_dir, "rf")

tmp_dir <- file.path(base_dir, "tmp_prediction_var")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
terraOptions(tempdir = tmp_dir, memfrac = 0.15)

zones <- c(1:7, 9:30, 31:50, 52:55)
soil_threshold <- 0.2
reuse_existing_climate <- TRUE
overwrite_dual <- TRUE
predict_cores <- min(4L, max(1L, parallel::detectCores() - 1L))

model_set <- data.table(
  method = c("rf_var", "mf_var"),
  clm_prefix = c("clm_rfVar_zone", "clm_mfVar_zone"),
  clm_object = c("clm_rfVar", "clm_mfVar"),
  soil_source_method = c("plain_rf", "plain_mf")
)

remove_raster_files <- function(filepath, retries = 5L, wait = 0.5,
                                stop_if_failed = TRUE) {
  base <- tools::file_path_sans_ext(filepath)
  
  files <- unique(c(
    filepath,
    paste0(filepath, ".aux.xml"),
    paste0(filepath, ".ovr"),
    paste0(filepath, ".msk"),
    paste0(base, ".aux.xml"),
    paste0(base, ".ovr"),
    paste0(base, ".tfw")
  ))
  
  files <- files[file.exists(files)]
  if (length(files) == 0) return(invisible(TRUE))
  
  remaining <- files
  
  for (attempt in seq_len(retries)) {
    gc()
    suppressWarnings(unlink(remaining, force = TRUE))
    remaining <- files[file.exists(files)]
    
    if (length(remaining) == 0) return(invisible(TRUE))
    if (attempt < retries) Sys.sleep(wait)
  }
  
  msg <- paste0(
    "Cannot remove existing raster file. The file may still be open or locked:\n",
    paste(remaining, collapse = "\n")
  )
  
  if (stop_if_failed) stop(msg, call. = FALSE)
  warning(msg, call. = FALSE)
  invisible(FALSE)
}

get_stack <- function(varlist, raster_dir) {
  tif_paths <- file.path(raster_dir, paste0(varlist, ".tif"))
  
  if (!all(file.exists(tif_paths))) {
    missing_files <- tif_paths[!file.exists(tif_paths)]
    cat("  [WARN] missing tif:",
        paste(basename(missing_files), collapse = ", "), "\n")
    return(NULL)
  }
  
  s <- rast(tif_paths)
  names(s) <- varlist
  s
}

load_rf <- function(file, object_name) {
  e <- new.env()
  load(file, envir = e)
  
  if (!exists(object_name, envir = e)) {
    stop("Object ", object_name, " not found in model file: ", file)
  }
  
  get(object_name, envir = e)
}

get_varlist <- function(model) {
  if (!is.null(model$varlist) && length(model$varlist) > 0) {
    return(model$varlist)
  }
  
  if (!is.null(model$importance) && nrow(model$importance) > 0) {
    varlist <- rownames(model$importance)
    if (!is.null(varlist) && length(varlist) > 0) return(varlist)
  }
  
  stop("No predictor names found in model object.")
}

rf_prob1 <- function(model, data) {
  data <- as.data.frame(data)
  predict(model, data, type = "prob")[, "1"]
}

mcPredict_terra <- function(x, y, filename, compress = FALSE,
                            cores = predict_cores, ...) {
  remove_raster_files(filename)
  
  predict(
    x,
    y,
    fun = rf_prob1,
    filename = filename,
    overwrite = TRUE,
    cores = cores,
    cpkgs = "randomForest",
    na.rm = TRUE,
    wopt = list(
      datatype = "FLT4S",
      gdal = if (compress) "COMPRESS=LZW" else "COMPRESS=NONE"
    ),
    ...
  )
}

calculate_dual_suitability <- function(climate, soil,
                                       threshold = soil_threshold) {
  ifel(
    is.na(soil),
    NA,
    ifel(soil > threshold, climate, 0)
  )
}

predict_one_climate_model <- function(model_file, object_name, raster_dir,
                                      out_file, template, label) {
  if (file.exists(out_file) && reuse_existing_climate) {
    p <- rast(out_file)
    same_geom <- compareGeom(p, template, stopOnError = FALSE)
    
    cat("  existing ", label, " geometry matches template: ",
        same_geom, "\n", sep = "")
    
    if (same_geom) {
      cat("  [USE EXISTING]", out_file, "\n")
      return(p)
    }
    
    rm(p)
    gc()
  }
  
  if (!file.exists(model_file)) {
    cat("[SKIP] model not found:", model_file, "\n")
    return(NULL)
  }
  
  m <- load_rf(model_file, object_name)
  rfVar <- get_varlist(m)
  
  cat("  ", label, " vars: ", length(rfVar), "\n", sep = "")
  
  s <- get_stack(rfVar, raster_dir)
  if (is.null(s)) {
    rm(m)
    gc()
    return(NULL)
  }
  
  same_geom <- compareGeom(s[[1]], template, stopOnError = FALSE)
  cat("  geometry matches template:", same_geom, "\n")
  
  if (same_geom) {
    p <- mcPredict_terra(
      s,
      m,
      filename = out_file,
      compress = TRUE
    )
  } else {
    tmp_file <- tempfile(
      pattern = paste0(
        "tmp_", label, "_",
        tools::file_path_sans_ext(basename(out_file)), "_"
      ),
      tmpdir = tmp_dir,
      fileext = ".tif"
    )
    
    on.exit(
      suppressWarnings(
        remove_raster_files(tmp_file, stop_if_failed = FALSE)
      ),
      add = TRUE
    )
    
    p0 <- mcPredict_terra(
      s,
      m,
      filename = tmp_file,
      compress = FALSE
    )
    
    remove_raster_files(out_file)
    
    p <- resample(
      p0,
      template,
      method = "bilinear",
      filename = out_file,
      overwrite = TRUE,
      wopt = list(
        datatype = "FLT4S",
        gdal = "COMPRESS=LZW"
      )
    )
    
    rm(p0)
    gc()
    remove_raster_files(tmp_file)
  }
  
  names(p) <- label
  rm(m, s, p)
  gc()
  
  rast(out_file)
}

write_dual <- function(climate, soil_file, dual_file, template, label) {
  if (!file.exists(soil_file)) {
    cat("[SKIP] existing soil suitability not found:", soil_file, "\n")
    return(FALSE)
  }
  
  psoil <- rast(soil_file)
  psoil <- mask(psoil, template)
  
  if (!compareGeom(climate, psoil, stopOnError = FALSE)) {
    psoil <- resample(psoil, climate, method = "bilinear")
  }
  
  dual <- calculate_dual_suitability(
    climate,
    psoil,
    threshold = soil_threshold
  )
  
  names(dual) <- label
  
  if (overwrite_dual) remove_raster_files(dual_file)
  
  writeRaster(
    dual,
    dual_file,
    overwrite = overwrite_dual,
    wopt = list(
      datatype = "FLT4S",
      gdal = "COMPRESS=LZW"
    )
  )
  
  rm(psoil, dual)
  gc()
  TRUE
}

# 1. Normal-period climate and dual suitability =================================

for (mrow in seq_len(nrow(model_set))) {
  method <- model_set$method[mrow]
  soil_source_method <- model_set$soil_source_method[mrow]
  
  cat(
    "\n==============================\n",
    "MODEL VERSION: ", method, "\n",
    "SOIL SOURCE: ", soil_source_method, "\n",
    "==============================\n",
    sep = ""
  )
  
  out_dir_clm <- file.path(base_dir, "clim suitability", method, "normal")
  out_dir_dual <- file.path(base_dir, "dual suit", method, "normal")
  soil_dir <- file.path(
    base_dir,
    "soil suitability",
    soil_source_method,
    "normal"
  )
  
  dir.create(out_dir_clm, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir_dual, recursive = TRUE, showWarnings = FALSE)
  
  for (i in zones) {
    cat("\n[RUN]", method, "normal | zone", i, "\n")
    
    clim_file <- file.path(
      out_dir_clm,
      paste0("clim_suit_zone", i, ".tif")
    )
    
    soil_file <- file.path(
      soil_dir,
      paste0("soil_suit_zone", i, ".tif")
    )
    
    dual_file <- file.path(
      out_dir_dual,
      paste0("dual_suitability_zone", i, ".tif")
    )
    
    pclim <- predict_one_climate_model(
      model_file = file.path(
        clm_mod_dir,
        paste0(model_set$clm_prefix[mrow], i, ".Rdata")
      ),
      object_name = model_set$clm_object[mrow],
      raster_dir = climdir_normal,
      out_file = clim_file,
      template = r,
      label = "climate_suitability"
    )
    
    if (is.null(pclim)) next
    
    completed <- write_dual(
      climate = pclim,
      soil_file = soil_file,
      dual_file = dual_file,
      template = r,
      label = "dual_suitability"
    )
    
    if (completed) cat("[DONE]", method, "normal | zone", i, "\n")
    
    rm(pclim)
    gc()
  }
}

# 2. Future climate and dual suitability ========================================

scenarios <- c("ssp245", "ssp585")
time_periods <- c("2011-2040", "2041-2070", "2071-2100")

for (mrow in seq_len(nrow(model_set))) {
  method <- model_set$method[mrow]
  soil_source_method <- model_set$soil_source_method[mrow]
  
  soil_dir <- file.path(
    base_dir,
    "soil suitability",
    soil_source_method,
    "normal"
  )
  
  for (scenario in scenarios) {
    for (period in time_periods) {
      scen_name <- paste0(period, toupper(scenario))
      
      cat("\n===", method, "future:", scen_name, "===\n")
      
      climdir_fut <- paste0(
        "H:/Jing/ecoChina/play/China/ClimateData/CN/800m/8GCMs_ensemble_",
        scenario,
        "_",
        period
      )
      
      out_dir_clm_fut <- file.path(
        base_dir,
        "clim suitability",
        method,
        scen_name
      )
      
      out_dir_dual_fut <- file.path(
        base_dir,
        "dual suit",
        method,
        scen_name
      )
      
      dir.create(out_dir_clm_fut, recursive = TRUE, showWarnings = FALSE)
      dir.create(out_dir_dual_fut, recursive = TRUE, showWarnings = FALSE)
      
      for (i in zones) {
        cat("\n[RUN]", method, scen_name, "| zone", i, "\n")
        
        clim_file <- file.path(
          out_dir_clm_fut,
          paste0("clim_suit_zone", i, ".tif")
        )
        
        soil_file <- file.path(
          soil_dir,
          paste0("soil_suit_zone", i, ".tif")
        )
        
        dual_file <- file.path(
          out_dir_dual_fut,
          paste0("dual_suitability_zone", i, ".tif")
        )
        
        pclim <- predict_one_climate_model(
          model_file = file.path(
            clm_mod_dir,
            paste0(model_set$clm_prefix[mrow], i, ".Rdata")
          ),
          object_name = model_set$clm_object[mrow],
          raster_dir = climdir_fut,
          out_file = clim_file,
          template = r,
          label = "future_climate_suitability"
        )
        
        if (is.null(pclim)) next
        
        completed <- write_dual(
          climate = pclim,
          soil_file = soil_file,
          dual_file = dual_file,
          template = r,
          label = "future_dual_suitability"
        )
        
        if (completed) cat("[DONE]", method, scen_name, "| zone", i, "\n")
        
        rm(pclim)
        gc()
      }
    }
  }
}

cat(
  "\nCOMPLETE\n",
  "Methods: ", paste(model_set$method, collapse = ", "), "\n",
  "Climate suitability root: ",
  file.path(base_dir, "clim suitability"), "\n",
  "Dual suitability root: ",
  file.path(base_dir, "dual suit"), "\n",
  sep = ""
)
