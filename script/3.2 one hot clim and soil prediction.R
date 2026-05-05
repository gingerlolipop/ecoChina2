# Dual Suitability Prediction — Part 1: Normal dual suitability
# ============================================================
library(terra)
library(randomForest)
library(CEMT)

setwd("H:/Jing/ecoChina/play/China")

original_raster <- rast('raster/veg_3')

climdir_normal <- 'H:/Jing/ecoChina/play/China/ClimateData/CN/800m/Normal_1961_1990'
soilDir        <- 'H:/Jing/soil rasters/tif2'

clm_mod_dir    <- 'H:/Jing/ecoChina2/rf_final'
soil_mod_dir   <- 'H:/Jing/ecoChina2/rf_final_soil'

out_dir_clm    <- 'H:/Jing/ecoChina2/clim suitability/normal'
out_dir_soil   <- 'H:/Jing/ecoChina2/soil suitability/normal'
out_dir_dual   <- 'H:/Jing/ecoChina2/dual suit/normal'
tmp_dir        <- 'H:/Jing/ecoChina2/tmp_prediction'

dir.create(out_dir_clm,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_soil, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_dual, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir,      recursive = TRUE, showWarnings = FALSE)

terraOptions(tempdir = tmp_dir, memfrac = 0.15)

zones     <- setdiff(1:55, 8)
threshold <- 0.2


# ── helper: remove old temporary raster files ────────────────────────────────
remove_tmp_files <- function(filepath) {
  base <- tools::file_path_sans_ext(filepath)
  files <- c(filepath,
             paste0(base, ".aux.xml"),
             paste0(base, ".ovr"),
             paste0(base, ".tfw"))
  files <- files[file.exists(files)]
  if (length(files) > 0) unlink(files, force = TRUE)
}


# ── helper: predict rfEnsemble onto SpatRaster stack by small row chunks ─────
predict_ensemble_raster <- function(ens, raster_stack, filename,
                                    overwrite = TRUE, chunk_nrows = 50L) {
  
  remove_tmp_files(filename)
  
  out        <- raster_stack[[1]]
  names(out) <- "suitability"
  
  # Important: do not assign writeStart() back to out.
  # writeStart() returns write/block metadata; out must remain a SpatRaster.
  writeStart(out, filename = filename, overwrite = overwrite)
  
  nr <- nrow(raster_stack)
  
  readStart(raster_stack)
  
  for (row_start in seq(1, nr, by = chunk_nrows)) {
    
    nrows_i <- min(chunk_nrows, nr - row_start + 1L)
    
    vals <- readValues(raster_stack,
                       row   = row_start,
                       nrows = nrows_i,
                       mat   = TRUE)
    
    vals[vals == -9999] <- NA
    
    vals <- as.data.frame(vals)
    names(vals) <- names(raster_stack)
    vals <- vals[, ens$varlist, drop = FALSE]
    
    complete_idx <- complete.cases(vals)
    pred_vals    <- rep(NA_real_, nrow(vals))
    
    if (sum(complete_idx) > 0) {
      X        <- vals[complete_idx, , drop = FALSE]
      pred_mat <- do.call(cbind, lapply(
        ens$models,
        function(m) predict(m, X, type = "prob")[, "1"]
      ))
      pred_vals[complete_idx] <- as.numeric(pred_mat %*% ens$weights)
    }
    
    writeValues(out, pred_vals, start = row_start, nrows = nrows_i)
    
    rm(vals, complete_idx, pred_vals)
    gc()
  }
  
  readStop(raster_stack)
  out <- writeStop(out)
  
  out
}


# ── helper: stack rasters ────────────────────────────────────────────────────
get_stack <- function(varlist, raster_dir) {
  
  tif_paths <- file.path(raster_dir, paste0(varlist, ".tif"))
  
  if (!all(file.exists(tif_paths))) {
    missing_files <- tif_paths[!file.exists(tif_paths)]
    cat("  [WARN] missing tif:", paste(basename(missing_files), collapse = ", "), "\n")
    return(NULL)
  }
  
  raster_stack        <- rast(tif_paths)
  names(raster_stack) <- varlist
  raster_stack
}


# ── helper: dual suitability ─────────────────────────────────────────────────
calculate_dual_suitability <- function(climate, soil, threshold = 0.2) {
  ifel(!is.na(climate) & !is.na(soil) & soil > threshold, climate, 0)
}


# 1. Predict normal climate + soil + dual suitability ========================

for (i in zones) {
  
  cat("\n[RUN] zone", i, "\n")
  
  clim_file <- file.path(out_dir_clm,  paste0("clim_suit_zone", i, ".tif"))
  soil_file <- file.path(out_dir_soil, paste0("soil_suit_zone", i, ".tif"))
  dual_file <- file.path(out_dir_dual, paste0("dual_suitability_zone", i, ".tif"))
  
  if (file.exists(dual_file)) {
    cat("[SKIP] dual already exists: zone", i, "\n"); next
  }
  
  
  # 1.1 Climate suitability ---------------------------------------------------
  if (!file.exists(clim_file)) {
    
    mod_file <- file.path(clm_mod_dir, paste0("clm_ens_zone", i, ".Rdata"))
    if (!file.exists(mod_file)) {
      cat("[SKIP] climate model not found: zone", i, "\n"); next
    }
    
    load(mod_file)   # loads clim_ens
    rfVar <- clim_ens$varlist
    cat("  climate vars:", length(rfVar), "\n")
    
    raster_stack <- get_stack(rfVar, climdir_normal)
    if (is.null(raster_stack)) {
      cat("[SKIP] missing climate rasters: zone", i, "\n")
      rm(clim_ens); gc(); next
    }
    
    tmp_clim <- file.path(tmp_dir, paste0("tmp_clim_zone", i, ".tif"))
    
    pclim <- predict_ensemble_raster(
      clim_ens,
      raster_stack,
      filename    = tmp_clim,
      overwrite   = TRUE,
      chunk_nrows = 50L
    )
    
    pclim        <- resample(pclim, original_raster, method = "bilinear")
    names(pclim) <- "climate_suitability"
    
    writeRaster(pclim, clim_file, overwrite = TRUE)
    
    remove_tmp_files(tmp_clim)
    
    rm(clim_ens, raster_stack); gc()
    
  } else {
    cat("  [SKIP] climate suitability already exists\n")
    pclim <- rast(clim_file)
  }
  
  
  # 1.2 Soil suitability ------------------------------------------------------
  if (!file.exists(soil_file)) {
    
    mod_file <- file.path(soil_mod_dir, paste0("soil_ens_zone", i, ".Rdata"))
    if (!file.exists(mod_file)) {
      cat("[SKIP] soil model not found: zone", i, "\n"); next
    }
    
    load(mod_file)   # loads soil_ens
    rfVar <- soil_ens$varlist
    cat("  soil vars:", length(rfVar), "\n")
    
    raster_stack <- get_stack(rfVar, soilDir)
    if (is.null(raster_stack)) {
      cat("[SKIP] missing soil rasters: zone", i, "\n")
      rm(soil_ens); gc(); next
    }
    
    tmp_soil <- file.path(tmp_dir, paste0("tmp_soil_zone", i, ".tif"))
    
    psoil <- predict_ensemble_raster(
      soil_ens,
      raster_stack,
      filename    = tmp_soil,
      overwrite   = TRUE,
      chunk_nrows = 50L
    )
    
    psoil        <- crop(psoil, ext(original_raster))
    psoil        <- resample(psoil, original_raster, method = "bilinear")
    psoil        <- mask(psoil, original_raster)
    names(psoil) <- "soil_suitability"
    
    writeRaster(psoil, soil_file, overwrite = TRUE)
    
    remove_tmp_files(tmp_soil)
    
    rm(soil_ens, raster_stack); gc()
    
  } else {
    cat("  [SKIP] soil suitability already exists\n")
    psoil <- rast(soil_file)
  }
  
  
  # 1.3 Dual suitability ------------------------------------------------------
  if (!compareGeom(pclim, psoil, stopOnError = FALSE)) {
    psoil <- resample(psoil, pclim, method = "bilinear")
  }
  
  dual_suitability        <- calculate_dual_suitability(pclim, psoil, threshold = threshold)
  names(dual_suitability) <- "dual_suitability"
  
  writeRaster(dual_suitability, dual_file, overwrite = TRUE)
  
  cat("[DONE] zone", i, "\n")
  
  rm(pclim, psoil, dual_suitability); gc()
}



# 2. Predict future climate + overlay with normal soil suitability ============
# Future dual suitability = future climate suitability filtered by normal soil suitability

scenarios    <- c("ssp245", "ssp585")
time_periods <- c("2011-2040", "2041-2070", "2071-2100")

clim_mod_dir <- 'H:/Jing/ecoChina2/rf_final'
soil_dir     <- 'H:/Jing/ecoChina2/soil suitability/normal'
tmp_dir      <- 'H:/Jing/ecoChina2/tmp_prediction'

dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

for (scenario in scenarios) {
  for (period in time_periods) {
    
    cat("\n=== Future dual suitability:", toupper(scenario), period, "===\n")
    
    climdir_fut <- paste0(
      'H:/Jing/ecoChina/play/China/ClimateData/CN/800m/8GCMs_ensemble_',
      scenario, '_', period
    )
    
    out_dir_clm_fut <- paste0(
      'H:/Jing/ecoChina2/clim suitability/',
      period, toupper(scenario)
    )
    
    out_dir_dual_fut <- paste0(
      'H:/Jing/ecoChina2/dual suit/',
      period, toupper(scenario)
    )
    
    dir.create(out_dir_clm_fut,  recursive = TRUE, showWarnings = FALSE)
    dir.create(out_dir_dual_fut, recursive = TRUE, showWarnings = FALSE)
    
    
    for (i in zones) {
      
      cat("\n[RUN]", toupper(scenario), period, "| zone", i, "\n")
      
      clim_file <- file.path(out_dir_clm_fut,  paste0("clim_suit_zone", i, ".tif"))
      soil_file <- file.path(soil_dir,         paste0("soil_suit_zone", i, ".tif"))
      dual_file <- file.path(out_dir_dual_fut, paste0("dual_suitability_zone", i, ".tif"))
      
      if (file.exists(dual_file)) {
        cat("[SKIP] future dual already exists: zone", i, "\n"); next
      }
      
      
      # 2.1 Future climate suitability ---------------------------------------
      if (!file.exists(clim_file)) {
        
        mod_file <- file.path(clim_mod_dir, paste0("clm_ens_zone", i, ".Rdata"))
        if (!file.exists(mod_file)) {
          cat("[SKIP] climate model not found: zone", i, "\n"); next
        }
        
        load(mod_file)   # loads clim_ens
        rfVar <- clim_ens$varlist
        cat("  future climate vars:", length(rfVar), "\n")
        
        raster_stack <- get_stack(rfVar, climdir_fut)
        if (is.null(raster_stack)) {
          cat("[SKIP] missing future climate rasters: zone", i, "\n")
          rm(clim_ens); gc(); next
        }
        
        tmp_clim <- file.path(
          tmp_dir,
          paste0("tmp_clim_", scenario, "_", gsub("-", "", period), "_zone", i, ".tif")
        )
        
        pclim <- predict_ensemble_raster(
          clim_ens,
          raster_stack,
          filename    = tmp_clim,
          overwrite   = TRUE,
          chunk_nrows = 50L
        )
        
        pclim        <- resample(pclim, original_raster, method = "bilinear")
        names(pclim) <- "future_climate_suitability"
        
        writeRaster(pclim, clim_file, overwrite = TRUE)
        
        remove_tmp_files(tmp_clim)
        
        rm(clim_ens, raster_stack); gc()
        
      } else {
        cat("  [SKIP] future climate suitability already exists\n")
        pclim <- rast(clim_file)
      }
      
      
      # 2.2 Load normal soil suitability --------------------------------------
      if (!file.exists(soil_file)) {
        cat("[SKIP] normal soil suitability not found: zone", i, "\n")
        rm(pclim); gc(); next
      }
      
      psoil <- rast(soil_file)
      
      
      # 2.3 Future dual suitability ------------------------------------------
      if (!compareGeom(pclim, psoil, stopOnError = FALSE)) {
        psoil <- resample(psoil, pclim, method = "bilinear")
      }
      
      dual_suitability        <- calculate_dual_suitability(pclim, psoil, threshold = threshold)
      names(dual_suitability) <- "future_dual_suitability"
      
      writeRaster(dual_suitability, dual_file, overwrite = TRUE)
      
      cat("[DONE]", toupper(scenario), period, "| zone", i, "\n")
      
      rm(pclim, psoil, dual_suitability); gc()
    }
  }
}