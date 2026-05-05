# Dual Suitability Prediction — Part 1: Normal dual suitability
# ============================================================
library(terra)
library(randomForest)
library(CEMT)
library(data.table)

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

zones <- setdiff(1:55, 8)
threshold <- 0.2


# ── helper: predict rfEnsemble onto SpatRaster stack by blocks ─────────────
predict_ensemble_raster <- function(ens, raster_stack, filename, overwrite = TRUE) {
  
  out <- raster_stack[[1]]
  names(out) <- "suitability"
  
  bs <- terra::blocks(raster_stack)
  
  out <- terra::writeStart(out, filename = filename, overwrite = overwrite)
  
  for (b in 1:nrow(bs)) {
    
    vals <- terra::readValues(
      raster_stack,
      row   = bs$row[b],
      nrows = bs$nrows[b],
      mat   = TRUE
    )
    
    vals <- as.data.frame(vals)
    names(vals) <- names(raster_stack)
    vals <- vals[, ens$varlist, drop = FALSE]
    
    complete_idx <- complete.cases(vals)
    pred_vals <- rep(NA_real_, nrow(vals))
    
    if (sum(complete_idx) > 0) {
      
      X <- vals[complete_idx, , drop = FALSE]
      
      pred_mat <- do.call(cbind, lapply(
        ens$models,
        function(m) predict(m, X, type = "prob")[, "1"]
      ))
      
      pred_vals[complete_idx] <- as.numeric(pred_mat %*% ens$weights)
    }
    
    terra::writeValues(out, pred_vals, bs$row[b])
    
    rm(vals, complete_idx, pred_vals)
    gc()
  }
  
  out <- terra::writeStop(out)
  out
}


# ── helper: stack rasters with cache ─────────────────────────────
get_stack_cached <- function(varlist, raster_dir, cache_env) {
  
  raster_layers <- list()
  
  for (var in varlist) {
    
    if (!exists(var, envir = cache_env, inherits = FALSE)) {
      tif_path <- file.path(raster_dir, paste0(var, ".tif"))
      
      if (file.exists(tif_path)) {
        r <- rast(tif_path)
        r[r == -9999] <- NA
        names(r) <- var
        assign(var, r, envir = cache_env)
      } else {
        cat("  [WARN] tif not found:", var, "\n")
      }
    }
    
    if (exists(var, envir = cache_env, inherits = FALSE)) {
      raster_layers[[var]] <- get(var, envir = cache_env)
    }
  }
  
  if (length(raster_layers) < length(varlist)) return(NULL)
  
  raster_stack <- raster_layers[[1]]
  if (length(raster_layers) > 1) {
    for (j in 2:length(raster_layers)) {
      raster_stack <- c(raster_stack, raster_layers[[j]])
    }
  }
  
  names(raster_stack) <- varlist
  
  raster_stack
}


# ── helper: dual suitability function ───────────────────────────
calculate_dual_suitability <- function(climate, soil, threshold = 0.2) {
  ifel(!is.na(climate) & !is.na(soil) & soil > threshold, climate, 0)
}


# 1. Predict normal climate + soil + dual suitability ========================

clim_cache <- new.env(parent = emptyenv())
soil_cache <- new.env(parent = emptyenv())

for (i in zones) {
  
  cat("\n[RUN] zone", i, "\n")
  
  clim_file <- file.path(out_dir_clm,  paste0("clim_suit_zone", i, ".tif"))
  soil_file <- file.path(out_dir_soil, paste0("soil_suit_zone", i, ".tif"))
  dual_file <- file.path(out_dir_dual, paste0("dual_suitability_zone", i, ".rds"))
  
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
    
    raster_stack <- get_stack_cached(rfVar, climdir_normal, clim_cache)
    if (is.null(raster_stack)) {
      cat("[SKIP] missing climate rasters: zone", i, "\n")
      rm(clim_ens); gc(); next
    }
    
    tmp_clim <- file.path(tmp_dir, paste0("tmp_clim_zone", i, ".tif"))
    
    pclim <- predict_ensemble_raster(
      clim_ens,
      raster_stack,
      filename  = tmp_clim,
      overwrite = TRUE
    )
    
    pclim <- resample(pclim, original_raster, method = "bilinear")
    names(pclim) <- "climate_suitability"
    
    writeRaster(pclim, clim_file, overwrite = TRUE)
    
    if (file.exists(tmp_clim)) file.remove(tmp_clim)
    
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
    
    raster_stack <- get_stack_cached(rfVar, soilDir, soil_cache)
    if (is.null(raster_stack)) {
      cat("[SKIP] missing soil rasters: zone", i, "\n")
      rm(soil_ens); gc(); next
    }
    
    tmp_soil <- file.path(tmp_dir, paste0("tmp_soil_zone", i, ".tif"))
    
    psoil <- predict_ensemble_raster(
      soil_ens,
      raster_stack,
      filename  = tmp_soil,
      overwrite = TRUE
    )
    
    psoil <- crop(psoil, ext(original_raster))
    psoil <- resample(psoil, original_raster, method = "bilinear")
    psoil <- mask(psoil, original_raster)
    names(psoil) <- "soil_suitability"
    
    writeRaster(psoil, soil_file, overwrite = TRUE)
    
    if (file.exists(tmp_soil)) file.remove(tmp_soil)
    
    rm(soil_ens, raster_stack); gc()
    
  } else {
    cat("  [SKIP] soil suitability already exists\n")
    psoil <- rast(soil_file)
  }
  
  
  # 1.3 Dual suitability ------------------------------------------------------
  if (!compareGeom(pclim, psoil, stopOnError = FALSE)) {
    psoil <- resample(psoil, pclim, method = "bilinear")
  }
  
  dual_suitability <- calculate_dual_suitability(pclim, psoil, threshold = threshold)
  names(dual_suitability) <- "dual_suitability"
  
  saveRDS(dual_suitability, dual_file)
  
  cat("[DONE] zone", i, "\n")
  
  rm(pclim, psoil, dual_suitability); gc()
}