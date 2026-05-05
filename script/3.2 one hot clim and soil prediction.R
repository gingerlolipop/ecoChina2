#3.2 one hot clim and soil prediction
# Dual Suitability Prediction — Part 1: Climate suitability maps
# ============================================================
library(raster)
library(terra)
library(randomForest)
library(CEMT)

setwd("H:/Jing/ecoChina/play/China")

original_raster <- raster('raster/veg_3')
climdir_normal  <- 'H:/Jing/ecoChina/play/China/ClimateData/CN/800m/Normal_1961_1990'
clm_mod_dir     <- 'H:/Jing/ecoChina2/rf_final'
out_dir_clm     <- 'H:/Jing/ecoChina/play/China/clim suitability/normal'

dir.create(out_dir_clm, recursive = TRUE, showWarnings = FALSE)

zones <- setdiff(1:55, 8)


# ── helper: predict rfEnsemble onto raster stack ─────────────
predict_ensemble_raster <- function(ens, raster_stack, template) {
  vals         <- as.data.frame(values(raster_stack))
  complete_idx <- complete.cases(vals[, ens$varlist, drop = FALSE])
  pred_vals    <- rep(NA_real_, nrow(vals))
  
  if (sum(complete_idx) > 0) {
    X     <- vals[complete_idx, ens$varlist, drop = FALSE]
    plist <- lapply(ens$models,
                    function(m) predict(m, X, type = "prob")[, "1"])
    pred_vals[complete_idx] <- as.numeric(do.call(cbind, plist) %*% ens$weights)
  }
  
  r_out         <- raster(template)
  values(r_out) <- pred_vals
  r_out
}


# 1. Predict climate suitability — Normal 1961-1990 ==========================

for (i in zones) {
  
  out_file <- file.path(out_dir_clm, paste0("clim_suit_zone", i, ".tif"))
  if (file.exists(out_file)) {
    cat("[SKIP] already exists: zone", i, "\n"); next
  }
  
  setwd("H:/Jing/ecoChina/play/China")
  
  # Load climate ensemble model
  mod_file <- file.path(clm_mod_dir, paste0("clm_ens_zone", i, ".Rdata"))
  if (!file.exists(mod_file)) {
    cat("[SKIP] model not found: zone", i, "\n"); next
  }
  load(mod_file)   # loads clim_ens
  rfVar <- clim_ens$varlist
  cat("[RUN] zone", i, "| vars:", length(rfVar), "\n")
  
  # Stack climate rasters
  raster_layers <- list()
  for (var in rfVar) {
    tif_path <- file.path(climdir_normal, paste0(var, ".tif"))
    if (file.exists(tif_path)) {
      r        <- raster(tif_path)
      r[r == -9999] <- NA
      names(r) <- var
      raster_layers[[var]] <- r
    } else {
      cat("  [WARN] tif not found:", var, "\n")
    }
  }
  
  if (length(raster_layers) < length(rfVar)) {
    cat("[SKIP] missing raster files: zone", i, "\n")
    rm(clim_ens, raster_layers); gc(); next
  }
  
  raster_stack <- stack(raster_layers)
  
  # Predict
  pclim <- predict_ensemble_raster(clim_ens, raster_stack, raster_stack[[1]])
  pclim_resampled <- resample(pclim, original_raster)
  
  plot(pclim_resampled, main = paste0("Zone ", i, " Climate Suitability (Normal)"))
  
  # Save
  writeRaster(pclim_resampled, out_file, format = "GTiff", overwrite = TRUE)
  cat("[DONE] zone", i, "\n")
  
  rm(clim_ens, raster_layers, raster_stack, pclim, pclim_resampled); gc()
}


# 2. Predict climate suitability — SSP245 & SSP585, 3 periods ================

scenarios    <- c("ssp245", "ssp585")
time_periods <- c("2011-2040", "2041-2070", "2071-2100")

for (scenario in scenarios) {
  for (period in time_periods) {
    
    climdir_fut <- paste0('H:/Jing/ecoChina/play/China/ClimateData/CN/800m/8GCMs_ensemble_',
                          scenario, '_', period)
    out_dir_fut <- paste0('H:/Jing/ecoChina/play/China/clim suitability/',
                          period, toupper(scenario))
    dir.create(out_dir_fut, recursive = TRUE, showWarnings = FALSE)
    
    cat("\n=== Scenario:", toupper(scenario), "| Period:", period, "===\n")
    
    for (i in zones) {
      
      out_file <- file.path(out_dir_fut, paste0("clim_suit_zone", i, ".tif"))
      if (file.exists(out_file)) {
        cat("[SKIP] already exists: zone", i, "\n"); next
      }
      
      setwd("H:/Jing/ecoChina/play/China")
      
      mod_file <- file.path(clm_mod_dir, paste0("clm_ens_zone", i, ".Rdata"))
      if (!file.exists(mod_file)) {
        cat("[SKIP] model not found: zone", i, "\n"); next
      }
      load(mod_file)
      rfVar <- clim_ens$varlist
      cat("[RUN] zone", i, "| vars:", length(rfVar), "\n")
      
      raster_layers <- list()
      for (var in rfVar) {
        tif_path <- file.path(climdir_fut, paste0(var, ".tif"))
        if (file.exists(tif_path)) {
          r        <- raster(tif_path)
          r[r == -9999] <- NA
          names(r) <- var
          raster_layers[[var]] <- r
        } else {
          cat("  [WARN] tif not found:", var, "\n")
        }
      }
      
      if (length(raster_layers) < length(rfVar)) {
        cat("[SKIP] missing raster files: zone", i, "\n")
        rm(clim_ens, raster_layers); gc(); next
      }
      
      raster_stack <- stack(raster_layers)
      
      pclim <- predict_ensemble_raster(clim_ens, raster_stack, raster_stack[[1]])
      pclim_resampled <- resample(pclim, original_raster)
      
      writeRaster(pclim_resampled, out_file, format = "GTiff", overwrite = TRUE)
      cat("[DONE] zone", i, "\n")
      
      rm(clim_ens, raster_layers, raster_stack, pclim, pclim_resampled); gc()
    }
  }
}
# check
r <- 
  