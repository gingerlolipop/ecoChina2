# Zone 34 soil RF bug fix
# 1. Generate the missing soil variable-selection file and four soil RF models.
# 2. Predict zone 34 soil suitability for all four model versions.
# 3. Generate zone 34 dual suitability wherever climate suitability exists.
# ============================================================================

library(CEMT)
library(terra)
library(randomForest)
library(data.table)
library(foreach)
library(doSNOW)

rm(list = ls())
gc()

base_dir <- "H:/Jing/ecoChina2"
soil_dir <- "H:/Jing/soil rasters/tif2"

source(file.path(base_dir, "functions", "mcRFop_cls3.R"))
source(file.path(base_dir, "functions", "proportional sampling with print.R"))

# Run controls. After the models and soil rasters are created, training can be
# disabled when this file is reused only to repair later dual-suitability maps.
run_training <- TRUE
run_soil_prediction <- TRUE
run_dual_prediction <- TRUE

# All four soil RF files are trained. For the current repair, only optimized_mf
# needs soil and dual prediction because the other climate predictions are not ready.
methods_to_predict <- "optimized_mf"
# methods_to_predict <- c("optimized_mf", "optimized_rf", "plain_mf", "plain_rf")

i <- 34L
colname <- paste0("zone", i)
soil_threshold <- 0.2
predict_cores <- min(4L, max(1L, parallel::detectCores() - 1L))

VAR_ROW      <- 8L
NTREE_OPLIST <- 100L
NTREE_PLAIN  <- 500L
NTREE_MF     <- 100L
NFOREST      <- 10L
MF_NR        <- 1.3
NTREE1       <- 100L
NTREE2       <- 500L
NOP          <- 3L
THD          <- 0.75
BASE_SEED    <- 49L

OUT_DIR <- file.path(base_dir, "results")
MOD_DIR <- file.path(base_dir, "rf_soil")
TMP_DIR <- file.path(base_dir, "tmp_prediction")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MOD_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TMP_DIR, recursive = TRUE, showWarnings = FALSE)
terraOptions(tempdir = TMP_DIR, memfrac = 0.15)

# Multi-Forest functions ======================================================

mcmfRF2 <- function(xy_y, xy_n, nr = 1.2, varList, yCol,
                    reg = FALSE, nTree = 100, nForest = 10) {
  nCore <- min(max(1L, parallel::detectCores() - 1L), nTree)
  ntree_vec <- rep(floor(nTree / nCore), nCore)
  if (nTree %% nCore > 0) {
    ntree_vec[seq_len(nTree %% nCore)] <-
      ntree_vec[seq_len(nTree %% nCore)] + 1L
  }
  ntree_vec <- ntree_vec[ntree_vec > 0]
  
  cl <- makeCluster(length(ntree_vec), type = "SOCK")
  registerDoSNOW(cl)
  on.exit(stopCluster(cl), add = TRUE)
  
  n_prs <- nrow(xy_y)
  n_abs <- min(nrow(xy_n), floor(n_prs * nr))
  
  for (f in seq_len(nForest)) {
    train_abs <- xy_n[sample(seq_len(nrow(xy_n)), n_abs), ]
    train <- rbind(xy_y, train_abs)
    x2 <- train[, varList, drop = FALSE]
    y2 <- if (reg) train[[yCol]] else factor(train[[yCol]], levels = c(0, 1))
    
    rf2 <- foreach(
      ntree = ntree_vec,
      .combine = combine,
      .packages = "randomForest"
    ) %dopar% randomForest(x2, y2, ntree = ntree, importance = TRUE)
    
    if (f == 1L) rfC <- rf2 else rfC <- combine(rfC, rf2)
  }
  
  rfC
}

mcmfRFop <- function(xy_y, xy_n, nr = 1.2, varList, yCol,
                     nTree = 100, nForest = 10, nP = 10, thd = 0.8) {
  nCore <- min(max(1L, parallel::detectCores() - 1L), nTree)
  ntree_vec <- rep(floor(nTree / nCore), nCore)
  if (nTree %% nCore > 0) {
    ntree_vec[seq_len(nTree %% nCore)] <-
      ntree_vec[seq_len(nTree %% nCore)] + 1L
  }
  ntree_vec <- ntree_vec[ntree_vec > 0]
  
  cl <- makeCluster(length(ntree_vec), type = "SOCK")
  registerDoSNOW(cl)
  on.exit(stopCluster(cl), add = TRUE)
  
  n_prs <- nrow(xy_y)
  n_abs <- min(nrow(xy_n), floor(n_prs * nr))
  
  for (f in seq_len(nForest)) {
    train_abs <- xy_n[sample(seq_len(nrow(xy_n)), n_abs), ]
    train <- rbind(xy_y, train_abs)
    x2 <- train[, varList, drop = FALSE]
    y2 <- factor(train[[yCol]], levels = c(0, 1))
    
    Op <- classOP(x2, y2, nTree1 = 5, nTree2 = 10, nOP = nP, thd = thd)
    
    rf2 <- foreach(
      ntree = ntree_vec,
      .combine = combine,
      .packages = "randomForest"
    ) %dopar% randomForest(Op$x, Op$y, ntree = ntree, importance = TRUE)
    
    if (f == 1L) rfC <- rf2 else rfC <- combine(rfC, rf2)
  }
  
  rfC
}

# Raster helpers ==============================================================

remove_raster_files <- function(filepath) {
  base <- tools::file_path_sans_ext(filepath)
  files <- c(
    filepath,
    paste0(filepath, ".aux.xml"),
    paste0(filepath, ".ovr"),
    paste0(base, ".aux.xml"),
    paste0(base, ".ovr"),
    paste0(base, ".tfw")
  )
  files <- files[file.exists(files)]
  if (length(files) > 0) unlink(files, force = TRUE)
}

load_rf <- function(file, object_name) {
  e <- new.env()
  loaded <- load(file, envir = e)
  if (!(object_name %in% loaded)) {
    stop("Object not found in ", basename(file), ": ", object_name)
  }
  get(object_name, envir = e)
}

get_varlist <- function(model) {
  if (!is.null(model$varlist) && length(model$varlist) > 0) {
    return(model$varlist)
  }
  if (!is.null(model$importance) && nrow(model$importance) > 0) {
    return(rownames(model$importance))
  }
  stop("No predictor names found in model object")
}

rf_prob1 <- function(model, data) {
  predict(model, as.data.frame(data), type = "prob")[, "1"]
}

predict_soil <- function(model_file, object_name, out_file, template) {
  model <- load_rf(model_file, object_name)
  varlist <- get_varlist(model)
  tif <- file.path(soil_dir, paste0(varlist, ".tif"))
  
  missing <- tif[!file.exists(tif)]
  if (length(missing) > 0) {
    stop("Missing soil predictor rasters: ", paste(basename(missing), collapse = ", "))
  }
  
  s <- rast(tif)
  names(s) <- varlist
  same_geom <- compareGeom(s[[1]], template, stopOnError = FALSE)
  
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  remove_raster_files(out_file)
  
  if (same_geom) {
    p <- predict(
      s, model,
      fun = rf_prob1,
      filename = out_file,
      overwrite = TRUE,
      cores = predict_cores,
      cpkgs = "randomForest",
      na.rm = TRUE,
      wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW")
    )
  } else {
    tmp_file <- file.path(TMP_DIR, paste0("soil_zone34_", basename(out_file)))
    remove_raster_files(tmp_file)
    
    p0 <- predict(
      s, model,
      fun = rf_prob1,
      filename = tmp_file,
      overwrite = TRUE,
      cores = predict_cores,
      cpkgs = "randomForest",
      na.rm = TRUE,
      wopt = list(datatype = "FLT4S", gdal = "COMPRESS=NONE")
    )
    
    p <- resample(
      p0, template,
      method = "bilinear",
      filename = out_file,
      overwrite = TRUE,
      wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW")
    )
    
    remove_raster_files(tmp_file)
    rm(p0)
  }
  
  rm(model, s, p)
  gc()
  invisible(out_file)
}

calculate_dual_suitability <- function(climate, soil,
                                       threshold = soil_threshold) {
  ifel(is.na(soil), NA, ifel(soil > threshold, climate, 0))
}

# 1. Train all four zone-34 soil models =======================================

if (run_training) {
  cat("\n=== PREPARE ZONE 34 SOIL DATA ===\n")
  
  soil_dat <- as.data.table(fRead("H:/Jing/soil rasters/new_soil_raster.csv"))
  soil_dat <- soil_dat[, c(1:6, 8:22), with = FALSE]
  
  # Keep the same zoneID and 15 soil predictors used by the main soil script.
  soil_train <- soil_dat[, -c(1, 3:6), with = FALSE]
  xlist <- names(soil_train)[2:16]
  soil_train[, (colname) := as.integer(zoneID == i)]
  
  n_presence_raw <- sum(soil_train[[colname]] == 1, na.rm = TRUE)
  valid <- complete.cases(soil_train[, ..xlist])
  n_presence_valid <- sum(soil_train[[colname]][valid] == 1, na.rm = TRUE)
  
  cat("Zone 34 presence before soil-NA removal:", n_presence_raw, "\n")
  cat("Zone 34 presence after soil-NA removal:", n_presence_valid, "\n")
  
  if (n_presence_raw == 0L) {
    stop("Zone 34 is absent from new_soil_raster.csv")
  }
  if (n_presence_valid < 2L) {
    miss <- colSums(is.na(soil_train[soil_train[[colname]] == 1, ..xlist]))
    miss <- sort(miss[miss > 0], decreasing = TRUE)
    print(miss)
    stop("Too few complete zone-34 soil records. Review the missingness printed above.")
  }
  
  # proSysSmpl()/CEMT::meanFrq() uses legacy data.frame column selection.
  # Passing a data.table makes it search for a literal column named "byCols".
  soil_hot <- as.data.frame(soil_train[valid])
  
  # Variable selection: one proportional sample from the original complete dataset.
  # This matches the one-round workflow used by the other successfully trained zones.
  set.seed(BASE_SEED + i)
  soil_op <- proSysSmpl(
    soil_hot,
    byCol = which(names(soil_hot) == colname),
    minSz = 2000
  )
  soil_op <- as.data.frame(soil_op)
  
  # Balance presence/absence before variable selection.
  n_pos <- sum(soil_op[[colname]] == 1)
  n_neg <- sum(soil_op[[colname]] == 0)
  if (n_neg > n_pos * 1.5) {
    idx_pos <- which(soil_op[[colname]] == 1)
    idx_neg <- which(soil_op[[colname]] == 0)
    idx_neg <- sample(idx_neg, min(length(idx_neg), round(n_pos * 1.3)))
    soil_op <- soil_op[c(idx_pos, idx_neg), ]
  }
  print(table(soil_op[[colname]]))
  
  if (length(unique(soil_op[[colname]])) < 2L) {
    stop("Variable-selection sample contains only one class")
  }
  
  soilopList <- mcRFop_cls(
    soil_op[, xlist, drop = FALSE],
    factor(soil_op[[colname]], levels = c(0, 1)),
    nTree = NTREE_OPLIST
  )
  
  opfile <- file.path(OUT_DIR, paste0("soilopList_zone", i, ".csv"))
  write.csv(as.data.frame(soilopList), opfile, row.names = FALSE)
  cat("[SAVED]", opfile, "\n")
  
  if (nrow(soilopList) < VAR_ROW + 1L || soilopList[VAR_ROW + 1L, 2] == "0") {
    stop("No valid optimized variable set at VAR_ROW = ", VAR_ROW)
  }
  
  varlist <- trimws(unlist(strsplit(soilopList[VAR_ROW + 1L, 2], ",")))
  varlist <- intersect(varlist, xlist)
  if (length(varlist) < 2L) stop("Too few optimized soil predictors")
  
  cat("Optimized soil predictors:", length(varlist), "\n")
  print(varlist)
  
  # RF training: repeat the same one-round proportional sampling from the original data,
  # then cap absences to about 1.3 times presences, as in the main soil script.
  set.seed(BASE_SEED + i)
  soil_ps <- proSysSmpl(
    soil_hot,
    byCol = which(names(soil_hot) == colname),
    minSz = 2000
  )
  soil_ps <- as.data.frame(soil_ps)
  
  n_pos <- sum(soil_ps[[colname]] == 1)
  n_neg <- sum(soil_ps[[colname]] == 0)
  if (n_neg > n_pos * 1.5) {
    idx_pos <- which(soil_ps[[colname]] == 1)
    idx_neg <- which(soil_ps[[colname]] == 0)
    idx_neg <- sample(idx_neg, min(length(idx_neg), round(n_pos * 1.3)))
    soil_ps <- soil_ps[c(idx_pos, idx_neg), ]
  }
  print(table(soil_ps[[colname]]))
  
  xy_y <- soil_ps[soil_ps[[colname]] == 1, , drop = FALSE]
  xy_n <- soil_ps[soil_ps[[colname]] == 0, , drop = FALSE]
  soil_y <- factor(soil_ps[[colname]], levels = c(0, 1))
  soil_x0 <- soil_ps[, xlist, drop = FALSE]
  soil_x <- soil_ps[, varlist, drop = FALSE]
  
  cat("\n[TRAIN] plain RF\n")
  soil_plain <- randomForest(
    soil_x0, soil_y,
    ntree = NTREE_PLAIN,
    importance = TRUE
  )
  soil_plain$varlist <- xlist
  save(soil_plain, file = file.path(MOD_DIR, "soil_plain_zone34.Rdata"))
  
  cat("[TRAIN] plain multi-Forest RF\n")
  soil_mf <- mcmfRF2(
    xy_y, xy_n,
    nr = MF_NR,
    varList = xlist,
    yCol = colname,
    reg = FALSE,
    nTree = NTREE_MF,
    nForest = NFOREST
  )
  soil_mf$varlist <- xlist
  save(soil_mf, file = file.path(MOD_DIR, "soil_mf_zone34.Rdata"))
  
  cat("[TRAIN] optimized RF\n")
  soil_zOp <- classOP(
    soil_x, soil_y,
    nTree1 = NTREE1,
    nTree2 = NTREE2,
    nOP = NOP,
    thd = THD
  )
  soil_zOp$varlist <- varlist
  save(soil_zOp, file = file.path(MOD_DIR, "soil_zOp_zone34.Rdata"))
  
  cat("[TRAIN] optimized multi-Forest RF\n")
  soil_mfOp <- mcmfRFop(
    xy_y, xy_n,
    nr = MF_NR,
    varList = varlist,
    yCol = colname,
    nTree = NTREE_MF,
    nForest = NFOREST,
    nP = NOP,
    thd = THD
  )
  soil_mfOp$varlist <- varlist
  save(soil_mfOp, file = file.path(MOD_DIR, "soil_mfOp_zone34.Rdata"))
  
  cat("\n[SAVED] all four zone-34 soil RF files\n")
  print(file.info(file.path(
    MOD_DIR,
    c(
      "soil_plain_zone34.Rdata",
      "soil_mf_zone34.Rdata",
      "soil_zOp_zone34.Rdata",
      "soil_mfOp_zone34.Rdata"
    )
  ))[, c("size", "mtime")])
  
  rm(
    soil_dat, soil_train, soil_hot, soil_op, soilopList, soil_ps,
    xy_y, xy_n, soil_y, soil_x0, soil_x,
    soil_plain, soil_mf, soil_zOp, soil_mfOp
  )
  gc()
}

# 2. Predict zone-34 soil suitability for all model versions ==================

model_set <- data.table(
  method = c("optimized_mf", "optimized_rf", "plain_mf", "plain_rf"),
  soil_prefix = c("soil_mfOp_zone", "soil_zOp_zone", "soil_mf_zone", "soil_plain_zone"),
  soil_object = c("soil_mfOp", "soil_zOp", "soil_mf", "soil_plain")
)

model_set_pred <- model_set[method %in% methods_to_predict]

r <- rast(file.path(base_dir, "raster/ecosys_ori.tif"))

if (run_soil_prediction) {
  cat("\n=== PREDICT ZONE 34 SOIL SUITABILITY ===\n")
  
  for (j in seq_len(nrow(model_set_pred))) {
    method <- model_set_pred$method[j]
    model_file <- file.path(
      MOD_DIR,
      paste0(model_set_pred$soil_prefix[j], i, ".Rdata")
    )
    soil_file <- file.path(
      base_dir, "soil suitability", method, "normal",
      paste0("soil_suit_zone", i, ".tif")
    )
    
    cat("\n[PREDICT SOIL]", method, "| zone", i, "\n")
    predict_soil(
      model_file = model_file,
      object_name = model_set_pred$soil_object[j],
      out_file = soil_file,
      template = r
    )
    cat("[SAVED]", soil_file, "\n")
  }
}

# 3. Generate zone-34 dual suitability wherever climate rasters exist =========

if (run_dual_prediction) {
  cat("\n=== GENERATE ZONE 34 DUAL SUITABILITY ===\n")
  
  periods <- c(
    "normal",
    "2011-2040SSP245", "2041-2070SSP245", "2071-2100SSP245",
    "2011-2040SSP585", "2041-2070SSP585", "2071-2100SSP585"
  )
  
  for (method in model_set_pred$method) {
    soil_file <- file.path(
      base_dir, "soil suitability", method, "normal",
      paste0("soil_suit_zone", i, ".tif")
    )
    
    if (!file.exists(soil_file)) {
      cat("[SKIP] missing soil suitability:", method, "\n")
      next
    }
    
    for (period in periods) {
      clim_file <- file.path(
        base_dir, "clim suitability", method, period,
        paste0("clim_suit_zone", i, ".tif")
      )
      
      if (!file.exists(clim_file)) {
        cat("[SKIP] climate not available:", method, "|", period, "\n")
        next
      }
      
      dual_file <- file.path(
        base_dir, "dual suit", method, period,
        paste0("dual_suitability_zone", i, ".tif")
      )
      dir.create(dirname(dual_file), recursive = TRUE, showWarnings = FALSE)
      
      pclim <- rast(clim_file)
      psoil <- rast(soil_file)
      
      if (!compareGeom(pclim, psoil, stopOnError = FALSE)) {
        psoil <- resample(psoil, pclim, method = "bilinear")
      }
      
      if (compareGeom(psoil, r, stopOnError = FALSE)) {
        psoil <- mask(psoil, r)
      }
      
      dual <- calculate_dual_suitability(pclim, psoil)
      names(dual) <- if (period == "normal") {
        "dual_suitability"
      } else {
        "future_dual_suitability"
      }
      
      remove_raster_files(dual_file)
      writeRaster(
        dual,
        dual_file,
        overwrite = TRUE,
        wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW")
      )
      
      cat("[SAVED]", method, "|", period, "| zone", i, "\n")
      rm(pclim, psoil, dual)
      gc()
    }
  }
}

cat("\nZONE 34 BUG FIX COMPLETE\n")
