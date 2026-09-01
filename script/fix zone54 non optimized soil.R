library(CEMT)
library(terra)
library(randomForest)
library(data.table)
library(foreach)
library(doSNOW)

base_dir <- "H:/Jing/ecoChina2"
soilDir <- "H:/Jing/soil rasters/tif2"

source(file.path(base_dir, "functions", "mcRFop_cls3.R"))

i <- 54
ycol <- paste0("zone", i)

VAR_ROW <- 8L
NTREE_OPLIST <- 100L
NTREE_PLAIN <- 500L
NTREE_MF <- 100L
NFOREST <- 10L
MF_NR <- 1.3
BASE_SEED <- 49L
soil_threshold <- 0.2

result_dir <- file.path(base_dir, "results")
model_dir <- file.path(base_dir, "rf_soil")
test_dir <- file.path(base_dir, "testing data")
tmp_dir <- file.path(base_dir, "tmp_prediction")

dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

terraOptions(tempdir = tmp_dir, memfrac = 0.15)

r <- rast(file.path(base_dir, "raster/ecosys_ori.tif"))

# Use the original 70/30 split.
train <- fRead(file.path(result_dir, "soil_train_data.csv"))
test <- fRead(file.path(result_dir, "soil_test_data.csv"))

xlist <- names(train)[2:16]

train[[ycol]] <- as.integer(train$zoneID == i)
test[[ycol]] <- as.integer(test$zoneID == i)

train54 <- train[, c("zoneID", ycol, xlist), drop = FALSE]
test54 <- test[, c("zoneID", ycol, xlist), drop = FALSE]

fWrite(
  train54,
  file.path(test_dir, paste0("soil_train_zone", i, ".csv"))
)

fWrite(
  test54,
  file.path(test_dir, paste0("soil_test_zone", i, ".csv"))
)

cat("[TRAIN]\n")
print(table(train54[[ycol]]))

cat("[TEST]\n")
print(table(test54[[ycol]]))

if (sum(train54[[ycol]] == 1) < 2L) {
  stop("Too few zone-54 presence observations in the training set.")
}

if (sum(test54[[ycol]] == 1) < 1L) {
  warning("The original 30% testing set contains no zone-54 presence.")
}

# Keep all presences from the 70% training set and sample absences.
sample_pa54 <- function(dat, vars, seed) {
  dat <- as.data.frame(dat)
  
  dat <- dat[
    complete.cases(dat[, c(ycol, vars), drop = FALSE]),
    ,
    drop = FALSE
  ]
  
  pos <- which(dat[[ycol]] == 1)
  neg <- which(dat[[ycol]] == 0)
  
  if (length(pos) < 2L || length(neg) < 2L) {
    stop(
      "Too few training observations: presence = ",
      length(pos), ", absence = ", length(neg)
    )
  }
  
  set.seed(seed)
  
  n_neg <- min(
    length(neg),
    max(2L, floor(length(pos) * MF_NR))
  )
  
  neg <- sample(neg, n_neg, replace = FALSE)
  
  dat[c(pos, neg), , drop = FALSE]
}

mcmfRF2 <- function(xy_y, xy_n, nr = 1.2, varList, yCol,
                    nTree = 100, nForest = 10) {
  nCore <- min(max(1L, parallel::detectCores() - 1L), nTree)
  ntree_vec <- rep(floor(nTree / nCore), nCore)
  
  if (nTree %% nCore > 0L) {
    ntree_vec[seq_len(nTree %% nCore)] <-
      ntree_vec[seq_len(nTree %% nCore)] + 1L
  }
  
  ntree_vec <- ntree_vec[ntree_vec > 0L]
  
  cl <- parallel::makeCluster(length(ntree_vec), type = "SOCK")
  registerDoSNOW(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  n_abs <- min(nrow(xy_n), floor(nrow(xy_y) * nr))
  
  for (f in seq_len(nForest)) {
    train_abs <- xy_n[
      sample(seq_len(nrow(xy_n)), n_abs, replace = FALSE),
      ,
      drop = FALSE
    ]
    
    train_f <- rbind(xy_y, train_abs)
    x2 <- train_f[, varList, drop = FALSE]
    y2 <- factor(train_f[[yCol]], levels = c(0, 1))
    
    rf2 <- foreach(
      ntree = ntree_vec,
      .combine = combine,
      .packages = "randomForest"
    ) %dopar% {
      randomForest(
        x2,
        y2,
        ntree = ntree,
        importance = TRUE
      )
    }
    
    if (f == 1L) rfC <- rf2 else rfC <- combine(rfC, rf2)
  }
  
  rfC
}

# Variable selection
opfile <- file.path(
  result_dir,
  paste0("soilopList_zone", i, ".csv")
)

if (!file.exists(opfile)) {
  soil_ps <- sample_pa54(
    train54,
    xlist,
    BASE_SEED + i
  )
  
  soil_x <- soil_ps[, xlist, drop = FALSE]
  soil_y <- factor(soil_ps[[ycol]], levels = c(0, 1))
  
  cat("[OPLIST SAMPLE]\n")
  print(table(soil_y))
  
  soilopList <- mcRFop_cls(
    soil_x,
    soil_y,
    nTree = NTREE_OPLIST
  )
  
  write.csv(
    as.data.frame(soilopList),
    opfile,
    row.names = FALSE
  )
  
  cat("[OPLIST SAVED]", opfile, "\n")
}

# Read selected variables
soilList <- read.csv(opfile)

if (
  nrow(soilList) < VAR_ROW + 1L ||
  is.na(soilList[VAR_ROW + 1L, 2]) ||
  soilList[VAR_ROW + 1L, 2] == "0"
) {
  stop("No valid variable set for zone 54.")
}

varlist <- trimws(
  unlist(strsplit(soilList[VAR_ROW + 1L, 2], ","))
)

varlist <- intersect(varlist, xlist)

if (length(varlist) < 2L) {
  stop("Too few selected variables for zone 54.")
}

# Formal training data
soil_ps <- sample_pa54(
  train54,
  varlist,
  BASE_SEED + i
)

soil_x <- soil_ps[, varlist, drop = FALSE]
soil_y <- factor(soil_ps[[ycol]], levels = c(0, 1))

xy_y <- soil_ps[
  soil_ps[[ycol]] == 1,
  ,
  drop = FALSE
]

xy_n <- soil_ps[
  soil_ps[[ycol]] == 0,
  ,
  drop = FALSE
]

cat("[RF SAMPLE]\n")
print(table(soil_y))
cat("Predictors:", paste(varlist, collapse = ", "), "\n")

# Plain RF
soil_plain <- randomForest(
  soil_x,
  soil_y,
  ntree = NTREE_PLAIN,
  importance = TRUE
)

soil_plain$varlist <- varlist

save(
  soil_plain,
  file = file.path(
    model_dir,
    paste0("soil_plain_zone", i, ".Rdata")
  )
)

# Plain multi-Forest
soil_mf <- mcmfRF2(
  xy_y,
  xy_n,
  nr = MF_NR,
  varList = varlist,
  yCol = ycol,
  nTree = NTREE_MF,
  nForest = NFOREST
)

soil_mf$varlist <- varlist

save(
  soil_mf,
  file = file.path(
    model_dir,
    paste0("soil_mf_zone", i, ".Rdata")
  )
)

cat("[MODELS SAVED] zone", i, "\n")

rf_prob1 <- function(model, data) {
  predict(
    model,
    as.data.frame(data),
    type = "prob"
  )[, "1"]
}

predict_soil <- function(model, out_file) {
  files <- file.path(
    soilDir,
    paste0(model$varlist, ".tif")
  )
  
  if (!all(file.exists(files))) {
    stop(
      "Missing soil rasters: ",
      paste(basename(files[!file.exists(files)]), collapse = ", ")
    )
  }
  
  s <- rast(files)
  names(s) <- model$varlist
  
  cores <- min(
    4L,
    max(1L, parallel::detectCores() - 1L)
  )
  
  if (compareGeom(s[[1]], r, stopOnError = FALSE)) {
    predict(
      s,
      model,
      fun = rf_prob1,
      filename = out_file,
      overwrite = TRUE,
      cores = cores,
      cpkgs = "randomForest",
      na.rm = TRUE,
      wopt = list(
        datatype = "FLT4S",
        gdal = "COMPRESS=LZW"
      )
    )
  } else {
    tmp_file <- tempfile(
      pattern = "soil54_",
      tmpdir = tmp_dir,
      fileext = ".tif"
    )
    
    p0 <- predict(
      s,
      model,
      fun = rf_prob1,
      filename = tmp_file,
      overwrite = TRUE,
      cores = cores,
      cpkgs = "randomForest",
      na.rm = TRUE,
      wopt = list(datatype = "FLT4S")
    )
    
    resample(
      p0,
      r,
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
    unlink(tmp_file, force = TRUE)
  }
  
  rast(out_file)
}

models <- list(
  plain_mf = soil_mf,
  plain_rf = soil_plain
)

# Normal soil and dual suitability
for (method in names(models)) {
  soil_out_dir <- file.path(
    base_dir,
    "soil suitability",
    method,
    "normal"
  )
  
  dual_out_dir <- file.path(
    base_dir,
    "dual suit",
    method,
    "normal"
  )
  
  dir.create(soil_out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dual_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  soil_file <- file.path(
    soil_out_dir,
    paste0("soil_suit_zone", i, ".tif")
  )
  
  clim_file <- file.path(
    base_dir,
    "clim suitability",
    method,
    "normal",
    paste0("clim_suit_zone", i, ".tif")
  )
  
  dual_file <- file.path(
    dual_out_dir,
    paste0("dual_suitability_zone", i, ".tif")
  )
  
  if (!file.exists(clim_file)) {
    stop("Missing climate suitability: ", clim_file)
  }
  
  cat("\n[PREDICT]", method, "| zone", i, "\n")
  
  psoil <- predict_soil(models[[method]], soil_file)
  psoil <- mask(psoil, r)
  pclim <- rast(clim_file)
  
  if (!compareGeom(psoil, pclim, stopOnError = FALSE)) {
    psoil <- resample(psoil, pclim, method = "bilinear")
  }
  
  dual <- ifel(
    is.na(psoil),
    NA,
    ifel(psoil > soil_threshold, pclim, 0)
  )
  
  names(dual) <- "dual_suitability"
  
  writeRaster(
    dual,
    dual_file,
    overwrite = TRUE,
    wopt = list(
      datatype = "FLT4S",
      gdal = "COMPRESS=LZW"
    )
  )
  
  cat("[DONE]", method, "normal | zone", i, "\n")
  
  rm(psoil, pclim, dual)
  gc()
}


