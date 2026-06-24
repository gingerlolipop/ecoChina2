# Soil One-Hot RF — no spatial CV + multi-Forest
# Zone 1 first, then loop over zones 2:55 except zone 8
# ============================================================
library(CEMT)
library(data.table)
library(foreach)
library(doSNOW)

rm(list = ls()); gc()

base_dir <- "H:/Jing/ecoChina2"

source(file.path(base_dir, "functions", "mcRFop_cls3.R"))
source(file.path(base_dir, "functions", "proportional sampling with print.R"))

# 0. Multi-Forest functions ===================================================

mcmfRF2 <- function(xy_y, xy_n, nr = 1.2, varList, yCol,
                    reg = FALSE, nTree = 100, nForest = 10) {
  library(foreach); library(doSNOW); library(randomForest)
  nCore <- min(max(1L, parallel::detectCores() - 1L), nTree)
  ntree_vec <- rep(floor(nTree / nCore), nCore)
  if (nTree %% nCore > 0) {
    ntree_vec[seq_len(nTree %% nCore)] <- ntree_vec[seq_len(nTree %% nCore)] + 1L
  }
  ntree_vec <- ntree_vec[ntree_vec > 0]
  cl <- makeCluster(length(ntree_vec), type = "SOCK")
  registerDoSNOW(cl)
  on.exit(stopCluster(cl), add = TRUE)
  
  n_prs <- nrow(xy_y)
  n_abs <- min(nrow(xy_n), floor(n_prs * nr))
  
  for (f in 1:nForest) {
    train_abs <- xy_n[sample(seq_len(nrow(xy_n)), n_abs, replace = FALSE), ]
    train <- rbind(xy_y, train_abs)
    x2 <- train[, varList, drop = FALSE]
    if (!reg) y2 <- factor(train[[yCol]], levels = c(0, 1))
    if (reg)  y2 <- train[[yCol]]
    
    rf2 <- foreach(
      ntree = ntree_vec,
      .combine = combine,
      .packages = "randomForest"
    ) %dopar% randomForest(x2, y2, ntree = ntree, importance = TRUE)
    
    if (f == 1) rfC <- rf2 else rfC <- combine(rfC, rf2)
  }
  rfC
}

mcmfRFop <- function(xy_y, xy_n, nr = 1.2, varList, yCol,
                     nTree = 100, nForest = 10, nP = 10, thd = 0.8) {
  library(foreach); library(doSNOW); library(randomForest)
  nCore <- min(max(1L, parallel::detectCores() - 1L), nTree)
  ntree_vec <- rep(floor(nTree / nCore), nCore)
  if (nTree %% nCore > 0) {
    ntree_vec[seq_len(nTree %% nCore)] <- ntree_vec[seq_len(nTree %% nCore)] + 1L
  }
  ntree_vec <- ntree_vec[ntree_vec > 0]
  cl <- makeCluster(length(ntree_vec), type = "SOCK")
  registerDoSNOW(cl)
  on.exit(stopCluster(cl), add = TRUE)
  
  n_prs <- nrow(xy_y)
  n_abs <- min(nrow(xy_n), floor(n_prs * nr))
  
  for (f in 1:nForest) {
    train_abs <- xy_n[sample(seq_len(nrow(xy_n)), n_abs, replace = FALSE), ]
    train <- rbind(xy_y, train_abs)
    x2 <- train[, varList, drop = FALSE]
    y2 <- factor(train[[yCol]], levels = c(0, 1))
    
    Op <- classOP(x2, y2, nTree1 = 5, nTree2 = 10, nOP = nP, thd = thd)
    x3 <- Op$x
    y3 <- Op$y
    
    rf2 <- foreach(
      ntree = ntree_vec,
      .combine = combine,
      .packages = "randomForest"
    ) %dopar% randomForest(x3, y3, ntree = ntree, importance = TRUE)
    
    if (f == 1) rfC <- rf2 else rfC <- combine(rfC, rf2)
  }
  rfC
}

rf_acc <- function(m, x, y, zone, model_type, out_dir) {
  y <- factor(y, levels = c(0, 1))
  x <- as.data.frame(x)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Training-set confusion matrix. This is always calculated directly.
  pred <- factor(predict(m, x, type = "response"), levels = c(0, 1))
  cm_train <- table(observed = y, predicted = pred)
  train_accuracy <- sum(diag(cm_train)) / sum(cm_train)
  
  fwrite(
    as.data.table(cm_train),
    file.path(out_dir, paste0(model_type, "_confusion_train_zone", zone, ".csv"))
  )
  
  # OOB confusion matrix, if the randomForest object still has it.
  # Combined multi-Forest objects may not always retain reliable OOB statistics.
  oob_accuracy <- NA_real_
  if (!is.null(m$confusion)) {
    cm_oob_raw <- as.matrix(m$confusion)
    cm_oob <- cm_oob_raw[, setdiff(colnames(cm_oob_raw), "class.error"), drop = FALSE]
    oob_accuracy <- sum(diag(cm_oob)) / sum(cm_oob)
    
    fwrite(
      as.data.table(cm_oob, keep.rownames = "observed"),
      file.path(out_dir, paste0(model_type, "_confusion_oob_zone", zone, ".csv"))
    )
  }
  
  data.table(
    zone = zone,
    model = model_type,
    n = length(y),
    n_presence = sum(y == 1),
    n_absence = sum(y == 0),
    n_predictors = ncol(x),
    train_accuracy = train_accuracy,
    oob_accuracy = oob_accuracy
  )
}

# ── params ─────────────────────────────────────
VAR_ROW       <- 8L
NTREE_OPLIST  <- 100L
NTREE_PLAIN   <- 500L
NTREE_MF      <- 100L
NFOREST       <- 10L
MF_NR         <- 1.3
NTREE1        <- 100L
NTREE2        <- 500L
NOP           <- 3L
THD           <- 0.75
BASE_SEED     <- 49L
OUT_DIR       <- file.path(base_dir, "results")
MOD_DIR       <- file.path(base_dir, "rf_soil")
ACC_DIR       <- file.path(base_dir, "accuracy_soil")
# ───────────────────────────────────────────────

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MOD_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(ACC_DIR, showWarnings = FALSE, recursive = TRUE)

# 1. Import & prep ============================================================

setwd("H:/Jing/soil rasters")
soil_dat2 <- fRead("new_soil_raster.csv"); hd(soil_dat2)
soil_dat2 <- soil_dat2[, c(1:6, 8:22)]
names(soil_dat2)

soil_train <- soil_dat2[complete.cases(soil_dat2[, c(2, 7:ncol(soil_dat2))]), ]
cat("Rows after NA removal:", nrow(soil_train), "\n")
cat("Zones present:", sort(unique(soil_train$zoneID)), "\n")

soil_coords <- soil_train[, c("x", "y")]
soil_train2 <- soil_train[, -c(1, 3:6)]


# randomly split agg_data into training/testing 7:3
set.seed(49)
train_indices <- sample(nrow(soil_train2), size = 0.7 * nrow(soil_train2))
soil_train_data <- soil_train2[train_indices, ]
soil_coords <- soil_coords[train_indices, ]
soil_test_data <- soil_train2[-train_indices, ]

fWrite(soil_train_data, file.path(OUT_DIR, "soil_train_data.csv"))
fWrite(soil_coords, file.path(OUT_DIR, "soil_train_coords.csv"))

fWrite(soil_test_data, file.path(OUT_DIR, "soil_test_data.csv"))

zone_counts <- table(soil_train_data$zoneID); print(zone_counts)

colnames(soil_train_data)
xlist <- names(soil_train_data)[2:16]
cat("Soil predictors:", xlist, "\n")

rm(soil_train_data,soil_train,soil_train2,soil_dat2,soil_test_data);gc()

# 2. One-hot encoding =========================================================
library(caret)
soil_train_data <- fRead(file.path(OUT_DIR, "soil_train_data.csv"))


soil_train_data$zone <- as.character(soil_train_data$zoneID)
soil_hot_mat <- predict(dummyVars(~zone, data = soil_train_data), newdata = soil_train_data)
soil_hot <- cbind(soil_train_data, soil_hot_mat, soil_coords)

zones_all <- sort(setdiff(unique(soil_train_data$zoneID), 8))
cat("Total zones to process:", length(zones_all), "\n")

fWrite(soil_hot, file.path(OUT_DIR, "soil_hot_encoded.csv"))

# 3. Variable selection: zone 1 first =========================================
rm();gc()
library(randomForest)


soil_hot <-fRead(file.path(OUT_DIR, "soil_hot_encoded.csv"))

i <- 1
colname <- paste0("zone", i)
opfile <- file.path(OUT_DIR, paste0("soilopList_zone", i, ".csv"))

n1 <- sum(soil_hot[[colname]] == 1, na.rm = TRUE)
n_rounds <- if (n1 >= 6000) 3L else if (n1 >= 3000) 2L else 1L
cat("[OPLIST] zone", i, "| n_presence =", n1, "| sampling rounds =", n_rounds, "\n")

if (!file.exists(opfile)) {
  set.seed(BASE_SEED + i)
  soil_ps <- soil_hot
  for (rr in seq_len(n_rounds)) {
    soil_ps <- proSysSmpl(soil_ps,
                          byCol = which(colnames(soil_ps) == colname),
                          minSz = 2000)
  }
  print(table(soil_ps[[colname]]))
  
  soil_y <- factor(soil_ps[[colname]], levels = c(0, 1))
  soil_x <- soil_ps[, xlist]
  
  soilopList <- mcRFop_cls(soil_x, soil_y, nTree = NTREE_OPLIST)
  write.csv(as.data.frame(soilopList), opfile, row.names = FALSE)
  
  rm(soil_ps, soil_y, soil_x, soilopList); gc()
} else {
  cat("[SKIP] opList already exists: zone", i, "\n")
}

# 4. Variable selection: zones 2:55, skipping zone 8 ==========================

for (i in setdiff(2:55, 8)) {
  colname <- paste0("zone", i)
  
  if (!(colname %in% names(soil_hot))) {
    cat("[SKIP] no col:", colname, "\n"); next
  }
  
  n1 <- sum(soil_hot[[colname]] == 1, na.rm = TRUE)
  if (is.na(n1) || n1 == 0) {
    cat("[SKIP] no presence:", colname, "\n"); next
  }
  
  opfile <- file.path(OUT_DIR, paste0("soilopList_zone", i, ".csv"))
  if (file.exists(opfile)) {
    cat("[SKIP] opList already exists:", colname, "\n"); next
  }
  
  tryCatch({
    n_rounds <- if (n1 >= 6000) 3L else if (n1 >= 3000) 2L else 1L
    cat("[OPLIST] zone", i, "| n_presence =", n1, "| sampling rounds =", n_rounds, "\n")
    
    set.seed(BASE_SEED + i)
    soil_ps <- soil_hot
    for (rr in seq_len(n_rounds)) {
      soil_ps <- proSysSmpl(soil_ps,
                            byCol = which(colnames(soil_ps) == colname),
                            minSz = 2000)
    }
    print(table(soil_ps[[colname]]))
    
    if (length(unique(soil_ps[[colname]])) < 2) {
      cat("[SKIP] single class:", colname, "\n"); next
    }
    
    soil_y <- factor(soil_ps[[colname]], levels = c(0, 1))
    soil_x <- soil_ps[, xlist]
    
    soilopList <- mcRFop_cls(soil_x, soil_y, nTree = NTREE_OPLIST)
    write.csv(as.data.frame(soilopList), opfile, row.names = FALSE)
    
    rm(soil_ps, soil_y, soil_x, soilopList); gc()
    
  }, error = function(e) {
    cat("[ERROR] opList zone", i, ":", conditionMessage(e), "\n"); gc()
  })
}

# 5. RF training: zone 1 first ================================================

acc_all <- list()

i <- 1
colname <- paste0("zone", i)
opfile <- file.path(OUT_DIR, paste0("soilopList_zone", i, ".csv"))

n1 <- sum(soil_hot[[colname]] == 1, na.rm = TRUE)
n_rounds <- if (n1 >= 6000) 3L else if (n1 >= 3000) 2L else 1L
cat("[RF] zone", i, "| n_presence =", n1, "| sampling rounds =", n_rounds, "\n")

set.seed(BASE_SEED + i)
soil_ps <- soil_hot
for (rr in seq_len(n_rounds)) {
  soil_ps <- proSysSmpl(soil_ps,
                        byCol = which(colnames(soil_ps) == colname),
                        minSz = 2000)
}
soilList <- read.csv(opfile)
varlist <- trimws(unlist(strsplit(soilList[VAR_ROW + 1, 2], ",")))
varlist <- intersect(varlist, names(soil_ps))
soil_ps <- soil_ps[complete.cases(soil_ps[, c(colname, varlist)]), ]

n_pos <- sum(soil_ps[[colname]] == 1)
n_neg <- sum(soil_ps[[colname]] == 0)
if (n_neg > n_pos * 1.5) {
  idx_pos <- which(soil_ps[[colname]] == 1)
  idx_neg <- which(soil_ps[[colname]] == 0)
  idx_neg_sub <- sample(idx_neg, min(length(idx_neg), round(n_pos * 1.3)))
  soil_ps <- soil_ps[c(idx_pos, idx_neg_sub), ]
}
print(table(soil_ps[[colname]]))

soil_ps <- as.data.frame(soil_ps)
xy_y <- soil_ps[soil_ps[[colname]] == 1, ]
xy_n <- soil_ps[soil_ps[[colname]] == 0, ]
soil_y <- factor(soil_ps[[colname]], levels = c(0, 1))
soil_x <- soil_ps[, varlist, drop = FALSE]

# 5.1 Plain single RF
soil_plain <- randomForest(soil_x, soil_y, ntree = NTREE_PLAIN, importance = TRUE)
soil_plain$varlist <- varlist
acc_all[[length(acc_all) + 1L]] <- rf_acc(soil_plain, soil_x, soil_y, i, "plain_rf", ACC_DIR)
save(soil_plain, file = file.path(MOD_DIR, paste0("soil_plain_zone", i, ".Rdata")))

# 5.2 Plain multi-Forest RF
soil_mf <- mcmfRF2(xy_y, xy_n, nr = MF_NR, varList = varlist, yCol = colname,
                   reg = FALSE, nTree = NTREE_MF, nForest = NFOREST)
soil_mf$varlist <- varlist
acc_all[[length(acc_all) + 1L]] <- rf_acc(soil_mf, soil_x, soil_y, i, "plain_mf_rf", ACC_DIR)
save(soil_mf, file = file.path(MOD_DIR, paste0("soil_mf_zone", i, ".Rdata")))

# 5.3 Optimized single RF
soil_zOp <- classOP(soil_x, soil_y, nTree1 = NTREE1, nTree2 = NTREE2, nOP = NOP, thd = THD)
soil_zOp$varlist <- varlist
acc_all[[length(acc_all) + 1L]] <- rf_acc(soil_zOp, soil_x, soil_y, i, "optimized_rf", ACC_DIR)
save(soil_zOp, file = file.path(MOD_DIR, paste0("soil_zOp_zone", i, ".Rdata")))

# 5.4 Optimized multi-Forest RF
soil_mfOp <- mcmfRFop(xy_y, xy_n, nr = MF_NR, varList = varlist, yCol = colname,
                      nTree = NTREE_MF, nForest = NFOREST, nP = NOP, thd = THD)
soil_mfOp$varlist <- varlist
acc_all[[length(acc_all) + 1L]] <- rf_acc(soil_mfOp, soil_x, soil_y, i, "optimized_mf_rf", ACC_DIR)
save(soil_mfOp, file = file.path(MOD_DIR, paste0("soil_mfOp_zone", i, ".Rdata")))

rm(soil_ps, xy_y, xy_n, soil_y, soil_x, soil_plain, soil_mf,
   soil_zOp, soil_mfOp, soilList, varlist); gc()

# 6. RF training: zones 2:55, skipping zone 8 =================================

for (i in setdiff(2:55, 8)) {
  colname <- paste0("zone", i)
  
  if (!(colname %in% names(soil_hot))) {
    cat("[SKIP] no col:", colname, "\n"); next
  }
  
  n1 <- sum(soil_hot[[colname]] == 1, na.rm = TRUE)
  if (is.na(n1) || n1 == 0) {
    cat("[SKIP] no presence:", colname, "\n"); next
  }
  
  opfile <- file.path(OUT_DIR, paste0("soilopList_zone", i, ".csv"))
  if (!file.exists(opfile)) {
    cat("[SKIP] no opList:", colname, "\n"); next
  }
  
  tryCatch({
    n_rounds <- if (n1 >= 6000) 3L else if (n1 >= 3000) 2L else 1L
    cat("[RF] zone", i, "| n_presence =", n1, "| sampling rounds =", n_rounds, "\n")
    
    set.seed(BASE_SEED + i)
    soil_ps <- soil_hot
    for (rr in seq_len(n_rounds)) {
      soil_ps <- proSysSmpl(soil_ps,
                            byCol = which(colnames(soil_ps) == colname),
                            minSz = 2000)
    }
    soilList <- read.csv(opfile)
    if (nrow(soilList) < VAR_ROW + 1 || soilList[VAR_ROW + 1, 2] == "0") {
      cat("[SKIP] no valid varset:", colname, "\n"); next
    }
    varlist <- trimws(unlist(strsplit(soilList[VAR_ROW + 1, 2], ",")))
    varlist <- intersect(varlist, names(soil_ps))
    if (length(varlist) < 2) {
      cat("[SKIP] too few vars:", colname, "\n"); next
    }
    soil_ps <- soil_ps[complete.cases(soil_ps[, c(colname, varlist)]), ]
    
    n_pos <- sum(soil_ps[[colname]] == 1)
    n_neg <- sum(soil_ps[[colname]] == 0)
    if (n_neg > n_pos * 1.5) {
      idx_pos <- which(soil_ps[[colname]] == 1)
      idx_neg <- which(soil_ps[[colname]] == 0)
      idx_neg_sub <- sample(idx_neg, min(length(idx_neg), round(n_pos * 1.3)))
      soil_ps <- soil_ps[c(idx_pos, idx_neg_sub), ]
    }
    print(table(soil_ps[[colname]]))
    
    if (nrow(soil_ps) < 100 || length(unique(soil_ps[[colname]])) < 2) {
      cat("[SKIP] too few obs:", colname, "\n"); next
    }
    
    soil_ps <- as.data.frame(soil_ps)
    xy_y <- soil_ps[soil_ps[[colname]] == 1, ]
    xy_n <- soil_ps[soil_ps[[colname]] == 0, ]
    if (nrow(xy_y) < 2 || nrow(xy_n) < 2) {
      cat("[SKIP] too few pres/abs:", colname, "\n"); next
    }
    
    soil_y <- factor(soil_ps[[colname]], levels = c(0, 1))
    soil_x <- soil_ps[, varlist, drop = FALSE]
    
    soil_plain <- randomForest(soil_x, soil_y, ntree = NTREE_PLAIN, importance = TRUE)
    soil_plain$varlist <- varlist
    acc_all[[length(acc_all) + 1L]] <- rf_acc(soil_plain, soil_x, soil_y, i, "plain_rf", ACC_DIR)
    save(soil_plain, file = file.path(MOD_DIR, paste0("soil_plain_zone", i, ".Rdata")))
    
    soil_mf <- mcmfRF2(xy_y, xy_n, nr = MF_NR, varList = varlist, yCol = colname,
                       reg = FALSE, nTree = NTREE_MF, nForest = NFOREST)
    soil_mf$varlist <- varlist
    acc_all[[length(acc_all) + 1L]] <- rf_acc(soil_mf, soil_x, soil_y, i, "plain_mf_rf", ACC_DIR)
    save(soil_mf, file = file.path(MOD_DIR, paste0("soil_mf_zone", i, ".Rdata")))
    
    soil_zOp <- classOP(soil_x, soil_y, nTree1 = NTREE1, nTree2 = NTREE2, nOP = NOP, thd = THD)
    soil_zOp$varlist <- varlist
    acc_all[[length(acc_all) + 1L]] <- rf_acc(soil_zOp, soil_x, soil_y, i, "optimized_rf", ACC_DIR)
    save(soil_zOp, file = file.path(MOD_DIR, paste0("soil_zOp_zone", i, ".Rdata")))
    
    soil_mfOp <- mcmfRFop(xy_y, xy_n, nr = MF_NR, varList = varlist, yCol = colname,
                          nTree = NTREE_MF, nForest = NFOREST, nP = NOP, thd = THD)
    soil_mfOp$varlist <- varlist
    acc_all[[length(acc_all) + 1L]] <- rf_acc(soil_mfOp, soil_x, soil_y, i, "optimized_mf_rf", ACC_DIR)
    save(soil_mfOp, file = file.path(MOD_DIR, paste0("soil_mfOp_zone", i, ".Rdata")))
    
    cat("[DONE] zone", i, "\n")
    rm(soil_ps, xy_y, xy_n, soil_y, soil_x, soil_plain, soil_mf,
       soil_zOp, soil_mfOp, soilList, varlist); gc()
    
  }, error = function(e) {
    cat("[ERROR] RF zone", i, ":", conditionMessage(e), "\n"); gc()
  })
}

if (length(acc_all) > 0) {
  acc_all <- rbindlist(acc_all, fill = TRUE)
  fwrite(acc_all, file.path(ACC_DIR, "soil_rf_accuracy_summary.csv"))
  print(acc_all)
}
