# Climate One-Hot RF — no spatial CV + multi-Forest
# Zone 1 first, then loop over zones 2:55
# ============================================================
library(CEMT)
library(ClimateNAr)
library(randomForest)
library(caret)
library(data.table)
library(foreach)
library(doSNOW)

rm(list = ls()); gc()

base_dir <- "H:/Jing/ecoChina2"
setwd(base_dir)

source(file.path(base_dir, "functions", "mcRFop_cls3.R"))

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

rf_acc <- function(m, x, y, zone, model_name, out_dir) {
  
  y <- factor(y, levels = c(0, 1))
  
  # OOB accuracy from randomForest object, if available
  oob_acc <- NA_real_
  if (!is.null(m$confusion)) {
    cm_oob <- as.matrix(m$confusion[, c("0", "1"), drop = FALSE])
    oob_acc <- sum(diag(cm_oob)) / sum(cm_oob)
    
    fwrite(
      as.data.table(cm_oob, keep.rownames = "observed"),
      file.path(out_dir, paste0(model_name, "_zone", zone, "_confusion_oob.csv"))
    )
  }
  
  # Training-set confusion matrix
  pred <- predict(m, x, type = "response")
  pred <- factor(pred, levels = c(0, 1))
  
  cm_train <- table(
    observed = y,
    predicted = pred
  )
  
  train_acc <- sum(diag(cm_train)) / sum(cm_train)
  
  fwrite(
    as.data.table(cm_train),
    file.path(out_dir, paste0(model_name, "_zone", zone, "_confusion_train.csv"))
  )
  
  data.table(
    zone = zone,
    model = model_name,
    n = length(y),
    oob_accuracy = oob_acc,
    train_accuracy = train_acc
  )
}


# ── params ─────────────────────────────────────
SMPL_POS_OP   <- 5000L
SMPL_POS_RF   <- 8000L
SMPL_PA       <- 1 / 1.3
SMPL_MAXN     <- 20000L
BASE_SEED     <- 49L
VAR_ROW       <- 20L
NTREE_OPLIST  <- 100L
NTREE_PLAIN   <- 500L
NTREE_MF      <- 100L
NFOREST       <- 10L
MF_NR         <- 1.3
NTREE1        <- 100L
NTREE2        <- 500L
NOP           <- 3L
THD           <- 0.75
OUT_DIR       <- "results"
MOD_DIR       <- "rf"
ACC_DIR       <- "accuracy_climate"
# ───────────────────────────────────────────────

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MOD_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(ACC_DIR, showWarnings = FALSE, recursive = TRUE)

# 1. Import climate dataframe ============================================================

dat <- fRead("data raw/1. zoneID_Clm_800m_Normal_1961_1990SY.csv"); hd(dat); names(dat)
dat[dat == -9999] <- NA
print(sort(unique(dat$zoneID)))

dat <- dat[complete.cases(dat$zoneID), ]
dat <- dat[complete.cases(dat[, 6:ncol(dat)]), ]
zone_counts <- table(dat$zoneID); print(zone_counts)

xlist <- colnames(dat)[6:ncol(dat)]
x <- dat[, xlist]; hd(x)

summary(x$AHM)
summary(x$Tave_sm)

y <- as.factor(dat$zoneID); levels(y)

agg_data <- aggregate(dat[, 6:ncol(dat)], by = list(zoneID = dat$zoneID), FUN = mean)
write.csv(agg_data, file.path(OUT_DIR, "summarize_climstat_by_zoneID.csv"), row.names = FALSE)

# 2. One-hot encoding ====================================================================

dat$zoneID <- as.factor(dat$zoneID)
one_hot <- model.matrix(~zoneID - 1, data = dat); hd(one_hot)
colnames(one_hot) <- gsub("zoneID", "zone", colnames(one_hot))
one_hot_df <- as.data.frame(one_hot)
combined_data <- cbind(dat, one_hot_df)
head(combined_data)

# 3. Variable selection: zone 1 first ====================================================

i <- 1
colname <- paste0("zone", i)
opfile <- file.path(OUT_DIR, paste0("clmhot_opList_zone", i, ".csv"))

n1_all <- sum(combined_data[[colname]] == 1, na.rm = TRUE)
pos_i <- min(SMPL_POS_OP, as.integer(n1_all))
noise_i <- max(1L, min(as.integer(round(0.10 * pos_i)), pos_i - 1L))
cols <- c(colname, xlist)

cat("[OPLIST] zone", i, "n1 =", n1_all, "\n")

if (!file.exists(opfile)) {
  dt_s <- smpl_pa(combined_data, colname,
                  cols = cols, pos = pos_i, noise = noise_i,
                  pa = SMPL_PA, max_n = SMPL_MAXN, seed = BASE_SEED + i)
  print(dt_s[, .N, by = get(colname)])
  
  x_df <- na.omit(as.data.frame(dt_s[, ..xlist]))
  y_fac <- factor(dt_s[[colname]][as.numeric(rownames(x_df))], levels = c(0, 1))
  
  clmopList <- mcRFop_cls(x_df, y_fac, nTree = NTREE_OPLIST)
  write.csv(as.data.frame(clmopList), opfile, row.names = FALSE)
  rm(dt_s, x_df, y_fac, clmopList); gc()
} else {
  cat("[SKIP] opList already exists: zone", i, "\n")
}

# 4. Variable selection: zones 2:55 ======================================================

for (i in 2:55) {
  colname <- paste0("zone", i)
  
  if (!(colname %in% names(combined_data))) {
    cat("[SKIP] no col:", colname, "\n"); next
  }
  
  n1_all <- sum(combined_data[[colname]] == 1, na.rm = TRUE)
  if (is.na(n1_all) || n1_all == 0) {
    cat("[SKIP] no presence:", colname, "\n"); next
  }
  
  opfile <- file.path(OUT_DIR, paste0("clmhot_opList_zone", i, ".csv"))
  if (file.exists(opfile)) {
    cat("[SKIP] opList already exists:", colname, "\n"); next
  }
  
  tryCatch({
    pos_i <- min(SMPL_POS_OP, as.integer(n1_all))
    noise_i <- max(1L, min(as.integer(round(0.10 * pos_i)), pos_i - 1L))
    cols <- c(colname, xlist)
    
    cat("[OPLIST] zone", i, "n1 =", n1_all, "\n")
    
    dt_s <- smpl_pa(combined_data, colname,
                    cols = cols, pos = pos_i, noise = noise_i,
                    pa = SMPL_PA, max_n = SMPL_MAXN, seed = BASE_SEED + i)
    print(dt_s[, .N, by = get(colname)])
    
    x_df <- na.omit(as.data.frame(dt_s[, ..xlist]))
    y_fac <- factor(dt_s[[colname]][as.numeric(rownames(x_df))], levels = c(0, 1))
    
    if (anyNA(y_fac) || length(unique(y_fac)) < 2) {
      cat("[SKIP] invalid y:", colname, "\n"); next
    }
    
    clmopList <- mcRFop_cls(x_df, y_fac, nTree = NTREE_OPLIST)
    write.csv(as.data.frame(clmopList), opfile, row.names = FALSE)
    rm(dt_s, x_df, y_fac, clmopList); gc()
    
  }, error = function(e) {
    cat("[ERROR] opList", colname, ":", conditionMessage(e), "\n"); gc()
  })
}

# 5. RF training: zone 1 first ===========================================================

acc_all <- list()

i <- 1
colname <- paste0("zone", i)
opfile <- file.path(OUT_DIR, paste0("clmhot_opList_zone", i, ".csv"))

n1_all <- sum(combined_data[[colname]] == 1, na.rm = TRUE)
pos_i <- min(SMPL_POS_RF, as.integer(n1_all))
noise_i <- max(1L, min(as.integer(round(0.10 * pos_i)), pos_i - 1L))
cols <- c(colname, xlist)

cat("[RF] zone", i, "\n")

dt_s <- smpl_pa(combined_data, colname,
                cols = cols, pos = pos_i, noise = noise_i,
                pa = SMPL_PA, max_n = SMPL_MAXN, seed = BASE_SEED + i)
print(dt_s[, .N, by = get(colname)])

dt_s <- as.data.frame(dt_s)
xy_y <- dt_s[dt_s[[colname]] == 1, ]
xy_n <- dt_s[dt_s[[colname]] == 0, ]
clm_y <- factor(dt_s[[colname]], levels = c(0, 1))
clm_x0 <- dt_s[, xlist, drop = FALSE]

# 5.1 Plain single RF
clm_plain <- randomForest(clm_x0, clm_y, ntree = NTREE_PLAIN, importance = TRUE)
clm_plain$varlist <- xlist
acc_all[[length(acc_all) + 1L]] <- rf_acc(clm_plain, clm_x0, clm_y, i, "plain_rf", ACC_DIR)
save(clm_plain, file = file.path(MOD_DIR, paste0("clm_plain_zone", i, ".Rdata")))

# 5.2 Plain multi-Forest RF
clm_mf <- mcmfRF2(xy_y, xy_n, nr = MF_NR, varList = xlist, yCol = colname,
                  reg = FALSE, nTree = NTREE_MF, nForest = NFOREST)
clm_mf$varlist <- xlist
acc_all[[length(acc_all) + 1L]] <- rf_acc(clm_mf, clm_x0, clm_y, i, "plain_mf_rf", ACC_DIR)
save(clm_mf, file = file.path(MOD_DIR, paste0("clm_mf_zone", i, ".Rdata")))

# 5.3 Optimized single RF
clmList <- read.csv(opfile)

varlist <- trimws(unlist(strsplit(clmList[VAR_ROW + 1, 2], ",")))
varlist <- intersect(varlist, names(dt_s))

clm_x <- dt_s[, varlist, drop = FALSE]
clim_zOp <- classOP(clm_x, clm_y, nTree1 = NTREE1, nTree2 = NTREE2, nOP = NOP, thd = THD)
clim_zOp$varlist <- varlist
acc_all[[length(acc_all) + 1L]] <- rf_acc(clim_zOp, clm_x, clm_y, i, "optimized_rf", ACC_DIR)
save(clim_zOp, file = file.path(MOD_DIR, paste0("clm_zOp_zone", i, ".Rdata")))

# 5.4 Optimized multi-Forest RF
clim_mfOp <- mcmfRFop(xy_y, xy_n, nr = MF_NR, varList = varlist, yCol = colname,
                      nTree = NTREE_MF, nForest = NFOREST, nP = NOP, thd = THD)
clim_mfOp$varlist <- varlist
acc_all[[length(acc_all) + 1L]] <- rf_acc(clim_mfOp, clm_x, clm_y, i, "optimized_mf_rf", ACC_DIR)
save(clim_mfOp, file = file.path(MOD_DIR, paste0("clm_mfOp_zone", i, ".Rdata")))

rm(dt_s, xy_y, xy_n, clm_y, clm_x0, clm_plain, clm_mf,
   clm_x, clim_zOp, clim_mfOp, clmList, varlist); gc()

# 6. RF training: zones 2:55 =============================================================

for (i in 2:55) {
  colname <- paste0("zone", i)
  
  if (!(colname %in% names(combined_data))) {
    cat("[SKIP] no col:", colname, "\n"); next
  }
  
  n1_all <- sum(combined_data[[colname]] == 1, na.rm = TRUE)
  if (is.na(n1_all) || n1_all == 0) {
    cat("[SKIP] no presence:", colname, "\n"); next
  }
  
  opfile <- file.path(OUT_DIR, paste0("clmhot_opList_zone", i, ".csv"))
  if (!file.exists(opfile)) {
    cat("[SKIP] no opList:", colname, "\n"); next
  }
  
  tryCatch({
    pos_i <- min(SMPL_POS_RF, as.integer(n1_all))
    noise_i <- max(1L, min(as.integer(round(0.10 * pos_i)), pos_i - 1L))
    cols <- c(colname, xlist)
    
    cat("[RF] zone", i, "\n")
    
    dt_s <- smpl_pa(combined_data, colname,
                    cols = cols, pos = pos_i, noise = noise_i,
                    pa = SMPL_PA, max_n = SMPL_MAXN, seed = BASE_SEED + i)
    print(dt_s[, .N, by = get(colname)])
    
    dt_s <- as.data.frame(dt_s)
    xy_y <- dt_s[dt_s[[colname]] == 1, ]
    xy_n <- dt_s[dt_s[[colname]] == 0, ]
    if (nrow(xy_y) < 2 || nrow(xy_n) < 2) {
      cat("[SKIP] too few pres/abs:", colname, "\n"); next
    }
    
    clm_y <- factor(dt_s[[colname]], levels = c(0, 1))
    clm_x0 <- dt_s[, xlist, drop = FALSE]
    
    clm_plain <- randomForest(clm_x0, clm_y, ntree = NTREE_PLAIN, importance = TRUE)
    clm_plain$varlist <- xlist
    acc_all[[length(acc_all) + 1L]] <- rf_acc(clm_plain, clm_x0, clm_y, i, "plain_rf", ACC_DIR)
    save(clm_plain, file = file.path(MOD_DIR, paste0("clm_plain_zone", i, ".Rdata")))
    
    clm_mf <- mcmfRF2(xy_y, xy_n, nr = MF_NR, varList = xlist, yCol = colname,
                      reg = FALSE, nTree = NTREE_MF, nForest = NFOREST)
    clm_mf$varlist <- xlist
    acc_all[[length(acc_all) + 1L]] <- rf_acc(clm_mf, clm_x0, clm_y, i, "plain_mf_rf", ACC_DIR)
    save(clm_mf, file = file.path(MOD_DIR, paste0("clm_mf_zone", i, ".Rdata")))
    
    clmList <- read.csv(opfile)
    if (nrow(clmList) < VAR_ROW + 1 || clmList[VAR_ROW + 1, 2] == "0") {
      cat("[SKIP] no valid varset:", colname, "\n"); next
    }
    varlist <- trimws(unlist(strsplit(clmList[VAR_ROW + 1, 2], ",")))
    varlist <- intersect(varlist, names(dt_s))
    if (length(varlist) < 2) {
      cat("[SKIP] too few vars:", colname, "\n"); next
    }
    
    clm_x <- dt_s[, varlist, drop = FALSE]
    clim_zOp <- classOP(clm_x, clm_y, nTree1 = NTREE1, nTree2 = NTREE2, nOP = NOP, thd = THD)
    clim_zOp$varlist <- varlist
    acc_all[[length(acc_all) + 1L]] <- rf_acc(clim_zOp, clm_x, clm_y, i, "optimized_rf", ACC_DIR)
    save(clim_zOp, file = file.path(MOD_DIR, paste0("clm_zOp_zone", i, ".Rdata")))
    
    clim_mfOp <- mcmfRFop(xy_y, xy_n, nr = MF_NR, varList = varlist, yCol = colname,
                          nTree = NTREE_MF, nForest = NFOREST, nP = NOP, thd = THD)
    clim_mfOp$varlist <- varlist
    acc_all[[length(acc_all) + 1L]] <- rf_acc(clim_mfOp, clm_x, clm_y, i, "optimized_mf_rf", ACC_DIR)
    save(clim_mfOp, file = file.path(MOD_DIR, paste0("clm_mfOp_zone", i, ".Rdata")))
    
    cat("[DONE]", colname, "\n")
    rm(dt_s, xy_y, xy_n, clm_y, clm_x0, clm_plain, clm_mf,
       clm_x, clim_zOp, clim_mfOp, clmList, varlist); gc()
    
  }, error = function(e) {
    cat("[ERROR] RF", colname, ":", conditionMessage(e), "\n"); gc()
  })
}

if (length(acc_all) > 0) {
  acc_all <- rbindlist(acc_all, fill = TRUE)
  fwrite(acc_all, file.path(ACC_DIR, "climate_rf_accuracy_summary.csv"))
  print(acc_all)
}
