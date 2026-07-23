# 2.16 Climate plain RFs using mcRFop variables
# ============================================================
# Purpose:
#   Train only the two corrected plain climate models:
#     1. rf_var: ordinary RF using the mcRFop variable set
#     2. mf_var: plain multi-Forest RF using the mcRFop variable set
#
# This script reuses outputs from script 2.15:
#   results/train_data.csv
#   results/train_combined_data_onehot.csv
#   results/clmhot_opList_zone*.csv
#
# Existing climate models are not overwritten.
# New model files:
#   rf/clm_rfVar_zone*.Rdata
#   rf/clm_mfVar_zone*.Rdata
# ============================================================

library(CEMT)
library(data.table)
library(foreach)
library(doSNOW)
library(randomForest)

rm(list = ls())
gc()

base_dir <- "H:/Jing/ecoChina2"
setwd(base_dir)


# 0. Plain multi-Forest function ================================================

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
    train_abs <- xy_n[
      sample(seq_len(nrow(xy_n)), n_abs, replace = FALSE),
      ,
      drop = FALSE
    ]
    
    train <- rbind(xy_y, train_abs)
    x2 <- train[, varList, drop = FALSE]
    
    if (!reg) y2 <- factor(train[[yCol]], levels = c(0, 1))
    if (reg)  y2 <- train[[yCol]]
    
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
    
    if (f == 1L) {
      rfC <- rf2
    } else {
      rfC <- combine(rfC, rf2)
    }
  }
  
  rfC
}


# 1. Accuracy helper =============================================================

rf_acc <- function(m, x, y, zone, model_name, out_dir) {
  y <- factor(y, levels = c(0, 1))
  
  oob_acc <- NA_real_
  
  if (!is.null(m$confusion)) {
    cm_oob <- as.matrix(
      m$confusion[, c("0", "1"), drop = FALSE]
    )
    
    oob_acc <- sum(diag(cm_oob)) / sum(cm_oob)
    
    fwrite(
      as.data.table(cm_oob, keep.rownames = "observed"),
      file.path(
        out_dir,
        paste0(
          model_name,
          "_zone",
          zone,
          "_confusion_oob.csv"
        )
      )
    )
  }
  
  pred <- predict(m, x, type = "response")
  pred <- factor(pred, levels = c(0, 1))
  
  cm_train <- table(
    observed = y,
    predicted = pred
  )
  
  train_acc <- sum(diag(cm_train)) / sum(cm_train)
  
  fwrite(
    as.data.table(cm_train),
    file.path(
      out_dir,
      paste0(
        model_name,
        "_zone",
        zone,
        "_confusion_train.csv"
      )
    )
  )
  
  data.table(
    zone = zone,
    model = model_name,
    n = length(y),
    oob_accuracy = oob_acc,
    train_accuracy = train_acc
  )
}


# 2. Parameters ==================================================================

SMPL_POS_RF <- 8000L
SMPL_PA      <- 1 / 1.3
SMPL_MAXN    <- 20000L
BASE_SEED    <- 49L
VAR_ROW      <- 20L
NTREE_PLAIN  <- 500L
NTREE_MF     <- 100L
NFOREST      <- 10L
MF_NR        <- 1.3

OUT_DIR <- "results"
MOD_DIR <- "rf"
ACC_DIR <- "accuracy_climate_var"

dir.create(MOD_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(ACC_DIR, showWarnings = FALSE, recursive = TRUE)


# 3. Reuse script 2.15 inputs ====================================================

train_file <- file.path(
  OUT_DIR,
  "train_data.csv"
)

onehot_file <- file.path(
  OUT_DIR,
  "train_combined_data_onehot.csv"
)

if (!file.exists(train_file)) {
  stop(
    "Missing: ",
    train_file,
    "\nRun script 2.15 first."
  )
}

if (!file.exists(onehot_file)) {
  stop(
    "Missing: ",
    onehot_file,
    "\nRun script 2.15 first."
  )
}

train_data <- fRead(train_file)
combined_data <- fRead(onehot_file)

xlist <- names(train_data)[6:ncol(train_data)]
zones <- 1:55

acc_all <- list()


# 4. Train only the two new models ==============================================

for (i in zones) {
  colname <- paste0("zone", i)
  
  if (!(colname %in% names(combined_data))) {
    cat("[SKIP] no one-hot column:", colname, "\n")
    next
  }
  
  n1_all <- sum(
    combined_data[[colname]] == 1,
    na.rm = TRUE
  )
  
  if (is.na(n1_all) || n1_all == 0) {
    cat("[SKIP] no presence:", colname, "\n")
    next
  }
  
  opfile <- file.path(
    OUT_DIR,
    paste0(
      "clmhot_opList_zone",
      i,
      ".csv"
    )
  )
  
  if (!file.exists(opfile)) {
    cat("[SKIP] no opList:", colname, "\n")
    next
  }
  
  tryCatch({
    pos_i <- min(
      SMPL_POS_RF,
      as.integer(n1_all)
    )
    
    noise_i <- max(
      1L,
      min(
        as.integer(round(0.10 * pos_i)),
        pos_i - 1L
      )
    )
    
    cols <- c(colname, xlist)
    
    cat("\n[RF VAR]", colname, "\n")
    
    dt_s <- smpl_pa(
      combined_data,
      colname,
      cols = cols,
      pos = pos_i,
      noise = noise_i,
      pa = SMPL_PA,
      max_n = SMPL_MAXN,
      seed = BASE_SEED + i
    )
    
    print(
      dt_s[, .N, by = get(colname)]
    )
    
    dt_s <- as.data.frame(dt_s)
    
    xy_y <- dt_s[
      dt_s[[colname]] == 1,
      ,
      drop = FALSE
    ]
    
    xy_n <- dt_s[
      dt_s[[colname]] == 0,
      ,
      drop = FALSE
    ]
    
    if (nrow(xy_y) < 2 || nrow(xy_n) < 2) {
      cat("[SKIP] too few pres/abs:", colname, "\n")
      next
    }
    
    clm_y <- factor(
      dt_s[[colname]],
      levels = c(0, 1)
    )
    
    clmList <- read.csv(opfile)
    
    if (
      nrow(clmList) < VAR_ROW + 1 ||
      clmList[VAR_ROW + 1, 2] == "0"
    ) {
      cat("[SKIP] no valid varset:", colname, "\n")
      next
    }
    
    varlist <- trimws(
      unlist(
        strsplit(
          clmList[VAR_ROW + 1, 2],
          ","
        )
      )
    )
    
    varlist <- intersect(
      varlist,
      names(dt_s)
    )
    
    if (length(varlist) < 2) {
      cat("[SKIP] too few vars:", colname, "\n")
      next
    }
    
    clm_x <- dt_s[
      ,
      varlist,
      drop = FALSE
    ]
    
    
    # 4.1 Ordinary RF using mcRFop variables ------------------------------------
    
    clm_rfVar <- randomForest(
      clm_x,
      clm_y,
      ntree = NTREE_PLAIN,
      importance = TRUE
    )
    
    clm_rfVar$varlist <- varlist
    
    acc_all[[length(acc_all) + 1L]] <- rf_acc(
      clm_rfVar,
      clm_x,
      clm_y,
      i,
      "rf_var",
      ACC_DIR
    )
    
    save(
      clm_rfVar,
      file = file.path(
        MOD_DIR,
        paste0(
          "clm_rfVar_zone",
          i,
          ".Rdata"
        )
      )
    )
    
    
    # 4.2 Plain multi-Forest using mcRFop variables -----------------------------
    
    clm_mfVar <- mcmfRF2(
      xy_y,
      xy_n,
      nr = MF_NR,
      varList = varlist,
      yCol = colname,
      reg = FALSE,
      nTree = NTREE_MF,
      nForest = NFOREST
    )
    
    clm_mfVar$varlist <- varlist
    
    acc_all[[length(acc_all) + 1L]] <- rf_acc(
      clm_mfVar,
      clm_x,
      clm_y,
      i,
      "mf_var",
      ACC_DIR
    )
    
    save(
      clm_mfVar,
      file = file.path(
        MOD_DIR,
        paste0(
          "clm_mfVar_zone",
          i,
          ".Rdata"
        )
      )
    )
    
    cat("[DONE]", colname, "\n")
    
    rm(
      dt_s,
      xy_y,
      xy_n,
      clm_y,
      clmList,
      varlist,
      clm_x,
      clm_rfVar,
      clm_mfVar
    )
    
    gc()
    
  }, error = function(e) {
    cat(
      "[ERROR]",
      colname,
      ":",
      conditionMessage(e),
      "\n"
    )
    
    gc()
  })
}


# 5. Save summary ================================================================

if (length(acc_all) > 0) {
  acc_all <- rbindlist(
    acc_all,
    fill = TRUE
  )
  
  fwrite(
    acc_all,
    file.path(
      ACC_DIR,
      "climate_rf_var_accuracy_summary.csv"
    )
  )
  
  print(acc_all)
}
