# Soil One-Hot RF — Part 1: Data prep + variable selection
# ============================================================
library(CEMT)
library(raster)
library(randomForest)
library(caret)
library(data.table)
library(jsonlite)
library(sf)


rm(list = ls()); gc()

# ── params (adjust here) ─────────────────────────────────────
VAR_ROW   <- 8L      # row in opList for variable set (→ row 9, 0-indexed)
NTREE_OP  <- 100L    # trees for mcRFop_cls
BASE_SEED <- 49L
OUT_DIR   <- "H:/Jing/ecoChina2/results"
# ─────────────────────────────────────────────────────────────

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# 1. Import & prep ============================================================

setwd("H:/Jing/soil rasters")
soil_dat2 <- fRead('new_soil_raster.csv'); hd(soil_dat2)
soil_dat2 <- soil_dat2[, c(1:6, 8:22)]  # keep topsoil cols only
names(soil_dat2)

soil_train <- soil_dat2[complete.cases(soil_dat2[, c(2, 7:ncol(soil_dat2))]), ]
cat("Rows after NA removal:", nrow(soil_train), "\n")
cat("Zones present:", sort(unique(soil_train$zoneID)), "\n")  # zone 8 absent

# soil vars only (strip ID + spatial cols); keep x,y separately for spatial CV
soil_train2 <- soil_train[, -c(1, 3:6)]   # col1 = zoneID, cols 2-16 = 15 soil vars
soil_coords  <- soil_train[, c("x", "y")] # saved for Part 2 spatial CV

xlist <- names(soil_train2)[2:16]          # 15 soil variable names
cat("Soil predictors:", xlist, "\n")


# 2. One-hot encoding =========================================================

soil_train2$zone <- as.character(soil_train2$zoneID)
soil_hot_mat <- predict(dummyVars(~zone, data = soil_train2), newdata = soil_train2)
soil_hot <- cbind(soil_train2, soil_hot_mat, soil_coords)  # x, y appended for later

zones_all <- sort(setdiff(unique(soil_train$zoneID), 8))   # 1:55 minus zone 8
cat("Total zones to process:", length(zones_all), "\n")


# 3. Variable selection — zone 1 test =========================================

i       <- 1
colname <- paste0("zone", i)

n1 <- sum(soil_hot[[colname]] == 1, na.rm = TRUE)
cat("[TEST] zone", i, "| n_presence =", n1, "\n")

set.seed(BASE_SEED)
n_rounds <- if (n1 >= 6000) 3L else if (n1 >= 3000) 2L else 1L
cat(" sampling rounds:", n_rounds, "\n")

set.seed(BASE_SEED + i)
soil_ps <- soil_hot
for (r in seq_len(n_rounds)) {
  soil_ps <- proSysSmpl(soil_ps,
                        byCol = which(colnames(soil_ps) == colname),
                        minSz = 2000)
}
print(table(soil_ps[[colname]]))

soil_y <- factor(soil_ps[[colname]], levels = c(0, 1))
soil_x <- soil_ps[, xlist]

soilopList <- mcRFop_cls(soil_x, soil_y, nTree = NTREE_OP)
print(soilopList)

write.csv(as.data.frame(soilopList),
          file.path(OUT_DIR, paste0("soilopList_zone", i, ".csv")),
          row.names = FALSE)

cat("[DONE] zone", i, "\n")
rm(soil_ps, soil_y, soil_x, soilopList); gc()


# 4. Variable selection — zones 2:55 (skip 8) =================================

for (i in setdiff(2:55, 8)) {
  
  colname <- paste0("zone", i)
  
  if (!(colname %in% names(soil_hot))) {
    cat("[SKIP] column not found:", colname, "\n"); next
  }
  
  n1 <- sum(soil_hot[[colname]] == 1, na.rm = TRUE)
  if (is.na(n1) || n1 == 0) {
    cat("[SKIP] no presence:", colname, "\n"); next
  }
  
  opfile <- file.path(OUT_DIR, paste0("soilopList_zone", i, ".csv"))
  if (file.exists(opfile)) {
    cat("[SKIP] opList already exists:", colname, "\n"); next
  }
  
  cat("[RUN] zone", i, "| n_presence =", n1, "\n")
  
  tryCatch({
    set.seed(BASE_SEED + i)
    n_rounds <- if (n1 >= 6000) 3L else if (n1 >= 3000) 2L else 1L
    cat(" sampling rounds:", n_rounds, "\n")
    
    soil_ps <- soil_hot
    for (r in seq_len(n_rounds)) {
      soil_ps <- proSysSmpl(soil_ps,
                            byCol = which(colnames(soil_ps) == colname),
                            minSz = 2000)
    }

    tab <- table(soil_ps[[colname]])
    if (length(tab) < 2) {
      cat("[SKIP] single class after sampling:", colname, "\n"); next
    }
    print(tab)
    
    soil_y <- factor(soil_ps[[colname]], levels = c(0, 1))
    soil_x <- soil_ps[, xlist]
    
    if (nrow(soil_x) < 10) {
      cat("[SKIP] too few rows:", colname, "\n"); next
    }
    
    soilopList <- mcRFop_cls(soil_x, soil_y, nTree = NTREE_OP)
    
    write.csv(as.data.frame(soilopList),
              file.path(OUT_DIR, paste0("soilopList_zone", i, ".csv")),
              row.names = FALSE)
    
    cat("[DONE] zone", i, "\n")
    rm(soil_ps, soil_y, soil_x, soilopList); gc()
    
  }, error = function(e) {
    cat("[ERROR] zone", i, ":", conditionMessage(e), "\n")
    gc()
  })
}


# Soil One-Hot RF — Part 2: Spatial block CV + weighted ensemble
# ============================================================
# soil_hot and soil_train2 are still in memory from Part 1.

# ── params ─────────────────────────────────────
VAR_ROW  <- 8L        # row in opList for variable set (row 9, 0-indexed)
NTREE1   <- 100L      # trees per classOP optimization round
NTREE2   <- 500L      # trees for classOP final run
NOP      <- 3L        # classOP optimization rounds
BASE_SEED <- 49L

K_FOLD   <- 5L
NBX      <- 8L        # longitude blocks
NBY      <- 6L        # latitude blocks
CV_SEED  <- 4901L
THD      <- 0.75      # classOP label-cleaning threshold
BACC_MIN <- 0.55
SENS_MIN <- 0.30
SPEC_MIN <- 0.30

RES_DIR  <- "H:/Jing/ecoChina2/results"           # soilopList files from Part 1
OUT_DIR  <- "H:/Jing/ecoChina2/results_cv_soil"   # CV stats output
MOD_DIR  <- "H:/Jing/ecoChina2/rf_final_soil"     # ensemble models output
# ─────────────────────────────────────────────────────────────

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MOD_DIR, showWarnings = FALSE, recursive = TRUE)

xlist   <- names(soil_train2)[2:16]    # 15 soil variable names
sum_csv <- file.path(OUT_DIR, "cv_summary_allzones.csv")
if (file.exists(sum_csv)) file.remove(sum_csv)


# ── helper functions ─────────────────────────────────────────

make_sp_folds <- function(df, xcol = "x", ycol = "y",
                          nbx = 8L, nby = 6L, K = 5L, seed = 49L) {
  bx  <- cut(df[[xcol]], nbx, labels = FALSE)
  by  <- cut(df[[ycol]], nby, labels = FALSE)
  blk <- (bx - 1L) * nby + by
  set.seed(seed)
  ublk <- sample(unique(blk))
  fmap <- setNames(rep_len(1:K, length(ublk)), ublk)
  as.integer(fmap[as.character(blk)])
}

fold_eval <- function(m, te, colname, varlist, thr = 0.5) {
  p    <- predict(m, te[, varlist, drop = FALSE], type = "prob")[, "1"]
  y01  <- as.integer(as.character(te[[colname]]) == "1")
  pred <- as.integer(p >= thr)
  tp   <- sum(pred == 1L & y01 == 1L)
  tn   <- sum(pred == 0L & y01 == 0L)
  fp   <- sum(pred == 1L & y01 == 0L)
  fn   <- sum(pred == 0L & y01 == 1L)
  sens <- tp / max(1L, tp + fn)
  spec <- tn / max(1L, tn + fp)
  list(n_pos = sum(y01),
       n_neg = length(y01) - sum(y01),
       acc   = (tp + tn) / max(1L, tp + tn + fp + fn),
       sens  = sens, spec = spec,
       bacc  = 0.5 * (sens + spec))
}

predict.rfEnsemble <- function(object, newdata, ...) {
  X     <- as.data.frame(newdata[, object$varlist, drop = FALSE])
  plist <- lapply(object$models,
                  function(m) predict(m, X, type = "prob")[, "1"])
  as.numeric(do.call(cbind, plist) %*% object$weights)
}


# 5. Spatial CV — zone 1 test =================================================

i       <- 1
colname <- paste0("zone", i)

# 5.1 Load varlist
opfile  <- file.path(RES_DIR, paste0("soilopList_zone", i, ".csv"))
soilList <- read.csv(opfile)
varlist  <- trimws(unlist(strsplit(soilList[VAR_ROW + 1, 2], ",")))
varlist  <- intersect(varlist, names(soil_hot))
cat("[TEST] zone", i, "| vars =", length(varlist), "\n")

# 5.2 Sample
set.seed(BASE_SEED)
n1 <- sum(soil_hot[[colname]] == 1, na.rm = TRUE)
n_rounds <- if (n1 >= 6000) 3L else if (n1 >= 3000) 2L else 1L
cat(" sampling rounds:", n_rounds, "\n")

soil_ps <- soil_hot
for (r in seq_len(n_rounds)) {
  soil_ps <- proSysSmpl(soil_ps,
                        byCol = which(colnames(soil_ps) == colname),
                        minSz = 2000)
}
soil_ps <- soil_ps[complete.cases(soil_ps[, c(colname, "x", "y", varlist)]), ]

n_pos <- sum(soil_ps[[colname]] == 1)
n_neg <- sum(soil_ps[[colname]] == 0)
if (n_neg > n_pos * 1.5) { # absence 超过 presence 的1.5倍才触发 
  idx_pos <- which(soil_ps[[colname]] == 1) 
  idx_neg <- which(soil_ps[[colname]] == 0) 
  idx_neg_sub <- sample(idx_neg, min(length(idx_neg), round(n_pos * 1.3))) 
  soil_ps <- soil_ps[c(idx_pos, idx_neg_sub), ]
  }
print(table(soil_ps[[colname]]))

# 5.3 Assign spatial folds
soil_ps$fold <- make_sp_folds(soil_ps, xcol = "x", ycol = "y",
                              nbx = NBX, nby = NBY, K = K_FOLD,
                              seed = CV_SEED + i)

# 5.4 Train + evaluate per fold
fold_models <- vector("list", K_FOLD)
fold_stats  <- vector("list", K_FOLD)

for (k in 1:K_FOLD) {
  tr <- soil_ps[soil_ps$fold != k, ]
  te <- soil_ps[soil_ps$fold == k, ]
  
  if (nrow(te) < 50 || length(unique(te[[colname]])) < 2) next
  
  m  <- classOP(tr[, varlist], factor(tr[[colname]], levels = c(0, 1)),
                nTree1 = NTREE1, nTree2 = NTREE2, nOP = NOP, thd = THD)
  ev <- fold_eval(m, te, colname, varlist)
  
  fold_models[[k]] <- m
  fold_stats[[k]]  <- data.frame(zone = i, fold = k,
                                 n_tr = nrow(tr), n_te = nrow(te),
                                 n_pos = ev$n_pos, n_neg = ev$n_neg,
                                 acc = ev$acc, sens = ev$sens,
                                 spec = ev$spec, bacc = ev$bacc)
  
  cat("  fold", k, "| te:", ev$n_pos, "+/", ev$n_neg,
      "- | sens =", round(ev$sens, 3),
      " spec =", round(ev$spec, 3),
      " bacc =", round(ev$bacc, 3), "\n")
}

# 5.5 Save CV stats
cv_res <- do.call(rbind, fold_stats)
write.csv(cv_res,
          file.path(OUT_DIR, paste0("cv_zone", i, ".csv")),
          row.names = FALSE)

cv_sum <- data.frame(zone = i,
                     n_te  = sum(cv_res$n_te),
                     acc   = mean(cv_res$acc),
                     sens  = mean(cv_res$sens),
                     spec  = mean(cv_res$spec),
                     bacc  = mean(cv_res$bacc))
cv_sum$flag <- ifelse(cv_sum$bacc < BACC_MIN, "WARN", "OK")
write.csv(cv_sum, sum_csv, row.names = FALSE)
print(cv_sum)

# 5.6 Build weighted ensemble
keep <- cv_res$fold[cv_res$sens >= SENS_MIN &
                      cv_res$spec >= SPEC_MIN &
                      cv_res$bacc >= BACC_MIN]
if (length(keep) == 0) {
  cat("[WARN] no fold passed thresholds; keeping best bacc fold\n")
  keep <- cv_res$fold[which.max(cv_res$bacc)]
}
keep <- sort(unique(keep))

mods <- fold_models[keep]
w    <- cv_res$bacc[match(keep, cv_res$fold)]
ok   <- !vapply(mods, is.null, logical(1))
mods <- mods[ok]; w <- w[ok]
w    <- w / sum(w)

soil_ens <- list(models = mods, weights = w, varlist = varlist,
                 cv = cv_res, keep = keep[ok])
class(soil_ens) <- "rfEnsemble"

save(soil_ens, file = file.path(MOD_DIR, paste0("soil_ens_zone", i, ".Rdata")))
cat("[DONE] zone", i, "| kept folds:", paste(keep[ok], collapse = ","),
    "| weights:", paste(round(w, 3), collapse = ","), "\n")

rm(soil_ps, cv_res, cv_sum, fold_models, fold_stats, mods, w, soil_ens, soilList, varlist)
gc()


# 6. Spatial CV — zones 2:55 (skip 8) =========================================

for (i in setdiff(2:55, 8)) {
  
  colname <- paste0("zone", i)
  
  if (!(colname %in% names(soil_hot)))                    { cat("[SKIP] no col:", colname, "\n"); next }
  
  n1 <- sum(soil_hot[[colname]] == 1, na.rm = TRUE)
  if (is.na(n1) || n1 == 0)                               { cat("[SKIP] no presence:", colname, "\n"); next }
  
  opfile <- file.path(RES_DIR, paste0("soilopList_zone", i, ".csv"))
  if (!file.exists(opfile))                               { cat("[SKIP] no opList:", colname, "\n"); next }
  
  soilList <- read.csv(opfile)
  if (nrow(soilList) < VAR_ROW + 1 || soilList[VAR_ROW + 1, 2] == "0") {
    cat("[SKIP] no valid varset:", colname, "\n"); next
  }
  
  varlist <- trimws(unlist(strsplit(soilList[VAR_ROW + 1, 2], ",")))
  varlist <- intersect(varlist, names(soil_hot))
  if (length(varlist) < 2)                                { cat("[SKIP] too few vars:", colname, "\n"); next }
  
  tryCatch({
    
    # 6.1 Sample
    set.seed(BASE_SEED)
    n_rounds <- if (n1 >= 6000) 3L else if (n1 >= 3000) 2L else 1L
    cat(" sampling rounds:", n_rounds, "\n")
    
    soil_ps <- soil_hot
    for (r in seq_len(n_rounds)) {
      soil_ps <- proSysSmpl(soil_ps,
                            byCol = which(colnames(soil_ps) == colname),
                            minSz = 2000)
    }
    soil_ps <- soil_ps[complete.cases(soil_ps[, c(colname, "x", "y", varlist)]), ]
    
    n_pos <- sum(soil_ps[[colname]] == 1)
    n_neg <- sum(soil_ps[[colname]] == 0)
    if (n_neg > n_pos * 1.5) { # absence 超过 presence 的1.5倍才触发 
      idx_pos <- which(soil_ps[[colname]] == 1) 
      idx_neg <- which(soil_ps[[colname]] == 0) 
      idx_neg_sub <- sample(idx_neg, min(length(idx_neg), round(n_pos * 1.3))) 
      soil_ps <- soil_ps[c(idx_pos, idx_neg_sub), ]
    }
    print(table(soil_ps[[colname]]))
    
    if (nrow(soil_ps) < 100 || length(unique(soil_ps[[colname]])) < 2) {
      cat("[SKIP] too few obs:", colname, "\n"); next
    }
    
    cat("[RUN] zone", i, "| n1 =", n1,
        "| sampled:", sum(soil_ps[[colname]] == 1), "+/",
        sum(soil_ps[[colname]] == 0), "-",
        "| vars =", length(varlist), "\n")
    
    # 6.2 Assign spatial folds
    soil_ps$fold <- make_sp_folds(soil_ps, xcol = "x", ycol = "y",
                                  nbx = NBX, nby = NBY, K = K_FOLD,
                                  seed = CV_SEED + i)
    
    # 6.3 Train + evaluate per fold
    fold_models <- vector("list", K_FOLD)
    fold_stats  <- vector("list", K_FOLD)
    
    for (k in 1:K_FOLD) {
      tr <- soil_ps[soil_ps$fold != k, ]
      te <- soil_ps[soil_ps$fold == k, ]
      
      if (nrow(te) < 50 || length(unique(te[[colname]])) < 2) next
      
      m  <- classOP(tr[, varlist], factor(tr[[colname]], levels = c(0, 1)),
                    nTree1 = NTREE1, nTree2 = NTREE2, nOP = NOP, thd = THD)
      ev <- fold_eval(m, te, colname, varlist)
      
      fold_models[[k]] <- m
      fold_stats[[k]]  <- data.frame(zone = i, fold = k,
                                     n_tr = nrow(tr), n_te = nrow(te),
                                     n_pos = ev$n_pos, n_neg = ev$n_neg,
                                     acc = ev$acc, sens = ev$sens,
                                     spec = ev$spec, bacc = ev$bacc)
      
      cat("  fold", k, "| te:", ev$n_pos, "+/", ev$n_neg,
          "- | sens =", round(ev$sens, 3),
          " spec =", round(ev$spec, 3),
          " bacc =", round(ev$bacc, 3), "\n")
    }
    
    # 6.4 Save CV stats
    cv_res <- do.call(rbind, fold_stats)
    if (is.null(cv_res) || nrow(cv_res) == 0) { cat("[SKIP] no valid folds:", colname, "\n"); next }
    
    write.csv(cv_res,
              file.path(OUT_DIR, paste0("cv_zone", i, ".csv")),
              row.names = FALSE)
    
    cv_sum <- data.frame(zone = i,
                         n_te  = sum(cv_res$n_te),
                         acc   = mean(cv_res$acc),
                         sens  = mean(cv_res$sens),
                         spec  = mean(cv_res$spec),
                         bacc  = mean(cv_res$bacc))
    cv_sum$flag <- ifelse(cv_sum$bacc < BACC_MIN, "WARN", "OK")
    write.csv(cv_sum, sum_csv, append = file.exists(sum_csv), row.names = FALSE)
    
    # 6.5 Build weighted ensemble
    keep <- cv_res$fold[cv_res$sens >= SENS_MIN &
                          cv_res$spec >= SPEC_MIN &
                          cv_res$bacc >= BACC_MIN]
    if (length(keep) == 0) {
      cat("[WARN] zone", i, "no fold passed; keeping best bacc fold\n")
      keep <- cv_res$fold[which.max(cv_res$bacc)]
    }
    keep <- sort(unique(keep))
    
    mods <- fold_models[keep]
    w    <- cv_res$bacc[match(keep, cv_res$fold)]
    ok   <- !vapply(mods, is.null, logical(1))
    mods <- mods[ok]; w <- w[ok]
    wsum <- sum(w, na.rm = TRUE)
    w    <- if (!is.finite(wsum) || wsum <= 0) rep(1 / length(w), length(w)) else w / wsum
    
    soil_ens <- list(models = mods, weights = w, varlist = varlist,
                     cv = cv_res, keep = keep[ok])
    class(soil_ens) <- "rfEnsemble"
    
    save(soil_ens,
         file = file.path(MOD_DIR, paste0("soil_ens_zone", i, ".Rdata")))
    
    cat("[DONE] zone", i,
        "| kept:", paste(keep[ok], collapse = ","),
        "| w:", paste(round(w, 3), collapse = ","), "\n")
    
    rm(soil_ps, cv_res, cv_sum, fold_models, fold_stats,
       mods, w, soil_ens, soilList, varlist)
    gc()
    
  }, error = function(e) {
    cat("[ERROR] zone", i, ":", conditionMessage(e), "\n")
    gc()
  })
}

