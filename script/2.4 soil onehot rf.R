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
OUT_DIR   <- "E:/Jing/ecoChina2/results"
# ─────────────────────────────────────────────────────────────

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# 1. Import & prep ============================================================

setwd("E:/Jing/soil rasters")
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
soil_ps <- proSysSmpl(soil_hot, byCol = which(colnames(soil_hot) == colname), minSz = 2000)
soil_ps <- proSysSmpl(soil_ps,  byCol = which(colnames(soil_ps)  == colname), minSz = 2000)
soil_ps <- proSysSmpl(soil_ps,  byCol = which(colnames(soil_ps)  == colname), minSz = 2000)
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
    soil_ps <- proSysSmpl(soil_hot, byCol = which(colnames(soil_hot) == colname), minSz = 2000)
    soil_ps <- proSysSmpl(soil_ps,  byCol = which(colnames(soil_ps)  == colname), minSz = 2000)
    soil_ps <- proSysSmpl(soil_ps,  byCol = which(colnames(soil_ps)  == colname), minSz = 2000)
    
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