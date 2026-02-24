# Climate One-Hot classification rf#
library(CEMT)
library(terra)
library(ClimateNAr)

library(randomForest)
library(caret)
library(jsonlite)
library(sf)
library(data.table)

rm();gc()
setwd("H:/Jing/ecoChina2")

# 1. Import climate dataframe ============================================================

#rf data
dat <- fRead("data raw/1. zoneID_Clm_800m_Normal_1961_1990SY.csv");hd(dat);names(dat)
dat[dat == -9999] <- NA
print(sort(unique(dat$zoneID))) #has zone 56

dat <- dat[complete.cases(dat$zoneID), ]
hd(dat) #13845898 obs

dat <- dat[complete.cases(dat[,6:ncol(dat)]), ];hd(dat) #13806550   obs
zone_counts <- table(dat$zoneID);print(zone_counts) #what's deleted does not matter

xlist <- colnames(dat)[6:ncol(dat)]
x <- dat[,xlist]; hd(x)

summary(x$AHM) #max 3561.30
summary(x$Tave_sm) #max 32.00

y <- as.factor(dat$zoneID);levels(y) # zone 56 disappeared due to missing climate data

# Aggregate data by 'zoneID'
agg_data <- aggregate(dat[, 6:ncol(dat)], by=list(zoneID=dat$zoneID), FUN=mean)

# Write the aggregated data to a CSV file (mean climate of each zone)
write.csv(agg_data, "results/summarize_climstat_by_zoneID.csv", row.names=FALSE) #summary all clim stats by zoneID



# 2. One-hot Encoding====================================================================
dat$zoneID <- as.factor(dat$zoneID)
one_hot <- model.matrix(~zoneID - 1, data = dat);hd(one_hot)
colnames(one_hot) <- gsub("zoneID", "zone", colnames(one_hot))
one_hot_df <- as.data.frame(one_hot);colnames(one_hot_df)

# Combine the original data with the one-hot encoded columns
combined_data <- cbind(dat, one_hot_df)
head(combined_data)



# 3. Compressing & sampling =============================================================
# smpl_pa params (adjust here for tuning) ---------
SMPL_POS   <- 5000L      # target presence (cap by available n1)
SMPL_PA    <- 1/1.3       # P/A ratio (presence/absence)
SMPL_MAXN  <- 20000L      # max total sample size per zone
BASE_SEED  <- 49L

# 3.1 zone1，presence ≈ 5000±10%，P:A=1:1.3 ---------
colname <- paste0("zone", 1)

# predictors
xlist <- names(combined_data)[6:74]
cols  <- c(colname, xlist)

# per-zone pos/noise (noise = 10% of pos)
n1_all  <- sum(combined_data[[colname]] == 1, na.rm = TRUE)
pos_i   <- min(SMPL_POS, as.integer(n1_all))
noise_i <- as.integer(round(0.10 * pos_i))
noise_i <- min(noise_i, pos_i - 1L)

dt_s <- smpl_pa(combined_data, colname,
                cols  = cols,
                pos   = pos_i,
                noise = noise_i,
                pa    = SMPL_PA,
                max_n = SMPL_MAXN,
                seed  = BASE_SEED)

dt_s[, .N, by = get(colname)]

clm_y <- factor(dt_s[[colname]], levels = c(0, 1))
print(table(clm_y))

# x for mcRFop_cls: need data.frame + no NA
x_df  <- na.omit(as.data.frame(dt_s[, ..xlist]))
y_fac <- factor(dt_s[[colname]][as.numeric(rownames(x_df))], levels = c(0, 1))

clmopList <- mcRFop_cls(x_df, y_fac, nTree=100)

# make sure output dir exists
dir.create("results", showWarnings = FALSE, recursive = TRUE)

# build a single filename string
outfile <- paste0("results/clmhot_opList_zone", 1, ".csv")

# safest: convert to data.frame and write.csv (always works)
write.csv(as.data.frame(clmopList), outfile, row.names = FALSE)

# 3.2 Climate var optimize, zone 2-55 ------

for (i in 2:55){
  
  colname <- paste0("zone", i)
  
  # skip if no one-hot column
  if(!(colname %in% names(combined_data))){
    print(paste("[SKIP] column not found:", colname))
    next
  }
  
  # skip if no presence
  n1_all <- sum(combined_data[[colname]] == 1, na.rm = TRUE)
  if(is.na(n1_all) || n1_all == 0){
    print(paste("[SKIP] no 1s in:", colname))
    next
  }
  
  # ---------------- sampling parameters ----------------
  pos_i   <- min(SMPL_POS, as.integer(n1_all))
  
  noise_i <- as.integer(round(0.10 * pos_i))
  noise_i <- max(1L, min(noise_i, pos_i - 1L))   
  
  cols <- c(colname, names(combined_data)[6:74])
  
  print(paste("[RUN]", colname,
              "n1=", n1_all,
              "pos=", pos_i,
              "noise~", noise_i,
              "pa=", SMPL_PA,
              "max_n=", SMPL_MAXN))
  
  tryCatch({
    
    # ---------------- smpl_pa ----------------
    dt_s <- smpl_pa(combined_data, colname,
                    cols  = cols,
                    pos   = pos_i,
                    noise = noise_i,
                    pa    = SMPL_PA,
                    max_n = SMPL_MAXN,
                    seed  = BASE_SEED + i)
    
    # check sampled class balance
    tab <- dt_s[, .N, by = get(colname)]
    if(nrow(tab) < 2){
      print(paste("[SKIP] sampled data has only one class:", colname))
      next
    }
    
    # ---------------- prepare y ----------------
    clm_y <- factor(dt_s[[colname]], levels = c(0, 1))
    print(table(clm_y))
    
    # ---------------- prepare x ----------------
    x_df  <- na.omit(as.data.frame(dt_s[, ..xlist]))
    if(nrow(x_df) < 10){
      print(paste("[SKIP] too few rows after na.omit:", colname))
      next
    }
    
    y_fac <- factor(dt_s[[colname]][as.numeric(rownames(x_df))],
                    levels = c(0,1))
    
    if(anyNA(y_fac) || length(unique(y_fac)) < 2){
      print(paste("[SKIP] y invalid / single class after NA removal:", colname))
      next
    }
    
    # ---------------- mcRF optimize ----------------
    clmopList <- mcRFop_cls(x_df, y_fac, nTree = 100)
    
    # ---------------- save result (FIXED) ----------------
    
    outfile <- paste0("results/clmhot_opList_zone", i, ".csv")
    write.csv(as.data.frame(clmopList), outfile, row.names = FALSE)
    
    print(paste("[DONE]", colname))
    
    rm(dt_s, tab, clm_y, x_df, y_fac, clmopList)
    gc()
    
  }, error = function(e){
    
    print(paste("[ERROR]", colname, ":", conditionMessage(e)))
    gc()
    NULL
    
  })
  
}





# 8. One-Hot Random Forest zone 1 =======================================================

# rf params (adjust here for tuning) ---------
SMPL_POS   <- 8000L
SMPL_PA    <- 1/1.3
SMPL_MAXN  <- 20000L
BASE_SEED  <- 49L
VAR_ROW    <- 20L          # row in opList to pick variable set
NTREE1     <- 100L         # trees per forest in classOP optimization rounds
NTREE2     <- 500L         # trees per forest in classOP final run
NOP        <- 3L           # number of classOP optimization rounds

# zone column
colname <- paste0("zone", 1)

# Select varlist from optimization results
clmList  <- read.csv(paste0("results/clmhot_opList_zone", 1, ".csv"))
varlist  <- trimws(unlist(strsplit(clmList[VAR_ROW+1, 2], ",")))

# 8.1 Sampling ---------
xlist <- names(combined_data)[6:74]
cols  <- c(colname, xlist)

n1_all  <- sum(combined_data[[colname]] == 1, na.rm = TRUE)
if (is.na(n1_all) || n1_all == 0L) stop("No presence (1) rows for ", colname)

pos_i   <- min(SMPL_POS, as.integer(n1_all))
noise_i <- as.integer(round(0.10 * pos_i))
noise_i <- min(noise_i, pos_i - 1L)

dt_s <- smpl_pa(combined_data, colname,
                cols  = cols,
                pos   = pos_i,
                noise = noise_i,
                pa    = SMPL_PA,
                max_n = SMPL_MAXN,
                seed  = BASE_SEED)

dt_s[, .N, by = get(colname)]

# 8.2 Train RF (classOP) ---------
clm_y <- factor(dt_s[[colname]], levels = c(0, 1))
print(table(clm_y))

varlist <- intersect(varlist, names(dt_s))
if (length(varlist) < 2) stop("Too few predictors after intersect(varlist, names(dt_s))")

clm_x <- dt_s[, ..varlist]

clim_zOp <- classOP(clm_x, clm_y, nTree1 = NTREE1, nTree2 = NTREE2, nOP = NOP)

# 8.3 Save OP clm model ---------
dir.create("rf", showWarnings = FALSE, recursive = TRUE)
filename <- paste0("rf/clm_zOp_zone", 1, ".Rdata")
save(clim_zOp, file = filename)
rm(clim_zOp); gc()

# 8.4 Load + check + importance ---------
load(paste0("rf/clm_zOp_zone", 1, ".Rdata"))
print(clim_zOp)
clim_zOp$mtry
clim_zOp$ntree

importance_vals <- clim_zOp$importance
importance_df <- data.frame(
  Feature = rownames(importance_vals),
  Importance = importance_vals[, "MeanDecreaseGini"]
)
rm(clim_zOp); gc()




# 9. One-Hot Random Forest zone 2-55 =====================================================

for (i in 2:55){
  
  colname <- paste0("zone", i)
  
  # skip if no one-hot column
  if (!(colname %in% names(combined_data))){
    print(paste("[SKIP] column not found:", colname))
    next
  }
  
  # skip if no presence
  n1_all <- sum(combined_data[[colname]] == 1, na.rm = TRUE)
  if (is.na(n1_all) || n1_all == 0){
    print(paste("[SKIP] no 1s in:", colname))
    next
  }
  
  # skip if opList file missing
  opfile <- paste0("results/clmhot_opList_zone", i, ".csv")
  if (!file.exists(opfile)){
    print(paste("[SKIP] opList not found:", opfile))
    next
  }
  
  tryCatch({
    
    # 9.1 Select varlist ---------
    clmList <- read.csv(opfile)
    if (nrow(clmList) < VAR_ROW+1 || clmList[VAR_ROW+1, 2] == "0"){
      print(paste("[SKIP] no valid variable set at row", VAR_ROW, "for", colname))
      next
    }
    varlist <- trimws(unlist(strsplit(clmList[VAR_ROW+1, 2], ",")))
    
    # 9.2 Sampling ---------
    xlist <- names(combined_data)[6:74]
    cols  <- c(colname, xlist)
    
    pos_i   <- min(SMPL_POS, as.integer(n1_all))
    noise_i <- as.integer(round(0.10 * pos_i))
    noise_i <- max(1L, min(noise_i, pos_i - 1L))
    
    print(paste("[RUN]", colname,
                "n1=", n1_all,
                "pos=", pos_i,
                "vars=", length(varlist)))
    
    dt_s <- smpl_pa(combined_data, colname,
                    cols  = cols,
                    pos   = pos_i,
                    noise = noise_i,
                    pa    = SMPL_PA,
                    max_n = SMPL_MAXN,
                    seed  = BASE_SEED + i)
    
    tab <- dt_s[, .N, by = get(colname)]
    if (nrow(tab) < 2){
      print(paste("[SKIP] sampled data has only one class:", colname))
      next
    }
    
    # 9.3 Train RF (classOP) ---------
    clm_y <- factor(dt_s[[colname]], levels = c(0, 1))
    print(table(clm_y))
    
    varlist <- intersect(varlist, names(dt_s))
    if (length(varlist) < 2){
      print(paste("[SKIP] too few predictors for", colname))
      next
    }
    
    clm_x <- dt_s[, ..varlist]
    
    clim_zOp <- classOP(clm_x, clm_y, nTree1 = NTREE1, nTree2 = NTREE2, nOP = NOP)
    
    # 9.4 Save OP clm model ---------
    filename <- paste0("rf/clm_zOp_zone", i, ".Rdata")
    save(clim_zOp, file = filename)
    
    print(paste("[DONE]", colname))
    
    rm(dt_s, tab, clm_y, clm_x, clim_zOp, clmList, varlist)
    gc()
    
  }, error = function(e){
    print(paste("[ERROR]", colname, ":", conditionMessage(e)))
    gc()
    NULL
  })
  
}