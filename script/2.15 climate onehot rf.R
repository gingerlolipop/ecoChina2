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

dat <- dat[complete.cases(dat$ZoneID), ]
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
SMPL_POS   <- 8000L      # target presence (cap by available n1)
SMPL_PA    <- 1/1.3       # P/A ratio (presence/absence)
SMPL_MAXN  <- 20000L      # max total sample size per zone
BASE_SEED  <- 49L



# 3.1 Climate var optimize, zone 1-55 ===================================================
for (i in 1:55){
  
  colname <- paste0("zone", i)
  
  #skip if no one-hot column
  if(!(colname %in% names(combined_data))){
    print(paste("[SKIP] column not found:", colname))
    next
  }
  
  #skip if no presence
  n1_all <- sum(combined_data[[colname]] == 1, na.rm = TRUE)
  if(is.na(n1_all) || n1_all == 0){
    print(paste("[SKIP] no 1s in:", colname))
    next
  }
  
  #noise = ~10% of target presence (small zones also have noise)
  pos_i   <- min(SMPL_POS, as.integer(n1_all))
  noise_i <- as.integer(round(0.10 * pos_i))
  noise_i <- min(noise_i, pos_i - 1L) #keep valid range
  
  #keep only y + climate predictors (faster)
  cols <- c(colname, names(combined_data)[6:74])
  
  print(paste("[RUN]", colname, "n1=", n1_all, "pos=", pos_i, "noise~", noise_i,
              "pa=", SMPL_PA, "max_n=", SMPL_MAXN))
  
  tryCatch({
    
    #run the smpl_pa function
    dt_s <- smpl_pa(combined_data, colname,
                    cols  = cols,
                    pos   = pos_i,
                    noise = noise_i,
                    pa    = SMPL_PA,
                    max_n = SMPL_MAXN,
                    seed  = BASE_SEED + i)
    
    #quick class check after sampling
    tab <- dt_s[, .N, by = get(colname)]
    if(nrow(tab) < 2){
      print(paste("[SKIP] sampled data has only one class:", colname))
      next
    }
    
    #y
    clm_y <- factor(dt_s[[colname]], levels = c(0, 1))
    print(table(clm_y))
    
    #x (use original mcRFop_cls: needs data.frame & no NA)
    x_df  <- na.omit(as.data.frame(dt_s[, ..xlist]))
    if(nrow(x_df) < 10){
      print(paste("[SKIP] too few rows after na.omit:", colname))
      next
    }
    y_fac <- factor(dt_s[[colname]][as.numeric(rownames(x_df))], levels=c(0,1))
    
    if(anyNA(y_fac) || length(unique(y_fac)) < 2){
      print(paste("[SKIP] y invalid / single class after NA removal:", colname))
      next
    }
    
    #mcRF optimize
    clmopList <- mcRFop_cls(x_df, y_fac, nTree=100)
    
    #save results
    fWrite(clmopList, 'results/clmhot_opList_zone', i, '.csv')
    
    print(paste("[DONE]", colname))
    
    rm(dt_s, tab, clm_y, x_df, y_fac, clmopList);gc()
    
  }, error = function(e){
    print(paste("[ERROR]", colname, ":", conditionMessage(e)))
    gc()
    NULL
  })
  
}