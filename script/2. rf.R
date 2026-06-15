
library(CEMT)
library(terra)
library(randomForest)

# Random Forest ============================================================
#rf data
getwd()
setwd("C:/Users/jillb/Documents/paper 2/China")
dat <- fRead('1. coordNormal_1961_1990SY.csv');hd(dat);names(dat)
summary(dat$ahm) #ahm max 35319
summary(dat$tave_jja) #also, the training Temperature data are 10 times it should be. Consistent to the raster values for prediction.
dat[dat == -9999] <- NA #I added this step to drop replace all -9999s, then re-trained the RF model and saved as rf.Rdata
dat <- dat[complete.cases(dat), ]

x <- dat[,c(6:29,38:61,64:69)]; hd(x)
summary(x$ahm)
y <- as.factor(dat$zoneID);levels(y)

#clm var optimize--
opList <- mcRFop_cls(x,y,nTree=100);opList
write.csv(opList,'rf/opList.csv',row.names=F)
opList <- read.csv('rf/opList.csv');opList
typeof(opList)
rfVar <- strToList(opList[16,2]);rfVar
typeof(rfVar)
x2 <- x[,rfVar];head(x2)
typeof(x2)

y2 <- y

# rf model
rf <- randomForest(x2,y2,ntree=200);print(rf) #accy=0.72 ; = 71.86 after re-running
save(rf, file='rf/rf.Rdata')
load('rf/rf.Rdata')

# Create confusion matrix
confusion_matrix <- rf$confusion
write.csv(confusion_matrix, file = "rf/climrf_conf_mat.csv")

#predict current=========================================================
setwd("C:/Users/jillb/Documents/paper 2/China")
clmDir <- "ClimateData/China4k/Normal_1961_1990SY/"
###stk <- rasterStack(clmDir,c('mat'),rType='grid');stk #added varlist mannualy

rlist <- list.files(clmDir,pattern='.tif',full.name=T); rlist[1:5]; #create a file list of climate variables
for(lyr in rlist){
  r <- raster(lyr)
#  if(grepl('MAT|MCMT|MWMT|TD|AHM|SHM|EMT|EXT|MAR',lyr)){r=r/10} #convert integer to decimal for temperature variables
  if(lyr==rlist[1]){stk6190=r} else{stk6190=stack(stk6190,r)} # “6190” representing 1961_1990
}; stk6190

library(randomForest)
rfVar <- rownames(randomForest::importance(rf));rfVar
clm <- subset(stk6190,rfVar);clm

load('rf/rf.Rdata')
p <- predict(clm,rf); plot(p)
crs(p) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

library(mapdata)
library(maps)
library(GISTools)
plot (p,
      main="Predicted Spatial Distribution of Ecosystems",
      sub = "Normal Period 1961-1990",
      xlab="Longitude", ylab="Latitude")

map('world',add=TRUE)
maps::map.scale(x = 112, y = 58, ratio = F,relwidth=0.16)
north.arrow(105,57,miles2ft(0.0001),col='grey')




p <- geoPrj(p);p
writeRaster(p,'raster/rf_p.tif',overwrite=T)

