install.packages("H:/Jing/CEMT.zip", repos = NULL, type = "source", dependencies = TRUE)
library(CEMT)
library(terra)
# library(rgdal)  # rgdal is retired; terra now handles CRS operations

setwd("H:/Jing/ecoChina2")

## 1. convert veg type data to table---------
r <- rast('raster/veg_3');r 
#r2 <- rast('H:/Jing/ecoChina2/raster/veg_chn_poly2ras/veg_chn_poly2ras.tif');r2 #new raster I downloaded: 1:1million
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 


writeRaster(r,'raster/ecosys_ori.tif',filetype="GTiff",overwrite=TRUE)

xyv <- as.data.frame(r, xy = TRUE, na.rm = TRUE);hd(xyv)
#Q:x = longitude, y = latitude, v = ?some row has v value == ecosys name, but others are empty

id <- as.data.frame(r, cells = TRUE, xy = TRUE, na.rm = TRUE);hd(id)
print(sort(unique(id$veg_3))) #levels = 1:56

xyv1 <- id
names(xyv1)[3:4] <- c('zoneID','zone');hd(xyv1)
xyv1$zoneID <- as.factor(xyv1$zoneID);str(xyv1)
xyv1 <- xyv1[complete.cases(xyv1[, 1:3]), ];hd(xyv1) #nothing deleted
#revision: not remoging "" zone #xyv2 <- droplevels(xyv1[xyv1$zone!="",]);hd(xyv2) #remove "" zone

#zone and zoneID ---
lvl <- data.frame(cats(r)[[1]]);lvl
# revision: no removing lvl2 <- droplevels(lvl[lvl$VEGETATI_3!="",]);hd(lvl2) #remove "" zone
fWrite(lvl,'data_raw/1. zoneID_zone_count.csv')

rm(lvl,xyv,id,r2);gc()


# get DEM----------------
dem <- CEMT::getDEM('chn',id[,c('x','y')],z=90) #this function does not work on my laptop, need to retrieve from 3333b

id$id2 <- NA;hd(id)
colnames(id) <- c("lon", "lat", "Vegetati_3", "id2");hd(id)
id <- id[, c("Vegetati_3", "id2", "lat", "lon")];hd(id)

dem <- cbind(id,dem);hd(dem)

#output ClimateAP input file----
crd <- data.frame(ID=row(dt)[,1],dt[,c(3,2,1)],dem);hd(crd)
fWrite(dem,'1. coord.csv')

fRead('1. coord.csv')
