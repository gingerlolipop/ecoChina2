install.packages("H:/Jing/CEMT.zip", repos = NULL, type = "source", dependencies = TRUE)
library(CEMT)
library(terra)
library(ClimateNAr)
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


# 2. get DEM----------------
dem <- fRead('data raw/1. coord.csv') #read dem, got before using the getDEM function from CEMT. now cannot use coz permission denied.

# 3. get Climate -----------
varList_Y=c("MAT","MWMT","MCMT","TD","MAP","MSP","AHM","SHM","bFFP","eFFP","FFP","CMD","CMI","DD_0","DD5","DD_18","DD18","DD1040","EMT","EXT",
            "Eref", "rsds","NFFD", "PAS","RH")
varList_S=c("Tmax_wt","Tmax_sp","Tmax_sm","Tmax_at","Tmin_wt","Tmin_sp","Tmin_sm","Tmin_at","Tave_wt","Tave_sp","Tave_sm","Tave_at",
            "PPT_wt","PPT_sp","PPT_sm","PPT_at","rsds_wt","rsds_sp","rsds_sm","rsds_at",
            "DD_0_wt","DD_0_sp","DD_0_sm","DD_0_at","DD5_wt","DD5_sp","DD5_sm","DD5_at","DD_18_wt",
            "DD_18_sp","DD_18_sm","DD_18_at","DD18_wt","DD18_sp","DD18_sm","DD18_at","NFFD_wt","NFFD_sp",
            "NFFD_sm","NFFD_at","PAS_wt","PAS_sp","PAS_sm","PAS_at","Eref_wt","Eref_sp","Eref_sm","Eref_at","CMD_wt",
            "CMD_sp","CMD_sm","CMD_at","RH_wt","RH_sp","RH_sm","RH_at","CMI_wt","CMI_sp","CMI_sm","CMI_at")


clm_vars <- c(varList_Y, varList_S); clm_vars

options(timeout = 300)
options(download.file.method = "libcurl")

clm_6190 <- ClimateNAr(dem, "Normal_1961_1990.nrm", clm_vars,
                       outDir = "raster/ClimateData/China4k/Normal_1961_1990SY/")
typeof(clm_6190)


