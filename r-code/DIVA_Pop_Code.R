## ---- Description ----------------------------------
# Computes segment-level coastal population at DIVA resolution for use in the CIAM model. 
# DIVA shapefiles are available from 
# https://gitlab.com/daniel.lincke.globalclimateforum.org/diva_published/-/tree/master/data/gis
# Gridded coastal population is from Jones and O'Neill (2016), downscaled to 1 km by Gao (2020)
# and available to download at 
# https://sedac.ciesin.columbia.edu/data/set/popdynamics-1-km-downscaled-pop-base-year-projection-ssp-2000-2100-rev01/
# Second gridded coastal population dataset is from Merkens et al (2016), from 
# https://figshare.com/s/9a94ae958d6a45684382
# NOTE 1/31/21: I ended up using Jones and O'Neill (2016) for population data in my runs. 
#   If we end up using Merkens, we should acknowledge their help in obtaining population data in paper
# ----------------------------------------------------

library(sf)
library(ncdf4)
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(rgdal)
library(tidyverse)
library(rnaturalearth)

# set which SSP you want (1-5)
SSP <- 5

setwd("/Users/aewsma/codes/CIAM_adaptation_regimes/R-code")
source("proc_Diva_pop_data.R")

# Read CIAM/DIVA data 
xsc <- read_csv("/Users/aewsma/.julia/dev/MimiCIAM/data/input/xsc.csv")
isl <- xsc$seg[xsc$island==1]
diva_db <- read_csv("/Users/aewsma/.julia/dev/MimiCIAM/data/input/data.csv")%>%group_by(X1)%>%
  mutate(area=sum(area1,area2,area3,area4,area5,area6,area7,area8,area9,area10,area11,area12,area13,area14,area15))
diva_latlon <- read_csv("/Users/aewsma/.julia/dev/MimiCIAM/data/diva_segment_latlon.csv")

diva_lines <- read_sf("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/COASTLINES/diva_published-master-data-gis/data/gis/cls.shp")
diva_lines <- diva_lines %>% mutate(segNum=gsub("[^0-9.]","",locationna),segStr=gsub(" ","",gsub("[0-9.]","",locationna)),
                                    segID=paste0(segStr,segNum), popdens_old = diva_db$popdens[match(segID,diva_db$X1)])
# Fix names 
ciam_Names <- xsc$seg
diva_names <- unique(diva_lines$segID)
excl_df <- diva_lines %>% 
  mutate(matchCIAM = diva_latlon$segments[match(paste0(round(longi,3),round(lati,3)),
                                                paste0(diva_latlon$longi,diva_latlon$lati))])

diva_lines <- diva_lines %>% mutate(ciamID=excl_df$matchCIAM[match(segID,excl_df$segID)])
diva_lines <- diva_lines %>% mutate(popdens_old = diva_db$popdens[match(ciamID,diva_db$X1)])
diva_lines <- diva_lines %>% mutate(ciamID=ifelse(segID=="Indonesia3287","Indonesia3288",
                                                  ifelse(segID=="Micronesia,Fe11797","MicronesiaFedStates11798",
                                                         ifelse(segID=="Micronesia,Fe11799","MicronesiaFedStates11800",
                                                                ifelse(segID=="FrenchPolynesi5866","FrenchPolynesia5867",
                                                                       ifelse(segID=="FrenchPolynesi6924","FrenchPolynesia6925",ciamID))))))
# "Indonesia3287","Indonesia3288",ciamID))
# "MicronesiaFedStates11797","MicronesiaFedStates11798",ciamID))
# "Micronesia,Fe11799","MicronesiaFedStates11800",ciamID))
# "FrenchPolynesi5866","FrenchPolynesia5867",ciamID))

# Read Area Data
area_data <- raster("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/GPW-LANDAREA-30/gpw-v4-land-water-area-rev11_landareakm_30_sec_tif/gpw_v4_land_water_area_rev11_landareakm_30_sec.tif")

# Read Merkens Data 
ssp_2060m <- raster(paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/MERKENS_2016/DIVA_SSP",SSP,"_2060.tif"))
ssp_2070m <- raster(paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/MERKENS_2016/DIVA_SSP",SSP,"_2070.tif"))
ssp_2080m <- raster(paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/MERKENS_2016/DIVA_SSP",SSP,"_2080.tif"))
ssp_2090m <- raster(paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/MERKENS_2016/DIVA_SSP",SSP,"_2090.tif"))
ssp_2100m <- raster(paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/MERKENS_2016/DIVA_SSP",SSP,"_2100.tif"))

ssp_2010m <- raster(paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/MERKENS_2016/DIVA_SSP",SSP,"_2010.tif"))
ssp_2020m <- raster(paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/MERKENS_2016/DIVA_SSP",SSP,"_2020.tif"))
ssp_2030m <- raster(paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/MERKENS_2016/DIVA_SSP",SSP,"_2030.tif"))
ssp_2040m <- raster(paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/MERKENS_2016/DIVA_SSP",SSP,"_2040.tif"))
ssp_2050m <- raster(paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/MERKENS_2016/DIVA_SSP",SSP,"_2050.tif"))


# Jones Data: Read coastal population data as raster. Unit: # of persons  
ncname <- paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/CIAM-GIS/POPDATA/popdynamics-1-km-downscaled-pop-base-year-projection-ssp-2000-2100-rev01-proj-ssp",SSP,"-netcdf/SSP",SSP,"_1km/ssp",SSP,"_total_2060.nc4")
ssp_2060j <- raster(paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/popdynamics-1-km-downscaled-pop-base-year-projection-ssp-2000-2100-rev01-proj-ssp",SSP,"-netcdf/SSP",SSP,"_1km/ssp",SSP,"_total_2060.nc4"),varname="Band1")
ssp_2070j <- raster(paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/popdynamics-1-km-downscaled-pop-base-year-projection-ssp-2000-2100-rev01-proj-ssp",SSP,"-netcdf/SSP",SSP,"_1km/ssp",SSP,"_total_2070.nc4"),varname="Band1")
ssp_2080j <- raster(paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/popdynamics-1-km-downscaled-pop-base-year-projection-ssp-2000-2100-rev01-proj-ssp",SSP,"-netcdf/SSP",SSP,"_1km/ssp",SSP,"_total_2080.nc4"),varname="Band1")
ssp_2090j <- raster(paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/popdynamics-1-km-downscaled-pop-base-year-projection-ssp-2000-2100-rev01-proj-ssp",SSP,"-netcdf/SSP",SSP,"_1km/ssp",SSP,"_total_2090.nc4"),varname="Band1")
ssp_2100j <- raster(paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/popdynamics-1-km-downscaled-pop-base-year-projection-ssp-2000-2100-rev01-proj-ssp",SSP,"-netcdf/SSP",SSP,"_1km/ssp",SSP,"_total_2100.nc4"),varname="Band1")

ssp_2010j <- raster(paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/popdynamics-1-km-downscaled-pop-base-year-projection-ssp-2000-2100-rev01-proj-ssp",SSP,"-netcdf/SSP",SSP,"_1km/ssp",SSP,"_total_2010.nc4"),varname="Band1")
ssp_2020j <- raster(paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/popdynamics-1-km-downscaled-pop-base-year-projection-ssp-2000-2100-rev01-proj-ssp",SSP,"-netcdf/SSP",SSP,"_1km/ssp",SSP,"_total_2020.nc4"),varname="Band1")
ssp_2030j <- raster(paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/popdynamics-1-km-downscaled-pop-base-year-projection-ssp-2000-2100-rev01-proj-ssp",SSP,"-netcdf/SSP",SSP,"_1km/ssp",SSP,"_total_2030.nc4"),varname="Band1")
ssp_2040j <- raster(paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/popdynamics-1-km-downscaled-pop-base-year-projection-ssp-2000-2100-rev01-proj-ssp",SSP,"-netcdf/SSP",SSP,"_1km/ssp",SSP,"_total_2040.nc4"),varname="Band1")
ssp_2050j <- raster(paste0("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/POPDATA/popdynamics-1-km-downscaled-pop-base-year-projection-ssp-2000-2100-rev01-proj-ssp",SSP,"-netcdf/SSP",SSP,"_1km/ssp",SSP,"_total_2050.nc4"),varname="Band1")


# Read coastal buffer data 
buffer <- read_sf("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/buffer/buffer2.shp")
buffer <- buffer %>% mutate(segNum=gsub("[^0-9.]","",locationna),segStr=gsub(" ","",gsub("[0-9.]","",locationna)),
segID=paste0(segStr,segNum), popdens_old = diva_db$popdens[match(segID,diva_db$X1)])

# Harmonize issues with CIAM Names
ciam_Names <- diva_db$X1
buffer_names <- buffer$segID
excl_df <- buffer %>% 
  mutate(matchCIAM = diva_latlon$segments[match(paste0(round(longi,3),round(lati,3)),
                                                paste0(diva_latlon$longi,diva_latlon$lati))])

buffer <- buffer %>% mutate(ciamID=excl_df$matchCIAM[match(segID,excl_df$segID)])
buffer <- buffer %>% mutate(popdens_old = diva_db$popdens[match(ciamID,diva_db$X1)])
buffer<- buffer %>% mutate(ciamID=ifelse(segID=="Indonesia3287","Indonesia3288",
                                    ifelse(segID=="Micronesia,Fe11797","MicronesiaFedStates11798",
                                           ifelse(segID=="Micronesia,Fe11799","MicronesiaFedStates11800",
                                                  ifelse(segID=="FrenchPolynesi5866","FrenchPolynesia5867",
                                                         ifelse(segID=="FrenchPolynesi6924","FrenchPolynesia6925",ciamID))))))

## what follows should really be a function over general time period `years`...

##
## first for 2010-2050
##

years <- c(2010,2020,2030,2040,2050)
outname <- paste("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/output_ssp",SSP,"_",min(years),"_",max(years),".csv", sep="")

## Step 1: Isolate DIVA Line Segments that are already polygons
outdf <- data.frame()
for (i in ciam_Names){
  seg <- diva_lines[diva_lines$ciamID==i,]
  
  try( {seg2 <- st_polygonize(seg)
    seg2 <- as(seg2,"Spatial")
    plot(seg2)
  
    crs(seg2)<- "+proj=longlat"
    #seg3 <- spTransform(seg2,"+proj=longlat +ellps=WGS84 +no_defs")
  
    area_seg <- crop(area_data,seg2)
    area_seg <- mask(area_seg,seg2)
  
    for (y in 1:length(years)){
      pstr <- paste0("ssp_",years[y],"m")
      jstr <- paste0("ssp_",years[y],"j")
      year <- years[y]
      pvar <- get(pstr)
      jvar=get(jstr)
      
      pcrop <- crop(pvar,seg2)
      jcrop <- crop(jvar,seg2)
      
      pcrop <- mask(pcrop,seg2)
      jcrop <- mask(jcrop,seg2)
      
      
      if (sum(area_seg@data@values,na.rm = T)==0){}
      else{outdf <- computeAvgPopdens(seg2,jcrop,pcrop,area_seg,outdf,year)}
      
    } # end `for`

  }) # end `try`
}


# Phase 2: Use Buffer data for non-polygon segments 
poly_names <- unique(outdf$ciamID)
nonpoly_names <- ciam_Names[!(ciam_Names %in% poly_names)]
buffer_nonpoly <- buffer %>% filter(ciamID %in% nonpoly_names)

for (i in unique(nonpoly_names)){
  seg<- buffer_nonpoly[buffer_nonpoly$ciamID==i,]
  
  try( {
  
    seg2 <- as(seg,"Spatial")
    plot(seg2)
    crs(seg2)<- "+proj=longlat"
    #seg3 <- spTransform(seg2,"+proj=longlat +ellps=WGS84 +no_defs")
    
    area_seg <- crop(area_data,seg2)
    area_seg <- mask(area_seg,seg2)
    for (y in 1:length(years)){
      pstr <- paste0("ssp_",years[y],"m")
      jstr <- paste0("ssp_",years[y],"j")
      year <- years[y]
      pvar <- get(pstr)
      jvar=get(jstr)
      
      pcrop <- crop(pvar,seg2)
      jcrop <- crop(jvar,seg2)
      
      pcrop <- mask(pcrop,seg2)
      jcrop <- mask(jcrop,seg2)
      
      
      if (sum(area_seg@data@values,na.rm = T)==0){}
      else{outdf <- computeAvgPopdens(seg2,jcrop,pcrop,area_seg,outdf,year)}
      
    } # end `for`
 
  }) # end `try`
    
}



# ID duplicate segments
dupes <- c()
for (i in ciam_Names) {
  for (year in years) { # TW: added - okay?
    tmp <- outdf %>% filter(year==year,ciamID==i) # TW: changed from 2010 (or 2060; first year in years) - okay?
    if (nrow(tmp)>1){
      dupes <- c(dupes,i)
    }
  }
}

dupe_segs <- outdf %>% filter(ciamID %in% dupes)
ok_segs <- unique(outdf$ciamID)[!(unique(outdf$ciamID)%in% dupes)]

### Redo Duplicate Segments 
outdf2 <- data.frame()
for (i in dupes){
  seg <- diva_lines[diva_lines$ciamID==i,]
  
  try( {seg2 <- st_polygonize(seg)
  seg2 <- as(seg2,"Spatial")
  plot(seg2)
  
  crs(seg2)<- "+proj=longlat"
  #seg3 <- spTransform(seg2,"+proj=longlat +ellps=WGS84 +no_defs")
  
  area_seg <- crop(area_data,seg2)
  area_seg <- mask(area_seg,seg2)
  
  p_m10 <- crop(ssp_2010m, seg2)
  p_j10 <- crop(ssp_2010j,seg2)
  p_m20 <- crop(ssp_2020m, seg2)
  p_j20 <- crop(ssp_2020j,seg2)
  p_m20 <- mask(p_m20,seg2)
  p_j20 <- mask(p_j20,seg2)
  p_m10 <- mask(p_m10,seg2)
  p_j10 <- mask(p_j10,seg2)
  
  p_m30<- crop(ssp_2030m, seg2)
  p_j30<- crop(ssp_2030j, seg2)
  p_m30<- mask(p_m30,seg2)
  p_j30<- mask(p_j30,seg2)
  
  p_m40<- crop(ssp_2040m, seg2)
  p_j40<- crop(ssp_2040j, seg2)
  p_m40<- mask(p_m40,seg2)
  p_j40<- mask(p_j40,seg2)
  
  p_m50<- crop(ssp_2050m, seg2)
  p_j50<- crop(ssp_2050j, seg2)
  p_m50<- mask(p_m50,seg2)
  p_j50<- mask(p_j50,seg2)
  
  jvars <- c("p_j10","p_j20","p_j30","p_j40","p_j50")
  pvars <- c("p_m10","p_m20","p_m30","p_m40","p_m50")
  yrs <- c(2010,2020,2030,2040,2050)
  for (y in 1:length(yrs)){
    year <- yrs[y]
    pvar <- get(pvars[y])
    jvar=get(jvars[y])
    
    
    if (sum(area_seg@data@values,na.rm = T)==0){}
    else{outdf2 <- computeAvgPopdens(seg2,jvar,pvar,area_seg,outdf2,year)}
    
  }
  })
  
}

redupe_segs <- unique(outdf2$ciamID)
excl3 <- dupes[!(dupes %in% redupe_segs)]

outdf_10_50_new <- outdf %>% filter(!(ciamID %in% dupes))
outdf_10_50_new <- rbind(outdf_10_50_new,outdf2)

psegs <- ciam_Names[!(ciam_Names%in%unique(outdf_10_50_new$ciamID))]
prob_segs <- c(psegs)

outdf_10_50_new2 <- outdf_10_50_new %>% filter(!(ciamID %in% prob_segs))

tmp <- data.frame()

### Other Issue Segs
iss_segs <- c("MarshallIslands11966","Maldives9594","Maldives9715","MarshallIslands11916","Fiji5688","MarshallIslands11964","MarshallIslands11884","MarshallIslands11910","Maldives9592")
iss_df <- outdf %>% filter(popdens_jones>20,popdens_merkens==0)
iss_df2 <- outdf %>% filter(popdens_jones==0,popdens_merkens>20)
iss_df3 <- outdf %>% filter(popdens_old==0,popdens_merkens>20)
iss_df4 <- outdf %>% filter(popdens_old==0,popdens_jones>20)
iss_segs <- c(iss_segs,unique(iss_df$ciamID),unique(iss_df2$ciamID),unique(iss_df3$ciamID),unique(iss_df4$ciamID))%>% unique()

for (i in iss_segs){
  
  try({ tmp <- compute_seg_v2(i,years,tmp)})
}

iss_segs2 <- iss_segs[!(iss_segs %in% unique(tmp$ciamID))]
for (i in iss_segs2){
  try({ tmp <- compute_seg_v3(i,years,tmp)})
}

outdf_new <- outdf %>% filter(!(ciamID %in% iss_segs))
outdf_new <- rbind(outdf_new,tmp)

# Try again with problem segments
prob_segs <- c(ciam_Names[!(ciam_Names %in% outdf_new$ciamID)],"UnitedKingdom8057")%>% unique()
diva_probs <- diva_lines[diva_lines$ciamID %in% prob_segs,]
#st_write(diva_probs,"Desktop/CIAM-GIS/divaProbs/divaProbs3.shp","divaProbs3.shp")

## Some prob segs are polygons; id these
##  These tend to be areas where area data does not think there is land or too small for area data
prob_polys <- c()
tmp3 <- data.frame()
for (i in prob_segs){
  seg <- diva_lines[diva_lines$ciamID==i,]
  
  try( {seg2 <- st_polygonize(seg)
  seg2 <- as(seg2,"Spatial")
  plot(seg2)
  
  crs(seg2)<- "+proj=longlat"
  prob_polys <- c(prob_polys,i)
  tmp3 <- compute_seg_v2(i,years,tmp3)
  })
}
prob_nonpolys <- prob_segs[!(prob_segs %in% prob_polys)]
buffer_prob <- read_sf("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/divaProbs/buffer_probs3.shp")
buffer_prob <- buffer_prob %>% filter(ciamID %in%prob_nonpolys)
for (j in prob_nonpolys){
  tmp3 <- compute_seg_v3(j,years,tmp3)
}

outdf_fin <- outdf_new %>% filter(!(ciamID %in% unique(tmp3$ciamID)))
outdf_fin <- rbind(outdf_fin,tmp3)
i<-ciam_Names[!(ciam_Names %in% outdf_fin$ciamID)]

seg <- diva_lines[diva_lines$ciamID=="Chile8809",]
seg2 <- st_polygonize(seg)
seg2<- as(seg2,"Spatial")
crs(seg2)<-"+proj=longlat"

area_seg <- crop(area_data,seg2)
area_seg <- mask(area_seg,seg2)
areaval <- sum(area_seg@data@values,na.rm=T)

p10m <- crop(ssp_2010m,seg2)
p10m <- mask(p10m,seg2)
p20m <- crop(ssp_2020m,seg2)
p20m <- mask(p20m,seg2)
p30m <- crop(ssp_2030m,seg2)
p30m <- mask(p30m,seg2)
p40m <- crop(ssp_2040m,seg2)
p40m <- mask(p40m,seg2)
p50m <- crop(ssp_2050m,seg2)
p50m <- mask(p50m,seg2)

pop10<-sum(p10m@data@values,na.rm=T)
pop20<-sum(p20m@data@values,na.rm=T)
pop30<-sum(p30m@data@values,na.rm=T)
pop40<-sum(p40m@data@values,na.rm=T)
pop50<-sum(p50m@data@values,na.rm=T)

for (y in years){
  tmp1 <- data.frame(segID=seg2@data$segID,year=y, areaKm2=areaval,
                     pop_jones = 0, pop_merkens = 0, 
                     popdens_jones = 0,
                     popdens_merkens = 0, 
                     popdens_old = seg2@data$popdens_old,ciamID=seg2@data$ciamID)
  outdf_fin<- rbind(outdf_fin,tmp1)
}

write_csv(outdf_fin,outname)

##
## same but for 2060-2100
##

years <- c(2060,2070,2080,2090,2100)
outname <- paste("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/output_ssp",SSP,"_",min(years),"_",max(years),".csv", sep="")
## Step 1: Isolate DIVA Line Segments that are already polygons
outdf <- data.frame()
for (i in ciam_Names){
  seg <- diva_lines[diva_lines$ciamID==i,]
  
  try( {seg2 <- st_polygonize(seg)
  seg2 <- as(seg2,"Spatial")
  plot(seg2)
  
  crs(seg2)<- "+proj=longlat"
  #seg3 <- spTransform(seg2,"+proj=longlat +ellps=WGS84 +no_defs")
  
  area_seg <- crop(area_data,seg2)
  area_seg <- mask(area_seg,seg2)
  
  for (y in 1:length(years)){
    pstr <- paste0("ssp_",years[y],"m")
    jstr <- paste0("ssp_",years[y],"j")
    year <- years[y]
    pvar <- get(pstr)
    jvar=get(jstr)
    
    pcrop <- crop(pvar,seg2)
    jcrop <- crop(jvar,seg2)
    
    pcrop <- mask(pcrop,seg2)
    jcrop <- mask(jcrop,seg2)
    
    
    if (sum(area_seg@data@values,na.rm = T)==0){}
    else{outdf <- computeAvgPopdens(seg2,jcrop,pcrop,area_seg,outdf,year)}
    
  } # end `for`
  
  }) # end `try`
}


# Phase 2: Use Buffer data for non-polygon segments 
poly_names <- unique(outdf$ciamID)
nonpoly_names <- ciam_Names[!(ciam_Names %in% poly_names)]
buffer_nonpoly <- buffer %>% filter(ciamID %in% nonpoly_names)

for (i in unique(nonpoly_names)){
  seg<- buffer_nonpoly[buffer_nonpoly$ciamID==i,]
  
  try( {
    
    seg2 <- as(seg,"Spatial")
    plot(seg2)
    crs(seg2)<- "+proj=longlat"
    #seg3 <- spTransform(seg2,"+proj=longlat +ellps=WGS84 +no_defs")
    
    area_seg <- crop(area_data,seg2)
    area_seg <- mask(area_seg,seg2)
    for (y in 1:length(years)){
      pstr <- paste0("ssp_",years[y],"m")
      jstr <- paste0("ssp_",years[y],"j")
      year <- years[y]
      pvar <- get(pstr)
      jvar=get(jstr)
      
      pcrop <- crop(pvar,seg2)
      jcrop <- crop(jvar,seg2)
      
      pcrop <- mask(pcrop,seg2)
      jcrop <- mask(jcrop,seg2)
      
      if (sum(area_seg@data@values,na.rm = T)==0){}
      else{outdf <- computeAvgPopdens(seg2,jcrop,pcrop,area_seg,outdf,year)}
      
    } # end `for`
    
  }) # end `try`
  
}



# ID duplicate segments
dupes <- c()
for (i in ciam_Names) {
  for (year in years) {
    tmp <- outdf %>% filter(year==year,ciamID==i)
    if (nrow(tmp)>1){
      dupes <- c(dupes,i)
    }
  }
}

dupe_segs <- outdf %>% filter(ciamID %in% dupes)
ok_segs <- unique(outdf$ciamID)[!(unique(outdf$ciamID)%in% dupes)]

### Redo Duplicate Segments 
outdf2 <- data.frame()
for (i in dupes){
  seg <- diva_lines[diva_lines$ciamID==i,]
  
  try( {seg2 <- st_polygonize(seg)
  seg2 <- as(seg2,"Spatial")
  plot(seg2)
  
  crs(seg2)<- "+proj=longlat"
  #seg3 <- spTransform(seg2,"+proj=longlat +ellps=WGS84 +no_defs")
  
  area_seg <- crop(area_data,seg2)
  area_seg <- mask(area_seg,seg2)
  
  p_m60 <- crop(ssp_2060m, seg2)
  p_j60 <- crop(ssp_2060j,seg2)
  p_m70 <- crop(ssp_2070m, seg2)
  p_j70 <- crop(ssp_2070j,seg2)
  p_m70 <- mask(p_m70,seg2)
  p_j70 <- mask(p_j70,seg2)
  p_m60 <- mask(p_m60,seg2)
  p_j60 <- mask(p_j60,seg2)
  
  p_m80<- crop(ssp_2080m, seg2)
  p_j80<- crop(ssp_2080j, seg2)
  p_m80<- mask(p_m80,seg2)
  p_j80<- mask(p_j80,seg2)
  
  p_m90<- crop(ssp_2090m, seg2)
  p_j90<- crop(ssp_2090j, seg2)
  p_m90<- mask(p_m90,seg2)
  p_j90<- mask(p_j90,seg2)
  
  p_m100<- crop(ssp_2100m, seg2)
  p_j100<- crop(ssp_2100j, seg2)
  p_m100<- mask(p_m100,seg2)
  p_j100<- mask(p_j100,seg2)
  
  jvars <- c("p_j60","p_j70","p_j80","p_j90","p_j100")
  pvars <- c("p_m60","p_m70","p_m80","p_m90","p_m100")
  yrs <- c(2060,2070,2080,2090,2100)
  for (y in 1:length(yrs)){
    year <- yrs[y]
    pvar <- get(pvars[y])
    jvar=get(jvars[y])
    
    if (sum(area_seg@data@values,na.rm = T)==0){}
    else{outdf2 <- computeAvgPopdens(seg2,jvar,pvar,area_seg,outdf2,year)}
    
  }
  })
  
}

redupe_segs <- unique(outdf2$ciamID)
excl3 <- dupes[!(dupes %in% redupe_segs)]

outdf_60_100_new <- outdf %>% filter(!(ciamID %in% dupes))
outdf_60_100_new <- rbind(outdf_60_100_new,outdf2)

psegs <- ciam_Names[!(ciam_Names%in%unique(outdf_60_100_new$ciamID))]
prob_segs <- c(psegs)

outdf_60_100_new2 <- outdf_60_100_new %>% filter(!(ciamID %in% prob_segs))

# this function doesn't need to be defined again - put in a separate file to source
#DON'T RUN THIS YET UNTIL YOU SORT OUT THE YEARS HARDCODING ISSUE IN THERE
#compute_seg_v2 <- function(segname,years,df){

tmp <- data.frame()

### Other Issue Segs
iss_segs <- c("MarshallIslands11966","Maldives9594","Maldives9715","MarshallIslands11916","Fiji5688","MarshallIslands11964","MarshallIslands11884","MarshallIslands11910","Maldives9592")
iss_df <- outdf %>% filter(popdens_jones>20,popdens_merkens==0)
iss_df2 <- outdf %>% filter(popdens_jones==0,popdens_merkens>20)
iss_df3 <- outdf %>% filter(popdens_old==0,popdens_merkens>20)
iss_df4 <- outdf %>% filter(popdens_old==0,popdens_jones>20)
iss_segs <- c(iss_segs,unique(iss_df$ciamID),unique(iss_df2$ciamID),unique(iss_df3$ciamID),unique(iss_df4$ciamID))%>% unique()

for (i in iss_segs){
  try({ tmp <- compute_seg_v2(i,years,tmp)})
  
}
iss_segs2 <- iss_segs[!(iss_segs %in% unique(tmp$ciamID))]
for (i in iss_segs2){
  try({ tmp <- compute_seg_v3(i,years,tmp)})
  
}

outdf_new <- outdf %>% filter(!(ciamID %in% iss_segs))
outdf_new <- rbind(outdf_new,tmp)

# Try again with problem segments
prob_segs <- c(ciam_Names[!(ciam_Names %in% outdf_new$ciamID)],"UnitedKingdom8057")%>% unique()
diva_probs <- diva_lines[diva_lines$ciamID %in% prob_segs,]
#st_write(diva_probs,"Desktop/CIAM-GIS/divaProbs/divaProbs3.shp","divaProbs3.shp")

## Some prob segs are polygons; id these
##  These tend to be areas where area data does not think there is land or too small for area data
prob_polys <- c()
tmp3 <- data.frame()
for (i in prob_segs){
  seg <- diva_lines[diva_lines$ciamID==i,]
  
  try( {seg2 <- st_polygonize(seg)
  seg2 <- as(seg2,"Spatial")
  plot(seg2)
  
  crs(seg2)<- "+proj=longlat"
  prob_polys <- c(prob_polys,i)
  tmp3 <- compute_seg_v2(i,years,tmp3)
  })
}
prob_nonpolys <- prob_segs[!(prob_segs %in% prob_polys)]
buffer_prob <- read_sf("/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/divaProbs/buffer_probs3.shp")
buffer_prob <- buffer_prob %>% filter(ciamID %in%prob_nonpolys)
for (j in prob_nonpolys){
  tmp3 <- compute_seg_v3(j,years,tmp3)
}

outdf_fin <- outdf_new %>% filter(!(ciamID %in% unique(tmp3$ciamID)))
outdf_fin <- rbind(outdf_fin,tmp3)
i<-ciam_Names[!(ciam_Names %in% outdf_fin$ciamID)]

seg <- diva_lines[diva_lines$ciamID=="Chile8809",]
seg2 <- st_polygonize(seg)
seg2<- as(seg2,"Spatial")
crs(seg2)<-"+proj=longlat"

area_seg <- crop(area_data,seg2)
area_seg <- mask(area_seg,seg2)
areaval <- sum(area_seg@data@values,na.rm=T)

p60m <- crop(ssp_2060m,seg2)
p60m <- mask(p60m,seg2)
p70m <- crop(ssp_2070m,seg2)
p70m <- mask(p70m,seg2)
p80m <- crop(ssp_2080m,seg2)
p80m <- mask(p80m,seg2)
p90m <- crop(ssp_2090m,seg2)
p90m <- mask(p90m,seg2)
p100m <- crop(ssp_2100m,seg2)
p100m <- mask(p100m,seg2)

pop60<-sum(p60m@data@values,na.rm=T)
pop70<-sum(p70m@data@values,na.rm=T)
pop80<-sum(p80m@data@values,na.rm=T)
pop90<-sum(p90m@data@values,na.rm=T)
pop100<-sum(p100m@data@values,na.rm=T)

for (y in years){
  tmp1 <- data.frame(segID=seg2@data$segID,year=y, areaKm2=areaval,
                     pop_jones = 0, pop_merkens = 0, 
                     popdens_jones = 0,
                     popdens_merkens = 0, 
                     popdens_old = seg2@data$popdens_old,ciamID=seg2@data$ciamID)
  outdf_fin<- rbind(outdf_fin,tmp1)
}

write_csv(outdf_fin,outname)
