computeAvgPopdens <- function(segdata,pop_j,pop_m,areamask,outdf,year){
  area_val <- sum(areamask@data@values,na.rm=T)
  pop_val_j <- sum(pop_j@data@values,na.rm=T)
  pop_val_m <- sum(pop_m@data@values,na.rm=T)
  
  pd_m <- pop_m@data@values / areamask@data@values
  pd_j <- pop_j@data@values / areamask@data@values
  pd_m <- sum(pd_m *areamask@data@values,na.rm = T)/sum(areamask@data@values,na.rm=T)
  pd_j <- sum(pd_j * areamask@data@values,na.rm = T)/sum(areamask@data@values,na.rm=T)
  
  temp <- data.frame(segID=segdata@data$segID,year=year, areaKm2=area_val,
                     pop_jones = pop_val_j, pop_merkens = pop_val_m, 
                     popdens_jones = pd_j,
                     popdens_merkens = pd_m, 
                     popdens_old = segdata@data$popdens_old,ciamID=segdata@data$ciamID)
  outdf <- rbind(outdf,temp)
  return(outdf)
}
compute_seg_v2 <- function(segname,years,df){
  seg_diva<-diva_lines[diva_lines$ciamID ==segname,]
  seg2 <- st_polygonize(seg_diva)
  seg2 <- as(seg2,"Spatial")
  plot(seg2)
  crs(seg2)<- "+proj=longlat"
  
  arval<- area(seg2)/1e6
  
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
    
    
    pop_val_j <- sum(jcrop@data@values,na.rm=T)
    pop_val_m <- sum(pcrop@data@values,na.rm=T)
    tmp1 <- data.frame(segID=seg2@data$segID,year=year, areaKm2=arval,
                       pop_jones = pop_val_j, pop_merkens = pop_val_m, 
                       popdens_jones = pop_val_j/arval,
                       popdens_merkens = pop_val_m/arval, 
                       popdens_old = seg2@data$popdens_old,ciamID=seg2@data$ciamID)
    
    df <- rbind(df,tmp1)
    
  }
  return(df)
}
compute_seg_v3 <- function(segname,years,df){
  seg_buffer<-buffer_nonpoly[buffer_nonpoly$ciamID ==segname,]
  seg2 <- as(seg_buffer,"Spatial")
  plot(seg2)
  crs(seg2)<- "+proj=longlat +ellps=WGS84"
  
  arval<- .5*area(seg2)/1e6
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
    
    
    pop_val_j <- sum(jcrop@data@values,na.rm=T)
    pop_val_m <- sum(pcrop@data@values,na.rm=T)
    tmp1 <- data.frame(segID=seg2@data$segID,year=year, areaKm2=arval,
                       pop_jones = pop_val_j, pop_merkens = pop_val_m, 
                       popdens_jones = pop_val_j/arval,
                       popdens_merkens = pop_val_m/arval, 
                       popdens_old = seg2@data$popdens_old,ciamID=seg2@data$ciamID)
    
    df <- rbind(df,tmp1)
    
  }
  return(df)
}
