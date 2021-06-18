### Compare results of coastal population density scenarios 
### and create forcing data files for CIAM
library(tidyverse)

### things to set up
ciampath <- "/Users/aewsma/.julia/dev/MimiCIAM/"
gispath <- "/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/"
SSP<- 5



# Read in Data 
pjm1 <- read_csv(paste0(gispath,"output_ssp",SSP,"_2010_2050.csv"))
pjm2 <- read_csv(paste0(gispath,"output_ssp",SSP,"_2060_2100.csv"))
pjm <- rbind(pjm1,pjm2)

# Process New Data For Output - 12148 x 10 and duplicate last row of data for 20 years
library(reshape2)

# Jones and O'Neill data
pd_jones <- pjm %>% dplyr::select(ciamID,year,popdens_jones) 
pd_jones_last <- pd_jones %>% filter(year==2100) %>% mutate(year=2110)

newyrs=c(2120,2130,2140,2150,2160,2170,2180,2190,2200)
tmp <- data.frame()
outdf <- pd_jones_last
for (j in 1:length(newyrs)){
  tmp <- pd_jones_last %>% mutate(year=newyrs[j])
  outdf <- rbind(outdf,tmp)
}

pd_jones <- rbind(pd_jones,outdf)
pd_jones <- pd_jones %>% mutate(time=(year-2000)/10) %>% dplyr::select(-year) %>%rename(variable=time)%>%
  reshape2::dcast(variable~ciamID,value.var=c("popdens_jones"), fun.aggregate = max)

write_csv(pd_jones,paste0(ciampath,"data/ssp/popdens_seg_jones_ssp",SSP,".csv"))

# Merkens et al data
pd_m <- pjm %>% dplyr::select(ciamID,year,popdens_merkens) 
pd_m_last <- pd_m %>% filter(year==2100) %>% mutate(year=2110)

newyrs=c(2120,2130,2140,2150,2160,2170,2180,2190,2200)
tmp <- data.frame()
outdf <- pd_m_last
for (j in 1:length(newyrs)){
  tmp <- pd_m_last %>% mutate(year=newyrs[j])
  outdf <- rbind(outdf,tmp)
}

pd_m <- rbind(pd_m,outdf)
pd_m <- pd_m %>% mutate(time=(year-2000)/10) %>% dplyr::select(-year) %>%rename(variable=time)%>%
  reshape2::dcast(variable~ciamID,value.var=c("popdens_merkens"), fun.aggregate = max)

write_csv(pd_m,paste0(ciampath,"data/ssp/popdens_seg_merkens_ssp",SSP,".csv"))

if(FALSE) {
# Check against originals from Catherine
# commented out because not generally going to be needed, but testing workflow
  dfM_0 <- read_csv(paste0(ciampath,"data/ssp/popdens_seg_merkens_ssp",SSP,"_origCL.csv"))
  dfM_1 <- read_csv(paste0(ciampath,"data/ssp/popdens_seg_merkens_ssp",SSP,".csv"))
  dfM_diff <- dfM_1 - dfM_0
  max_diff <- rep(0, dim(dfJ_diff)[2]-1)
  for (i in 2:(length(max_diff)+1)) {max_diff[i] = max(dfM_1[1:5,i] - dfM_0[1:5,i])}
  dfJ_0 <- read_csv(paste0(ciampath,"data/ssp/popdens_seg_jones_ssp",SSP,"_origCL.csv"))
  dfJ_1 <- read_csv(paste0(ciampath,"data/ssp/popdens_seg_jones_ssp",SSP,".csv"))
  dfJ_diff <- dfJ_1 - dfJ_0
  max_diff <- rep(0, dim(dfJ_diff)[2]-1)
  for (i in 2:(length(max_diff)+1)) {max_diff[i] = max(dfJ_1[,i] - dfJ_0[,i])}
  # those indeed are a little bit different. many segments in 2060-2100, and a few in 2010-2050.
  pjm1_1 <- read_csv(paste0(gispath,"output_ssp",SSP,"_2010_2050.csv"))
  pjm1_1 <- pjm1_1[1:60740,]
  pjm2_1 <- read_csv(paste0(gispath,"output_ssp",SSP,"_2060_2100.csv"))
  pjm2_1 <- pjm2_1[1:60740,]
  pjm1_0 <- read_csv(paste0(gispath,"output_ssp",SSP,"_2010_2050_origCL.csv"))
  pjm2_0 <- read_csv(paste0(gispath,"output_ssp",SSP,"_2060_2100_origCL.csv"))
  for (seg in ciamIDs) {
    idx <- match(seg, ciamIDs)
    tmp1 <- pjm1_1 %>% filter(ciamID==seg); tmp1 <- tmp1["popdens_jones"]
    tmp0 <- pjm1_0 %>% filter(ciamID==seg); tmp0 <- tmp0["popdens_jones"]
    max_diff[idx] <- max(tmp1-tmp0)
  }
}