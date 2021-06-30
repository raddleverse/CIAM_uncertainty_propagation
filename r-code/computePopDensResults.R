### Compare results of coastal population density scenarios 
### and create forcing data files for CIAM
library(tidyverse)

### things to set up
ciampath <- "/Users/aewsma/.julia/dev/MimiCIAM/"
gispath <- "/Users/aewsma/codes/CIAM-Work/code/CIAM-GIS/"
SSP<- 5

# Read in Data 
## OLD (CL) VERSION:  (two separate files)
#pjm1 <- read_csv(paste0(gispath,"output_ssp",SSP,"_2010_2050.csv"))
#pjm2 <- read_csv(paste0(gispath,"output_ssp",SSP,"_2060_2100.csv"))
#pjm <- rbind(pjm1,pjm2)
## NEW (TW) VERSION:   (single large file)
pjm <- read_csv(paste0(gispath,"output_ssp",SSP,"_2010_2100.csv"))

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


## TW: below here is just a check against different versions of the DIVA_Pop_Code.R production, and originals from Catherine
if(FALSE) {
# Check against originals from Catherine
# commented out because not generally going to be needed, but testing workflow
  dfM_0 <- read_csv(paste0(ciampath,"data/ssp/popdens_seg_merkens_ssp",SSP,"_origCL.csv"))
  dfM_1 <- read_csv(paste0(ciampath,"data/ssp/popdens_seg_merkens_ssp",SSP,".csv"))
  dfM_diff <- dfM_1 - dfM_0
  max_diff <- rep(0, dim(dfM_diff)[2]-1)
  for (i in 2:(length(max_diff)+1)) {max_diff[i] = max(dfM_1[1:5,i] - dfM_0[1:5,i])}
  dfJ_0 <- read_csv(paste0(ciampath,"data/ssp/popdens_seg_jones_ssp",SSP,"_origCL.csv"))
  dfJ_1 <- read_csv(paste0(ciampath,"data/ssp/popdens_seg_jones_ssp",SSP,".csv"))
  dfJ_diff <- dfJ_1 - dfJ_0
  max_diff <- rep(0, dim(dfJ_diff)[2]-1)
  for (i in 2:(length(max_diff)+1)) {max_diff[i] = max(dfJ_1[,i] - dfJ_0[,i])}
  # continued debugging/checking for Jones et al data set
  idx_diff <- vector("list", nrow(dfJ_diff))
  for (i in 1:nrow(dfJ_diff)) {idx_diff[[i]] <- which(abs(dfJ_diff[i,2:12149]) > 1)}
  # in the first time step it's just the UK Virgin Islands UnitedKingdom8055 segment, off by 3324.814
  # in 2060, jumps up to 229 segments that have a difference greater than 1
  
  #Do stuff like this to check:  t <- 1:10; i <- 3; cbind(dfJ_0[t,idx_diff[[max(t)]][i]+1], dfJ_1[t,idx_diff[[max(t)]][i]+1])
  pseg_samp <- c("Angola1216",          "Australia12036"  ,    "Bahamas7719"    ,     "Bahrain9365"     ,    "Bangladesh1872"   ,   "Canada5220"   ,      
    "Chile4000"    ,       "China10494"     ,     "China2372"       ,    "Egypt1448"       ,    "France9427"    ,      "France9428"   ,      
    "FrenchPolynesia5636", "FrenchPolynesia5765", "FrenchPolynesia5775", "FrenchPolynesia5783", "FrenchPolynesia5799", "UnitedKingdom8055")
  
  # seems odd that in CL's original version, things would just suddenly change in 2050-2060 time step
  
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