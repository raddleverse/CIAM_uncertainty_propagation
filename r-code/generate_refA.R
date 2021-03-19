#### Script to generate reference adaptation level (RefA) for CIAM
# RefA values are produced from GAMS version by running Rcp0p50 NoAdaptation Case for period 0
#   and selecting the optimal adaptation value
# This script prepares GAMS files for use in Julia version of CIAM

## Catherine Ledna, 1/29/2020

#--------------------------------------------------------------
# 1. Read in files
library(tidyverse)
fdir <- "/Users/catherineledna/Documents/gamsdir/projdir/output/noSLRreference/"
files <- list.files(fdir)
files <- files[grepl("refA",files) & !(grepl("htm",files))]
xsc <- read_csv("/Users/catherineledna/.julia/dev/MimiCIAM/data/input/xsc.csv")
segs <- unique(xsc$seg)

temp <- data.frame()
for (f in files){
  fn <- read_csv(paste0(fdir,f),col_names = F) 
  temp <- rbind(temp,fn)
}

temp <- temp %>% rename(segments=X1,variable=X2,value=X3) %>%
  mutate(segments=gsub("\'","",segments))
omitted <- segs[!(segs %in% unique(temp$segments))]
temp2 <- data.frame(segments=omitted,variable="R",value=0)

out <- rbind(temp,temp2)
write_csv(out,"/Users/catherineledna/.julia/dev/MimiCIAM/data/input/refA.csv")
