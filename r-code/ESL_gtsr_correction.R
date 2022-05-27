library(sf)
library(readxl)

# all of the files have the segments in the same order. max differences in 
# lat/lon is 1E-5

rdir <- "/Users/aewsma/codes/CIAM_uncertainty_propagation/R-code"
ciamdir <- "/Users/aewsma/.julia/dev/MimiCIAM"
setwd(rdir)

# read data from GTSR (Muis et al. 2016)
buffer <- read_sf("../ciam-code/data/12712469/data/gtsr_rp100.shp")
gtsr100 <- buffer %>% st_drop_geometry()
dfG <- data.frame(gtsr100)

# read data from DINAS-COAST, from Original CIAM
dfD <- read.csv(paste0(ciamdir,"/data/input/surgeexposure.csv"))
segnames <- dfD$segments
segmap <- read.csv(paste0(ciamdir,"/data/diva_segment_latlon.csv"))

# read data from original DINAS-COAST surge levels
dfS <- read_excel("../ciam-code/data/DIVA_spreadsheet.xls")
# want: SEGID	S1	S10	S100	S1000	BRF	SMAX	LONGI	LATI
dfS <- data.frame(dfS[c("SEGID","S1","S10","S100","S1000","SMAX","LONGI","LATI")])
# ^--- this one is also in the same segment order as the GTSR and DINAS-COAST data

# > lat_g <- dfgtsr$LATI; lat_d <- segmap$lati; lon_g <- dfgtsr$LONGI; lon_d <- segmap$longi
# > plot(abs(lat_g-lat_d))
# > plot(abs(lon_g-lon_d))
# And so is the DINAS-COAST data set:
# > all(segmap$segments == dfdc$segments)

# initialize & shift dfD surge levels by dfS$s1
dfDG <- dfD
dfDG[,2:6] <- dfDG[,2:6] + dfS$S1

dfDG

# multiplicative factor to bias-correct D-C to GTSR. new = factor*old
# will fill up the `s100` column straight from GTSR
factor <- dfG$rp00100/dfS$S100
dfDG$s100 <- dfDG$s100*factor
dfDG$s1 <- dfDG$s1*factor
dfDG$s10 <- dfDG$s10*factor
dfDG$s1000 <- dfDG$s1000*factor
#dfDG$smax <- dfDG$smax*factor # don't shift `smax` - it's an upper bound

# normalize relative to s1 again
if (!all(dfDG$s1[1:5]==0)) {
  dfDG$s10 <- dfDG$s10 - dfDG$s1
  dfDG$s100 <- dfDG$s100 - dfDG$s1
  dfDG$s1000 <- dfDG$s1000 - dfDG$s1
  dfDG$s1 <- rep(0, length(dfDG$s1))
}
# make sure we have the original smax
dfDG$smax = dfD$smax

# write an output file in the style of the original surgeexposure.csv file
# give it a helpful name
write.csv(dfDG, file=paste0(ciamdir,"/data/input/surgeexposure_dc-gtsr.csv"), row.names=FALSE)


###