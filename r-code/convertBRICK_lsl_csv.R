### Code to convert BRICK gsl projections to lsl for desired lon/lat points 
# Catherine Ledna
# 3-29-19

####-----------------------------------------------------------------------
# Step 1: Source necessary functions and load libraries 

####-----------------------------------------------------------------------
source("convertBRICK_lsl_functions.R")
DATADIR<- "../../data/brick/ci_results_BRICK_100kruns/" 


####-----------------------------------------------------------------------
#Step 2: Run the code for desired percentiles 
# Code is for CSV run on 1/29/2020
####-----------------------------------------------------------------------

### Define Control Variables
RCP<-c("RCP85")
comps <- c("GlobalSeaLevel","AIS","TE","GIS","GSIC") # order of brick components to load 
pctls <- c("p50")#c("p5","p10","p90","p95","p50") # percentiles to extract 
#time_index <- 131:251 # 1980-2100 from 1850-2100 
YSTART=2010
YEND=2200
yearnames <- YSTART:YEND
yearnames_dec <- seq(YSTART,YEND,10)


### Load CIAM lonlat data
ciam_lonlat <- read_csv("../../data/ciam/diva_segment_latlon.csv") %>% select(-segid)

### Extract BRICK global sea level percentile matrices 
loadBRICK_components_csv(RCP[1],comps)
BRICKpcts <- vector("list",length(pctls)*length(RCP))
names(BRICKpcts)<- apply(expand.grid(RCP,pctls),1,function(x) paste(x,collapse=""))

BRICKpctsgmsl <- vector("list",length(pctls)*length(RCP))
names(BRICKpctsgmsl)<- apply(expand.grid(RCP,names(pctls)),1,function(x) paste(x,collapse=""))

for (i in 1:length(pctls)){
  for (j in 1:length(RCP)){
    BRICKcomps <- vector("list",length(comps))
    names(BRICKcomps)<- comps
    
    for (k in 1:length(comps)){
      val <- get(paste0(comps[k],".var")) %>% select_("Year",pctls[i]) %>% filter(Year>=YSTART & Year<=YEND)
      compval <- val[[pctls[i]]]
      BRICKcomps[[comps[k]]] = compval
    }
    
    #val2 <- getPercentileGMSL(RCP[j],comps,pctls[i])
    #n <- names(val)
    #n2 <- names(val2)
    # Subset val to desired indices: 2010-2100 
    #val <- lapply(seq_along(val), function(x) val[[x]][time_index])
    #val2 <- lapply(seq_along(val2), function(x) val2[[x]][time_index])
    #names(val)<- n
    #names(val2)<- n2
    
    name <- paste0(RCP[j],pctls[i])
    BRICKpcts[[name]]<- BRICKcomps
    #BRICKpctsgmsl[[name]]<- val2
  }
}

### Save GMSL values
# BRICKgmsl <- lapply(BRICKpctsgmsl, "[[", "GlobalSeaLevel")
# BRICKgmsl <- as.data.frame(BRICKgmsl)
# BRICKgmsl$time <- yearnames
# write_csv(BRICKgmsl, "../processed-data/lslr/BRICKgmsl_wfd.csv")


### Convert Global Sea Level to Local Sea Level 
# Result will be 11 x 12148 series saved as CSV (all segments, 2010-2100; 10 year increments)
for (i in 1:length(BRICKpcts)){
  
  brickdata <- BRICKpcts[[i]]
  brickname <- names(BRICKpcts[i])
  
  # Create matrix to store output 
  outdata <- matrix(nrow=length(YSTART:YEND), ncol=nrow(ciam_lonlat))
  segs <- ciam_lonlat$segments
  rownames(outdata)<- yearnames
  colnames(outdata)<- segs
  na_segs <- c()

  for (j in 1:nrow(ciam_lonlat)){
    longi <- ciam_lonlat$longi[j]
    lati <- ciam_lonlat$lati[j]
    segname <- ciam_lonlat$segments[j]
    
    lsl_pt <- brick_lsl(lat.in = lati, lon.in = longi,global = F,n.time = length(yearnames),slr_gis = brickdata$GIS,
                         slr_gsic = brickdata$GSIC, slr_ais = brickdata$AIS, slr_te = brickdata$TE)
    lsl_pt <- as.vector(lsl_pt)

    if (is.na(lsl_pt[1])) na_segs <- c(na_segs,segname)
    outdata[,segname] <- lsl_pt
  }
  
  # Save output as a CSV 
  out_df <- as.data.frame(outdata) %>% mutate(variable=yearnames) %>% select(variable,everything()) %>% filter(variable %in% yearnames_dec)
  out_df_ann <- as.data.frame(outdata) %>% mutate(variable=yearnames) %>% select(variable,everything())
  
  write.csv(out_df,paste0("../../data/brick-lsl/BRICKsneasy-lsl_",brickname,".csv"),row.names=F)
  write.csv(out_df_ann, paste0("../../data/brick-lsl//BRICKsneasy-lsl_",brickname,"_annual.csv"),row.names=F)
}

### Median Local Sea Level Projection (for table purposes)
#   source initAnalysis.R first. 
med_df <- data.frame()
for (i in 1:length(names(BRICKpcts))){
  df <- read_csv(paste0("../processed-data/lslr/BRICKfd-lsl_",names(BRICKpcts)[i],".csv")) %>% filter(variable==2100)%>% rename(time=variable) %>% 
    melt(id.vars=c("time"))%>%mutate(variable=as.character(variable))%>% filter(variable %in% seg_names)
  med <- median(df$value)
  f <- quantile(df$value,.05)
  n <- quantile(df$value,.95)
  
  temp <- data.frame(name=names(BRICKpcts)[i],med=med,f=f,n=n)
  med_df <- rbind(med_df, temp)
}
