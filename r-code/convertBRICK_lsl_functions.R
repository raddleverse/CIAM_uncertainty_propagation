# Need GIS, GSIC, AIS, TE, LWS
# Dimensions are 251,10589 for 251 timesteps (1850 to 2100) and 10589 runs 
# Need to arrange by percentile 
library(tidyverse)
library(ncdf4)
library(RNetCDF)
source('BRICK_LSL.R') # function to convert brick gsl to lsl for given lon/lat point 

# Function to load BRICK NetCDF components and return them as a list of matrices 
loadBRICK_components_netcdf<- function(rcp,components){
  header <- paste0(DATADIR,"BRICK-")
  
  for (c in components){
    f<- nc_open(paste0(header,c,"_",rcp,".nc"))
    assign(paste0(c,".var"),ncvar_get(f, paste0(c,"_",rcp)))
    nc_close(f)
  }
  
  outnames <- sapply(components,paste0,".var")
  out <- lapply(outnames, get,envir = environment())
  return(out)
  
}

# Load CSV version of BRICK outputs
loadBRICK_components_csv <- function(rcp,components,ci=T){
  if (ci==T){
    header <- paste0(DATADIR,"BRICK-ci-")
  }else{
    header <- paste0(DATADIR,"BRICK-")
  }
  
  for (c in components){
    f<- read_csv(paste0(header,c,"_",rcp,".csv")) %>% rename(p50=Mean_Chain,p10=LowerConf_0.9,p90=UpperConf_0.9,
                                                             p5=LowerConf_0.95,p95=UpperConf_0.95)
    assign(paste0(c,".var"),f,envir = globalenv())
  }
  
}

# Function to extract given percentile 
# components argument must contain GlobalSeaLevel 
getPercentile <- function(rcp,components,pct){
  brickComps <- loadBRICK_components(rcp,components)
  
  # Get dimensions (assumes all are of same dim)
  t_dim <- dim(brickComps$GlobalSeaLevel)[1] # time dimension 
  r_dim <- dim(brickComps$GlobalSeaLevel)[2] # number of runs 
  
  # For base component, use Global Sea Level
  # Find the index of fifth, 95th and median percentile
  sort_order <- order(brickComps$GlobalSeaLevel[t_dim,]) # Order of the indices for sorted matrix
  
  position <- round(pct*r_dim) # Fifth percentile index for sort order; rounds down to integer  
  if (pct==0) position <- 1 else if (pct==1) position <- r_dim #edge cases 
  
  ind <- sort_order[position] # Fifth percentile index
  
  # Extract the percentiles from the main list 
  newcomps <- components[components!="GlobalSeaLevel"]
  brickCompsPct <- vector("list",length(newcomps))
  names(brickCompsPct) <- newcomps
  
  for (n in newcomps){
    brickCompsPct[[n]] <- brickComps[[n]][,ind]
  }
  
  return(brickCompsPct)
  
}

getPercentileGMSL <- function(rcp,components,pct){
  brickComps <- loadBRICK_components(rcp,components)
  
  # Get dimensions (assumes all are of same dim)
  t_dim <- dim(brickComps$GlobalSeaLevel)[1] # time dimension 
  r_dim <- dim(brickComps$GlobalSeaLevel)[2] # number of runs 
  
  # For base component, use Global Sea Level
  # Find the index of fifth, 95th and median percentile
  sort_order <- order(brickComps$GlobalSeaLevel[t_dim,]) # Order of the indices for sorted matrix
  
  position <- round(pct*r_dim) # Fifth percentile index for sort order; rounds down to integer  
  if (pct==0) position <- 1 else if (pct==1) position <- r_dim #edge cases 
  
  ind <- sort_order[position] # Fifth percentile index
  
  # Extract the percentiles from the main list 
  newcomps <- components
  brickCompsPct <- vector("list",length(newcomps))
  names(brickCompsPct) <- newcomps
  
  for (n in newcomps){
    brickCompsPct[[n]] <- brickComps[[n]][,ind]
  }
  
  return(brickCompsPct)
  
}
