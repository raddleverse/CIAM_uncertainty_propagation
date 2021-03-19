## Process Model Outputs
# Inputs: CIAM Model Parameters and Results (timeseries and NPV) 
# Outputs: 
#   Mean global annual cost in 2050 and 2100 (and 1, 5, 50, 95 and 99th percentiles)
#   Mean global NPV (plus descriptors)
#   Cost as percent of GDP 
#   Characteristics of distributions
# To Add: Worst and best case parameters across slr min, max 
#   Worst x Min SLR, Worst x Max SLR, Best x Min SLR, Best x Max SLR 

library(tidyverse)
library(reshape2)

outputdir1<-"/Volumes/MASTERS/ciammcs/CIAM2020-10-14 12-03-32MC2880Reg1/results/"
outputdir2 <- "/Volumes/MASTERS/ciammcs/CIAM2020-10-17 20-35-15MC2880Reg2/results/"
outputdir3 <- "/Volumes/MASTERS/ciammcs/CIAM2020-10-19 13-19-57MC2880Reg3/results/"
outputdir4 <- "/Volumes/MASTERS/ciammcs/CIAM2020-10-27 11-38-54MC2880Reg4/results/"

# Socioeconomic Parameters
gdp_ssp5_file <- "/Users/catherineledna/.julia/dev/MimiCIAM/data/ssp/ypcc_IIASAGDP_SSP5_v9_130219.csv"
pop_ssp5_file <- "/Users/catherineledna/.julia/dev/MimiCIAM/data/ssp/pop_IIASAGDP_SSP5_v9_130219.csv"

# Read Data 
ts1<-read_csv(paste0(outputdir1,"globalts_85.csv")) %>% mutate(fixed="Flexible",noRetreat="With Retreat", regime=1)
ts2 <- read_csv(paste0(outputdir2,"globalts_85.csv")) %>% mutate(fixed="Flexible",noRetreat="No Retreat", regime=2)
ts3 <- read_csv(paste0(outputdir3,"globalts_85.csv")) %>% mutate(fixed="Fixed",noRetreat="With Retreat", regime=3)
ts4 <- read_csv(paste0(outputdir4,"globalts_85.csv")) %>% mutate(fixed="Fixed",noRetreat="No Retreat", regime=4)

npv1 <- read_csv(paste0(outputdir1,"globalnpv.csv")) %>% mutate(fixed="Flexible",noRetreat="With Retreat",regime=1)
npv2 <- read_csv(paste0(outputdir2,"globalnpv.csv")) %>% mutate(fixed="Flexible",noRetreat="No Retreat",regime=2)
npv3 <- read_csv(paste0(outputdir3,"globalnpv.csv")) %>% mutate(fixed="Fixed",noRetreat="With Retreat",regime=3)
npv4 <- read_csv(paste0(outputdir4,"globalnpv.csv")) %>% mutate(fixed="Fixed",noRetreat="No Retreat",regime=4)

params1<-read_csv(paste0(outputdir1,"trials.csv"))
params2<-read_csv(paste0(outputdir2,"trials.csv"))
params3<-read_csv(paste0(outputdir3,"trials.csv"))
params4 <- read_csv(paste0(outputdir4,"trials.csv"))

npv1 <- cbind(npv1,params1)
npv2 <- cbind(npv2,params2)
npv3 <- cbind(npv3,params3)
npv4 <- cbind(npv4,params3)

gdp_ssp5 <- read_csv(gdp_ssp5_file) %>% rename(time=variable) %>% melt(id.vars=c("time")) %>% filter(time <=10) %>% mutate(variable=as.character(variable))
pop_ssp5 <- read_csv(pop_ssp5_file) %>% rename(time=variable) %>% melt(id.vars=c("time")) %>% filter(time <=10) %>% mutate(variable=as.character(variable)) %>% rename(pop=value)
ssp5_data <- pop_ssp5 %>% mutate(gdp_percap=gdp_ssp5$value[match(paste0(time,variable),paste0(gdp_ssp5$time,gdp_ssp5$variable))],
                                 gdp_M = gdp_percap * pop)

# Unit $Million 2015 USD
usa_2020_gdp <- ssp5_data$gdp_M[ssp5_data$variable=="USA"&ssp5_data$time==2]
global_2020_gdp <- sum(ssp5_data$gdp_M[ssp5_data$time==2]) 

# DIAZ - CIAM 2016 - Total US Cost 2050 - 0.08% of USA GDP 
# Global - 270billion-2.2 trillion annually in 2100 (0.2-1.5% global projected GDP (147.6 trillion))
# Also reports: median regional cost as pct of gdp (5th and 95th percentile and country), countries with largest burden of GDP
#   Countries with largest NPV 

# Plot 1: Time Series By Main Category
#   Step 1: ID quantile (by NPV)
process_ts_1a <- function(ts){
  npvdf <- ts %>% filter(regions=="global",costtype=="OptimalCost") %>% group_by(ens,time,regions,segments,segID,noRetreat,regime)%>%
    summarize(cost=sum(cost),pct_gdp=sum(pct_gdp)) %>% ungroup() %>% mutate(df=1/(1 + .04)^(10 * (time-1)),npv=cost*df) %>% group_by(ens,regions,segments,segID,noRetreat,regime)%>%
    mutate(npv=sum(npv))
  
  npvdf_filt <- npvdf %>% filter(time==10) %>% select(ens,time,regions,segments,segID,noRetreat,regime,npv)
  npvq <- quantile(npvdf_filt$npv,c(.05,.5,.95))
  inds <- sapply(npvq,FUN=function(x) which.min(abs(npvdf_filt$npv-x)))
  
  bounds_filt <- npvdf %>% filter(ens %in% inds) %>% mutate(pct=ifelse(ens==inds[1],.05,ifelse(ens==inds[2],.5,.95)))
  ts_filt <- ts %>% filter(ens %in% inds,costtype=="OptimalCost",regions=="global") %>% group_by(ens,time,regions,segments,segID,noRetreat,regime,category) %>%
    summarize(cost=sum(cost),pct_gdp=sum(pct_gdp)) %>% ungroup() %>% mutate(p05=bounds_filt$cost[match(paste0(time,regions,segments,segID,noRetreat,regime,0.05),
                                                                                                       paste0(bounds_filt$time,bounds_filt$regions,bounds_filt$segments,bounds_filt$segID,bounds_filt$noRetreat,bounds_filt$regime,bounds_filt$pct))],
                                                                            p95=bounds_filt$cost[match(paste0(time,regions,segments,segID,noRetreat,regime,0.95),
                                                                                                       paste0(bounds_filt$time,bounds_filt$regions,bounds_filt$segments,bounds_filt$segID,bounds_filt$noRetreat,bounds_filt$regime,bounds_filt$pct))]) %>%
    filter(ens==inds[2]) %>% mutate(year=2000+time*10)
  
  return(ts_filt)
  
}

process_ts_1b <- function(ts){
  npvdf <- ts %>% filter(regions=="global",costtype=="OptimalCost") %>% group_by(ens,time,regions,segments,segID,noRetreat,regime)%>%
    summarize(cost=sum(cost),pct_gdp=sum(pct_gdp)) %>% ungroup() %>% mutate(df=1/(1 + .04)^(10 * (time-1)),npv=cost*df) %>% group_by(ens,regions,segments,segID,noRetreat,regime)%>%
    mutate(npv=sum(npv))
  
  npvdf_filt <- npvdf %>% filter(time==10) %>% select(ens,time,regions,segments,segID,noRetreat,regime,npv)
  npvq <- quantile(npvdf_filt$npv,c(.05,.5,.95))
  inds <- sapply(npvq,FUN=function(x) which.min(abs(npvdf_filt$npv-x)))
  
  bounds_filt <- npvdf %>% filter(ens %in% inds) %>% mutate(pct=ifelse(ens==inds[1],.05,ifelse(ens==inds[2],.5,.95)))
  ts_filt <- ts %>% filter(ens %in% inds,regions=="global",costtype!="OptimalCost") %>% group_by(ens,time,regions,segments,segID,noRetreat,regime,costtype) %>%
    summarize(cost=sum(cost),pct_gdp=sum(pct_gdp)) %>% ungroup() %>% mutate(p05=bounds_filt$cost[match(paste0(time,regions,segments,segID,noRetreat,regime,0.05),
                                                                                                       paste0(bounds_filt$time,bounds_filt$regions,bounds_filt$segments,bounds_filt$segID,bounds_filt$noRetreat,bounds_filt$regime,bounds_filt$pct))],
                                                                            p95=bounds_filt$cost[match(paste0(time,regions,segments,segID,noRetreat,regime,0.95),
                                                                                                       paste0(bounds_filt$time,bounds_filt$regions,bounds_filt$segments,bounds_filt$segID,bounds_filt$noRetreat,bounds_filt$regime,bounds_filt$pct))]) %>%
    filter(ens==inds[2]) %>% mutate(year=2000+time*10)
  
  return(ts_filt)
}

process_totals <- function(ts){
  npvdf <- ts %>% filter(costtype=="OptimalCost") %>% group_by(ens,time,regions,segments,segID,noRetreat,regime)%>%
    summarize(cost=sum(cost),pct_gdp=sum(pct_gdp)) %>% ungroup() %>% mutate(df=1/(1 + .04)^(10 * (time-1)),npv=cost*df) %>% group_by(ens,regions,segments,segID,noRetreat,regime)%>%
    mutate(npv=sum(npv)) %>% ungroup()
  
  npvdf_filt <- npvdf %>% filter(time==10) %>% select(ens,time,regions,segments,segID,noRetreat,regime,npv,cost)
  npvdf_filt_usa <- npvdf_filt %>% filter(regions=="USA")
  npvdf_filt_glob <- npvdf_filt %>% filter(regions=="global")
  npvqu <- quantile(npvdf_filt_usa$npv,c(.05,.5,.95))
  indsu <- sapply(npvqu,FUN=function(x) which.min(abs(npvdf_filt_usa$npv-x)))
  npvqg <- quantile(npvdf_filt_glob$npv,c(.05,.5,.95))
  indsg <- sapply(npvqg,FUN=function(x) which.min(abs(npvdf_filt_glob$npv-x)))
  
  bounds_filt_usa <- npvdf %>% filter(ens %in% indsu,regions=="USA") %>% mutate(pct=ifelse(ens==indsu[1],.05,ifelse(ens==indsu[2],.5,.95)))
  bounds_filt_glob <- npvdf%>% filter(ens %in% indsg, regions=="global") %>% mutate(pct=ifelse(ens==indsg[1],.05,ifelse(ens==indsg[2],.5,.95)))

  
  ts_filt_usa <- ts %>% filter(ens %in% c(indsu["50%"]),costtype=="OptimalCost",regions=="USA") %>% group_by(ens,time,regions,segments,segID,noRetreat,regime) %>%
    summarize(cost=sum(cost),pct_gdp=sum(pct_gdp)) %>% ungroup() %>% mutate(p05=bounds_filt_usa$cost[match(paste0(time,regions,segments,segID,noRetreat,regime,0.05),
                                                                                                       paste0(bounds_filt_usa$time,bounds_filt_usa$regions,bounds_filt_usa$segments,bounds_filt_usa$segID,bounds_filt_usa$noRetreat,bounds_filt_usa$regime,bounds_filt_usa$pct))],
                                                                            pgdp05=bounds_filt_usa$pct_gdp[match(paste0(time,regions,segments,segID,noRetreat,regime,0.05),
                                                                                                           paste0(bounds_filt_usa$time,bounds_filt_usa$regions,bounds_filt_usa$segments,bounds_filt_usa$segID,bounds_filt_usa$noRetreat,bounds_filt_usa$regime,bounds_filt_usa$pct))],
                                                                            p95=bounds_filt_usa$cost[match(paste0(time,regions,segments,segID,noRetreat,regime,0.95),
                                                                                                       paste0(bounds_filt_usa$time,bounds_filt_usa$regions,bounds_filt_usa$segments,bounds_filt_usa$segID,bounds_filt_usa$noRetreat,bounds_filt_usa$regime,bounds_filt_usa$pct))],
                                                                            pgdp95=bounds_filt_usa$pct_gdp[match(paste0(time,regions,segments,segID,noRetreat,regime,0.95),
                                                                                                           paste0(bounds_filt_usa$time,bounds_filt_usa$regions,bounds_filt_usa$segments,bounds_filt_usa$segID,bounds_filt_usa$noRetreat,bounds_filt_usa$regime,bounds_filt_usa$pct))]) %>%
    mutate(year=2000+time*10)
  
  ts_filt_glob <- ts %>% filter(ens ==indsg["50%"],costtype=="OptimalCost",regions=="global") %>% group_by(ens,time,regions,segments,segID,noRetreat,regime) %>%
    summarize(cost=sum(cost),pct_gdp=sum(pct_gdp)) %>% ungroup() %>% mutate(p05=bounds_filt_glob$cost[match(paste0(time,regions,segments,segID,noRetreat,regime,0.05),
                                                                                                       paste0(bounds_filt_glob$time,bounds_filt_glob$regions,bounds_filt_glob$segments,bounds_filt_glob$segID,bounds_filt_glob$noRetreat,bounds_filt_glob$regime,bounds_filt_glob$pct))],
                                                                            pgdp05 = bounds_filt_glob$pct_gdp[match(paste0(time,regions,segments,segID,noRetreat,regime,0.05),
                                                                                                                 paste0(bounds_filt_glob$time,bounds_filt_glob$regions,bounds_filt_glob$segments,bounds_filt_glob$segID,bounds_filt_glob$noRetreat,bounds_filt_glob$regime,bounds_filt_glob$pct))],
                                                                            p95=bounds_filt_glob$cost[match(paste0(time,regions,segments,segID,noRetreat,regime,0.95),
                                                                                                       paste0(bounds_filt_glob$time,bounds_filt_glob$regions,bounds_filt_glob$segments,bounds_filt_glob$segID,bounds_filt_glob$noRetreat,bounds_filt_glob$regime,bounds_filt_glob$pct))],
                                                                            pgdp95=bounds_filt_glob$pct_gdp[match(paste0(time,regions,segments,segID,noRetreat,regime,0.95),
                                                                                                               paste0(bounds_filt_glob$time,bounds_filt_glob$regions,bounds_filt_glob$segments,bounds_filt_glob$segID,bounds_filt_glob$noRetreat,bounds_filt_glob$regime,bounds_filt_glob$pct))]) %>%
    mutate(year=2000+time*10)
  ts_filt <- rbind(ts_filt_usa,ts_filt_glob) %>% mutate(pgdp_2020_p50=ifelse(regions=="USA",cost/(usa_2020_gdp/1e3),cost/(global_2020_gdp/1e3)),
                                                        pgdp_2020_p5=ifelse(regions=="USA",p05/(usa_2020_gdp/1e3),p05/(global_2020_gdp/1e3)),
                                                        pgdp_2020_p95=ifelse(regions=="USA",p95/(usa_2020_gdp/1e3),p95/(global_2020_gdp/1e3)))
  
  return(ts_filt)
  
  
}

ts1_filt <- process_ts(ts1)
ts2_filt <- process_ts(ts2)
ts3_filt <- process_ts(ts3)
ts_filt <- rbind(ts1_filt,ts2_filt,ts3_filt) %>% mutate(regime=paste0("Regime ",regime))

ts1b_filt <- process_ts_1b(ts1)
ts2b_filt <- process_ts_1b(ts2)
ts3b_filt <- process_ts_1b(ts3)
ts_filt_b <- rbind(ts1b_filt,ts2b_filt,ts3b_filt)%>% mutate(regime=paste0("Regime ",regime))

# Output Results - Global and US 5, 50, 95th percentile costs in 2050 and 2100 
totals1 <- process_totals(ts1)
totals2 <- process_totals(ts2)
totals3 <- process_totals(ts3)
totals <- rbind(totals1,totals2,totals3) %>% mutate(df =1/(1 + .04)^(10 * (time-1)), npv_int = cost*10 * df,
                                                    npv_int_05 = p05*10*df, npv_int_95 =p95*10*df) %>% 
  group_by(ens,regions,segments,segID,noRetreat,regime) %>%
  mutate(npv = sum(npv_int), npv05=sum(npv_int_05), npv95=sum(npv_int_95)) %>% ungroup()
        
npv_df <- rbind(npv1,npv2,npv3) %>% mutate(npv=npv*10) 



npv_mean_fd <- rbind(npv1,npv2,npv3) %>% mutate(npv=npv*10) %>% group_by(brick,retreat,noRetreat,regime) %>%
  summarize(mean_npv=mean(npv),sd=sd(npv)) %>% ungroup() %>% group_by(regime,brick,retreat)%>%
  mutate(perc_greater_1sd = ifelse(regime==1, sum((npv1$npv*10)>(mean_npv+sd))/2880,
                                   ifelse(regime==2,sum((npv2$npv*10)>(mean_npv+sd))/2880,
                                          sum((npv3$npv*10)>(mean_npv+sd))/2880)))
npv_mean_regime <-  rbind(npv1,npv2,npv3) %>% mutate(npv=npv*10) %>% group_by(regime) %>%
  summarize(mean_npv=mean(npv),sd=sd(npv)) %>% ungroup() %>%
  group_by(regime) %>% mutate(perc_greater_1sd = ifelse(regime==1, sum((npv1$npv*10)>(mean_npv+sd))/2880,
                                                        ifelse(regime==2,sum((npv2$npv*10)>(mean_npv+sd))/2880,
                                                               sum((npv3$npv*10)>(mean_npv+sd))/2880)))
npv_mean_ret <- rbind(npv1,npv2,npv3) %>% mutate(npv=npv*10) %>% group_by(retreat) %>%
  summarize(mean_npv=mean(npv),sd=sd(npv)) %>% ungroup()


# Plot 1a: Comparing Adaptation Regimes 
p1 <- ggplot(ts_filt)+
  facet_grid(cols=vars(regime))+
  geom_bar(aes(year,cost,fill=factor(category)),stat="identity")+
  ylab("Annual Cost (Billion 2010USD/year)")+
  xlab("")+
  theme_bw()+
  labs(fill="Cost Type")+
  geom_errorbar(data=ts_filt,aes(x=year,ymin=p05,ymax=p95))

# Plot 1b - to do: color scheme 
p1b <- ggplot(ts_filt_b)+
  facet_grid(cols=vars(regime))+
  geom_bar(aes(year,cost,fill=factor(costtype)),stat="identity")+
  ylab("Annual Cost (Billion 2010USD/year)")+
  xlab("")+
  theme_bw()+
  labs(fill="Cost Type")+
  geom_errorbar(data=ts_filt,aes(x=year,ymipn=p05,ymax=p95))


# Plot 2: Evaluating the distribution of outcomes 


npv1 <- npv1 %>% mutate(ens=match(paste0(round(npv,digits = 4),noRetreat,regime),paste0(round(npvdf1_filt$npv,digits = 4),npvdf1_filt$noRetreat,npvdf1_filt$regime)))



ts1_filt <- ts1 %>% filter(regions=="global",ens %in% inds1,costtype=="OptimalCost") %>% group_by(ens,time,regions,segments,segID,category,noRetreat,regime)%>%
  summarize(cost=sum(cost),pct_gdp=sum(pct_gdp)) %>% ungroup() %>% group_by(ens,time,regions,segments,segID,noRetreat,regime)%>%
  mutate(totcost=sum(cost)) %>% ungroup() %>% mutate(df=1/(1 + .04)^(10 * (time-1)),npv=totcost*df)

# Plot 2: Time Series by Subcategory 


# Plot 3: Variance by Regime vs Variance by Geophysical Attribute 
p1 <- ggplot(npv_df)+
  facet_grid(cols=vars(regime))+
  geom_density(aes(npv))


# Diagnostic Plot: Distribution of NPV vs GMSL
p1 <- ggplot(npv1)+geom_point(aes(gmsl,npv,color=factor(brick)))
p2 <- ggplot(npv1)+geom_density(aes(npv))
res<- cor(npv1_corr) # Correlation Matrix - shows strongest positive correlation is with gmsl 

# Calculate 5th, 50th and 95th percentiles 
#   Question: As fraction of NPV or GMSL? - going with NPV for now 
npvq <- quantile(npv1$npv,c(.05,.5,.95))
gmslq <- quantile(npv1$gmsl,c(0.05,.5,.95))
inds<-sapply(npvq,FUN=function(x) which.min(abs(npv1$npv-x)))

inds2 <- sapply(gmslq,FUN=function(x) which.min(abs(npv1$gmsl-x))) # 2480,
npv1_q <- npv1[inds,] %>% mutate(percentile=c(5,50,95))
npv1_q2 <- npv1[inds2,] %>% mutate(gmslpercentile=c(5,50,95))

# Figure out timeseries 
ts_q <- ts1 %>% filter(regions=="global",costtype=="OptimalCost") %>% group_by(ens,time,regions,segments,segID,costtype)%>%
  summarize(cost=sum(cost),pct_gdp=sum(pct_gdp)) %>% ungroup() %>% mutate(df=1/(1 + .04)^(10 * (time-1))) %>% 
  group_by(ens,regions,segments,segID) %>% summarize(npv=sum(cost*df))


tsq_glob <- ts_q %>% filter(regions=="global") %>% group_by(ens,time,regions,segments,segID,category,costtype)%>%
  summarize(cost=sum(cost),pct_gdp=sum(pct_gdp)) %>% ungroup() %>% mutate(percentile=ifelse(ens==1517,50,ifelse(ens==2472,5,95)))

tsq_opcost <- tsq_glob %>% filter(costtype=="OptimalCost") %>% group_by(ens,time,regions,segments,segID) %>% mutate(globcost=sum(cost),globpctgdp=sum(pct_gdp))%>%
  ungroup() %>% mutate(disfac=1/(1 + .04)^(10 * (time-1))) %>% group_by(ens) %>% mutate(npv=sum(unique(globcost)*unique(disfac)))
