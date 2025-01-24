################################################################################
#
# Examples reproducing the results published in the article:
# The Association between Total Precipitation and Diarrhea Morbidity: A Multi-Country Study across Diverse # Climate Zones'
# Updated on 2025.01.24. by Rui PAN 
#
# Note:
# 1. The code only include the main steps for the main model and results (e.g., Figures 2, 3, 4);
# 2. The example data were simulate data. I show the main inforamtion of data in the R script for better understanding the code.
# 3. Please feel free contact the first author (ruipan@m.u-tokyo.ac.jp) if you have any question or comments about the code. Thank you very much!
#
################################################################################


################################################################################
# PREPARE THE PACKAGE

# Sys.setenv(LANG="en"); Sys.setlocale("LC_ALL","C") # day of week in English
gc()
rm(list=ls())

# Load Some Packages
library(dlnm) ; library(mixmeta) ; library(splines)
library(tsModel) ; library(mgcv)
library(foreach) ; library(doParallel)
library(MASS) ; library(abind)
library(dplyr) ; library(data.table)
library(magrittr) ; library(scales)
library(Epi)
library(ggplot2)
library(RcppRoll)
library(mvmeta)
library(gnm)
library(foreign)
library(tidyverse)
library(rworldmap)

################################################################################
# STEP1: PREPARE THE DATA FOR ANALYSIS

# Prepare the data for analysis -------------------------------------------

# Read data
# setwd()
setwd("C:/Users/45529/Desktop/课题相关/20230227 diarrhoea/20230512 tidy data_other countries/999all/06_new_data_std_01/04_prec_analysis_code/01_code_04_correct_rain_20241122/01_code_submission_01")
country_all <- read.csv("sim_diar_01.csv") # Weekly time series diarrhea data
meta_ind <- read.csv("sim_meta_01.csv") # indicator for meta-regression

head(country_all)
head(meta_ind)

# Important variables in main model for analysis
# city: location code
# cases: weekly diarrhea cases in each locations
# tmean: weekly mean temperature
# prec: natural logarithm of weekly total precipitation
# stratum: per 4 weeks within same year
# prec02: weekly total precipitation

# CREATE THE METADATA AND A LIST OF CITY-SPECIFIC DATASETS
cities <- unique(country_all[c("country","city")])
dlist <- split(country_all, country_all$city)[cities$city]

# ORDER THE LIST OF DATASETS
dlist <- dlist[cities$city]

# RENAME CITIES
rownames(cities) <- NULL


# Prepare the main parameters for the analysis ----------------------------

# SPECIFICATION OF THE EXPOSURE FUNCTION
varfun = "bs"
vardegree = 2
varper <- c(33,66)

cenper <- 50 

# SPECIFICATION OF THE LAG FUNCTION
lag <- 8




################################################################################
# STEP 2: FIRST STAGE MODEL


# NA excluding
options(na.action="na.exclude")

# Save coefficient and co-variance
coef <- matrix(NA,nrow(cities),length(varper)+vardegree, dimnames=list(cities$city))
vcov <- vector("list",nrow(cities))
names(vcov) <- cities$city

# Run the loop ------------------------------------------------------------

# LOOP
time <- proc.time()[3]

for(i in seq(length(dlist))) {
  
  # PRINT
  cat(i,"")
  
  # EXTRACT THE DATA
  data <- dlist[[i]]
  head(data)
  
  # DEFINE THE CROSSBASIS
  argvar <- list(fun=varfun,knots=quantile(data$prec,varper/100,na.rm=T),degree=vardegree)
  cb <- crossbasis(data$prec,lag=8, argvar=argvar,arglag=list(knots=logknots(8,2)))
  # summary(cb)
  
  # Exclude the stratum with no cases
  data %>% group_by(stratum) %>%
    mutate(tot=sum(cases, na.rm=T)) %>%
    mutate(keep=ifelse(tot>0, TRUE, FALSE)) -> data
  
  # Main model
  model <- gnm(cases ~ cb + ns(runMean(tmean,lags=0:4), df=3), data=data, family=quasipoisson,
               eliminate=factor(stratum), subset=keep)
  
  # REDUCTION TO OVERALL CUMULATIVE
  red <- crossreduce(cb,model,cen=quantile(data$prec,cenper/100,na.rm=T))
  coef[i,] <- coef(red)
  vcov[[i]] <- vcov(red)
  
}

proc.time()[3]-time



################################################################################
# STEP 3: SECOND STAGE MODEL


# Multi-level mixed meta_regression --------------------------------------------

# Mixed meta-regression
meta_cc <- mixmeta(coef,vcov,random=~1|country/city, meta_ind,control=list(showiter=T, igls.inititer=10)) # For country BULP
meta_clim <- mixmeta(coef~clim,vcov,random=~1|country/city, meta_ind,control=list(showiter=T, igls.inititer=10)) # For climate prediction

summary(meta_cc)
summary(meta_clim)

# Save meta-regression results: It take long time to run the model for the ture data
# meta_l2_list <- list(meta_cc,meta_clim)
# saveRDS(meta_l2_list,file="meta_level2_list_02.Rds")


# Main results at country level ------------------------------------------------------------

# Read saved meta-regression model
summary(meta_cc)

# BULP at country and location level
blup_country <- unique(blup(meta_cc,vcov=T,level=1)) # Country level
blup_city <- unique(blup(meta_cc,vcov=T,level=2)) # Location level


# Prepare to obtain the total precipitation at country level
country_name_01 <- unique(country_all[c("country","name_country")])

# Split data at country level
country_list <- unique(country_all$country) # Country variabales
dlist3 <- split(country_all, country_all$country)[country_list] # Split by country 

# Save average of total precipitation at country level
prec_country <- list()

# Run the loop
for(i in seq(length(country_list))) {
  
  country <- country_list[i]
  
  prec_all_obs <-unlist(sapply(dlist3[country], function(x) x$prec), use.names=F)
  prec_all_obs <- as.data.frame(prec_all_obs)
  prec_country[i] <- prec_all_obs
  names(prec_country)[i] <- country
  
}



# Predict the curves and save the RRs 
# Create the table to save the RR
tab_rr <- matrix(NA, length(country_list) ,16)
rownames(tab_rr) <- country_list
colnames(tab_rr) <- c("minperc", "minrain", "maxperc1", "maxrain1", "maxperc2", "maxrain2",
                      "RR1","RRlow1","RRhigh1","RR2","RRlow2","RRhigh2","fit1","se1","fit2","se2")

# Save the predict model at country level
cpmod2 <- list()

# save some predict numbers for ggplot to show the association curves
cnm <- c("ctry","loc","knots","lprec","prec","rrmean","rrlci","rruci")
dfknots <- data.frame(matrix(NA,nrow=0,ncol=length(cnm),dimnames=list(NULL,cnm)))
dfknots <- data.frame()

# Percentiles for predict the curves
predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))

# Run the loop
for(i in seq(length(country_list))) {
  
  cat(i,"")
  
  # Create the variable to predict
  prec_country_per <- quantile(prec_country[[i]],predper/100,na.rm=T)
  predvar <- quantile(prec_country[[i]],1:99/100,na.rm=T) # Excluding some extreme values at country level
  
  # REDEFINE THE FUNCTION USING ALL THE ARGUMENTS (BOUNDARY KNOTS INCLUDED)
  argvar <- list(x=predvar,fun=varfun,
                 knots=quantile(prec_country[[i]],varper/100,na.rm=T),
                 Bound=range(prec_country_per[c("0.0%","100.0%")],na.rm=T))
  if(!is.null(vardegree)) argvar$degree <- vardegree
  bvar <- do.call(onebasis,argvar)
  
  # Predict the curves
  cp1 <- crosspred(bvar,coef=blup_country[[i]]$blup,vcov=blup_country[[i]]$vcov,
                   at=prec_country[[i]],cen=quantile(prec_country[[i]],50/100,na.rm=T),
                   model.link="log",by=0.01)
  
  # Define the MRP from 1st to 99th
  minper <- (1:99)[which.min((bvar%*%coef(cp1))[1:99,])]
  mmt <- quantile(prec_country[[i]],minper/100,na.rm=T) # MRP
  
  # Get the 1st and 99th percentile for RR calculation
  maxper1 <- 99
  max01 <- quantile(prec_country[[i]],maxper1/100,na.rm=T)
  maxper2 <- 1
  max02 <- quantile(prec_country[[i]],maxper2/100,na.rm=T)
  
  # Center at MRP
  pred <- crosspred(bvar,coef=blup_country[[i]]$blup,vcov=blup_country[[i]]$vcov,
                    model.link="log",at=prec_country[[i]],cen=mmt,by=0.1)
  
  # Save results for ggplot
  sdf <- data.frame("country"=country_list[i],
                    "lprec"=as.numeric(names(pred$allRRfit)),"prec"=exp(as.numeric(names(pred$allRRfit))),
                    "rrmean"=pred$allRRfit,"rrlci"=pred$allRRlow,"rruci"=pred$allRRhigh)
  
  # Some steps for RR calculation
  # We could miss some RRs due to the predict infinity numbers, so I added these steps
  predvar2 <- as.numeric(pred$predvar)
  
  # Find the corresponding absolute value for RR prediction
  maxvar1 <- which.min(abs(predvar2-max01))
  maxvar2 <- which.min(abs(predvar2-max02))
  max1 <- predvar2[maxvar1]
  max2 <- predvar2[maxvar2]
  
  # Save to pred_list
  # pred_list[[i]] <- pred
  # names(pred_list)[i] <-  country_list[i]
  
  # Save the RR and coefficients in table
  tab_rr[i, 1:6] <- c(minper, mmt, maxper1, max1, maxper2, max2)
  
  sdf$mmt <- mmt
  sdf$low_p <- max2
  sdf$high_p <- max1
  dfknots <- rbind(dfknots,sdf)
  
  tab_rr[i, 7] <- pred$allRRfit[as.character(max1)]
  tab_rr[i, 8] <- pred$allRRlow[as.character(max1)]
  tab_rr[i, 9] <- pred$allRRhigh[as.character(max1)]
  tab_rr[i, 10] <- pred$allRRfit[as.character(max2)]
  tab_rr[i, 11] <- pred$allRRlow[as.character(max2)]
  tab_rr[i, 12] <- pred$allRRhigh[as.character(max2)]
  
  tab_rr[i, 13] <- pred$allfit[as.character(max1)]
  tab_rr[i, 14] <- pred$allse[as.character(max1)]
  tab_rr[i, 15] <- pred$allfit[as.character(max2)]
  tab_rr[i, 16] <- pred$allse[as.character(max2)]
  
  cpmod2[[i]] <- pred
  
} 

head(dfknots)



# Prepare to ggplot
meta_ind %>% dplyr::select(country, name_country, region) -> cnt_01
cnt_01 <- unique(cnt_01[c("country","name_country","region")])

dfknots %>% left_join(cnt_01, by="country") %>% 
  arrange(region, country) %>% 
  mutate(country=factor(country, levels=unique(country)),
         name_country=factor(name_country, levels=unique(name_country)))-> dfknots_02
head(dfknots_02)

# Prepare the color
RColorBrewer::display.brewer.all()
col <- brewer_pal(palette = "Set1")(8)
col2 <- brewer_pal(palette = "Set2")(8)
cnt_col <- c(rep(col2[2],2), rep(col[2],3), rep(col[3],1), rep(col[4],2))
reg_col <- c(col2[2], col[2], col[3], col[4])

# ggplot: Figure 2
p1 <- ggplot(dfknots_02, aes(x = lprec, y = rrmean, group= region, color = region)) +
  geom_line(size = 1) +
  geom_vline(aes(xintercept = mmt), color = "red", linetype = "dashed") +
  geom_ribbon(aes(ymin = rrlci, ymax = rruci, fill = region), alpha = 0.1, color=NA) +
  scale_color_manual(values = reg_col) +
  scale_fill_manual(values = reg_col) +
  facet_wrap(~ name_country,ncol=3) +
  expand_limits(y = c(min(dfknots_02$rrlci), max(dfknots_02$rrlci))) +
  scale_x_continuous(limits = c(min(dfknots_02$lprec), max(dfknots_02$lprec)), 
                     breaks = c(-5, -2, 0, 2, 4, 6, 8, 10), 
                     labels = c(-5, -2, 0, 2, 4, 6, 8, 10)) +
  scale_y_continuous(limits = c(min(dfknots_02$rrlci), max(dfknots_02$rruci))) +
  geom_hline(yintercept = 1, linetype=1, size=.4) +
  coord_cartesian(ylim = c(0.9, 2.5)) +  
  labs(x = "Natural logarithm of weekly total precipitation", y = "RR") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),   # Remove minor grid lines
        strip.text.x = element_text(size = 10, face = "bold"),
        strip.text.y = element_text(size = 10, face = "bold"),
        legend.position = "bottom")+
  guides(color = guide_legend(title = "Region"), fill = guide_legend(title = "Region"))


# save high resolution figure
# sjPlot::save_plot("ggplot_all_02.png",p1,height=25,width=25,dpi=600)



# Save RR in each country for Figure 3
tab2 <- tab_rr
head(tab_rr)
tab_rr <- as.data.frame(tab_rr)
tab_rr %>%  mutate(exp_min= round(exp(minrain),2),
                   exp_max1= round(exp(maxrain1),2),
                   exp_max2= round(exp(maxrain2),2)) -> tab_rr

# FORMAT
tab_format_rr <- cbind(name=row.names(tab_rr),
                       minperc0= tab_rr$minperc,
                       minrain0= formatC(tab_rr$minrain, format="f", big.mark=",", digits=2),
                       maxperchigh= tab_rr$maxperc1,
                       maxrainhigh= formatC(tab_rr$maxrain1, format="f", big.mark=",", digits=2),
                       maxperclow= tab_rr$maxperc2,
                       maxrainlow= formatC(tab_rr$maxrain2, format="f", big.mark=",", digits=2),
                       RRhigh = paste0(
                         formatC(tab_rr$RR1, format="f", big.mark=",", digits=2), " (",
                         formatC(tab_rr$RRlow1, format="f", big.mark=",", digits=2), "-",
                         formatC(tab_rr$RRhigh1, format="f", big.mark=",", digits=2), ")"),
                       RRlow = paste0(
                         formatC(tab_rr$RR2, format="f", big.mark=",", digits=2), " (",
                         formatC(tab_rr$RRlow2, format="f", big.mark=",", digits=2), "-",
                         formatC(tab_rr$RRhigh2, format="f", big.mark=",", digits=2), ")"),
                       tab_rr[1:19])

head(tab_format_rr)
# Excluding the RR = 1
tab_format_rr %>% mutate(RRhigh=ifelse(minperc0==99, NA,  RRhigh),
                         RR1=ifelse(minperc0==99, NA,  RR1),
                         RRlow=ifelse(minperc0==1, NA,  RRlow),
                         RR2=ifelse(minperc0==1, NA,  RR2)) -> tab_format_rr
# Save RR table
# write.csv(tab_format_rr, file="01rr_country_tab01_meta_rr_99th.csv")


# Main results at climate zone level ----------------------------------------------

# Read mete-regression predicted by climate zone
summary(meta_clim)

# Predict by climate zone predictor: Later we exclude the polar results due to limit locations
clim_01 <- c("Tropical","Arid","Temperate","Cold")
fclim <- factor(meta_ind$clim, levels=c("Tropical","Arid","Temperate","Cold"))
datanew2 <- data.frame(clim=clim_01,row.names=clim_01)
mvpred2 <- predict(meta_clim,datanew2,vcov=T,format="list")

# Split the data at climate zone level
dlist2 <- split(country_all, country_all$clim)[clim_01]
prec_clim <- list()

# Run the loop: obtain the exposure for predict
for(i in seq(levels(fclim))) {
  
  clim <- levels(fclim)[i]
  
  prec_all_obs <-unlist(sapply(dlist2[clim], function(x) x$prec), use.names=F)
  prec_all_obs <- as.data.frame(prec_all_obs)
  prec_clim[i] <- prec_all_obs
  names(prec_clim)[i] <- clim
  
}



# Predict curves and save RR
# Create the table to save the RR by climate zone
tab_rr <- matrix(NA, 5 ,16)
rownames(tab_rr) <- c("Tropical","Arid","Temperate","Cold","Polar")
colnames(tab_rr) <- c("minperc", "minrain", "maxperc1", "maxrain1", "maxperc2", "maxrain2",
                      "RR1","RRlow1","RRhigh1","RR2","RRlow2","RRhigh2","fit1","se1","fit2","se2")

# Save predict curves
cpmod2 <- list()

# For ggplot
cnm <- c("ctry","loc","knots","lprec","prec","rrmean","rrlci","rruci")
dfknots <- data.frame(matrix(NA,nrow=0,ncol=length(cnm),dimnames=list(NULL,cnm)))
dfknots <- data.frame()

# Percentiles for predict the curves
predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))

# Run the loop
for(i in seq(levels(fclim))) {
  
  # Create the variable to predict
  prec_clim_per <- quantile(prec_clim[[i]],predper/100,na.rm=T)
  predvar <- quantile(prec_clim[[i]],1:99/100,na.rm=T)
  
  # REDEFINE THE FUNCTION USING ALL THE ARGUMENTS (BOUNDARY KNOTS INCLUDED)
  argvar <- list(x=predvar,fun=varfun,
                 knots=quantile(prec_clim[[i]],varper/100,na.rm=T),
                 Bound=range(prec_clim_per[c("0.0%","100.0%")],na.rm=T))
  if(!is.null(vardegree)) argvar$degree <- vardegree
  bvar <- do.call(onebasis,argvar)
  
  # Predict the curves
  cp1 <- crosspred(bvar,coef=mvpred2[[i]]$fit,vcov=mvpred2[[i]]$vcov,model.link="log",
                   at=prec_clim[[i]],cen=quantile(prec_clim[[i]],50/100,na.rm=T))
  
  # MRP, 1st and 99th
  minper <- (1:99)[which.min((bvar%*%coef(cp1))[1:99,])]
  mmt <- quantile(prec_clim[[i]],minper/100,na.rm=T)
  maxper1 <- 99
  max01 <- quantile(prec_clim[[i]],maxper1/100,na.rm=T)
  maxper2 <- 1
  max02 <- quantile(prec_clim[[i]],maxper2/100,na.rm=T)
  
  # Center at MRP
  pred <- crosspred(bvar,coef=mvpred2[[i]]$fit,vcov=mvpred2[[i]]$vcov,
                    model.link="log",at=prec_clim[[i]],cen=mmt, by=0.1)
  
  # Save for ggplot
  sdf <- data.frame("clim"=clim_01[i],
                    "lprec"=as.numeric(names(pred$allRRfit)),"prec"=exp(as.numeric(names(pred$allRRfit))),
                    "rrmean"=pred$allRRfit,"rrlci"=pred$allRRlow,"rruci"=pred$allRRhigh)
  
  # For RR calculation
  predvar2 <- as.numeric(pred$predvar)
  
  # Find the corresponding absolute value for RR prediction
  maxvar1 <- which.min(abs(predvar2-max01))
  maxvar2 <- which.min(abs(predvar2-max02))
  max1 <- predvar2[maxvar1]
  max2 <- predvar2[maxvar2]
  
  sdf$mmt <- mmt
  sdf$low_p <- max2
  sdf$high_p <- max1
  dfknots <- rbind(dfknots,sdf)
  
  # Save the RR in table
  tab_rr[i, 1:6] <- c(minper, mmt, maxper1, max1, maxper2, max2)
  
  tab_rr[i, 7] <- pred$allRRfit[as.character(max1)]
  tab_rr[i, 8] <- pred$allRRlow[as.character(max1)]
  tab_rr[i, 9] <- pred$allRRhigh[as.character(max1)]
  tab_rr[i, 10] <- pred$allRRfit[as.character(max2)]
  tab_rr[i, 11] <- pred$allRRlow[as.character(max2)]
  tab_rr[i, 12] <- pred$allRRhigh[as.character(max2)]
  
  tab_rr[i, 13] <- pred$allfit[as.character(max1)]
  tab_rr[i, 14] <- pred$allse[as.character(max1)]
  tab_rr[i, 15] <- pred$allfit[as.character(max2)]
  tab_rr[i, 16] <- pred$allse[as.character(max2)]
  
  cpmod2[[i]] <- pred
  
} 



# Plot the curves with ggplot version
head(dfknots)
dfknots %>% 
  mutate(clim=factor(clim, levels=unique(clim_01)))-> dfknots_02

# Prepare the color
RColorBrewer::display.brewer.all()
col <- brewer_pal(palette = "Set1")(8)
col2 <- brewer_pal(palette = "Set2")(8)
clim_col <- c(col[7], col[5], col2[5], col2[3])

# ggplot: Figure 2
p2 <- ggplot(dfknots_02, aes(x = lprec, y = rrmean, group=clim, color = clim)) +
  geom_line(size = 1) +
  geom_vline(aes(xintercept = mmt), color = "red", linetype = "dashed") +
  geom_ribbon(aes(ymin = rrlci, ymax = rruci, fill = clim), alpha = 0.1, color=NA) +
  scale_color_manual(values = clim_col) +
  scale_fill_manual(values = clim_col) +
  facet_wrap(~ clim,ncol=2) +
  expand_limits(y = c(min(dfknots_02$rrlci), max(dfknots_02$rrlci))) +
  scale_x_continuous(limits = c(min(dfknots_02$lprec), max(dfknots_02$lprec)), 
                     breaks = c(-5, -2, 0, 2, 4, 6, 8, 10), 
                     labels = c(-5, -2, 0, 2, 4, 6, 8, 10)) +
  scale_y_continuous(limits = c(min(dfknots_02$rrlci), max(dfknots_02$rruci))) +
  geom_hline(yintercept = 1, linetype=1, size=.4) +
  coord_cartesian(ylim = c(0.9, 1.8)) +  
  labs(x = "Natural logarithm of weekly total precipitation", y = "RR") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),   # Remove minor grid lines
        strip.text.x = element_text(size = 15, face = "bold"),
        strip.text.y = element_text(size = 15, face = "bold"),
        legend.position = "none")

# save high resolution figure
# sjPlot::save_plot("ggplot_clim_01.png",p2,height=20,width=20,dpi=600)



# Save RR by climate zone
head(tab_rr)
tab_rr <- as.data.frame(tab_rr)
tab_rr %>%  mutate(exp_min= round(exp(minrain),2),
                   exp_max1= round(exp(maxrain1),2),
                   exp_max2= round(exp(maxrain1),2)) -> tab_rr

# FORMAT
tab_format_rr <- cbind(name=row.names(tab_rr),
                       minperc0= tab_rr$minperc,
                       minrain0= formatC(tab_rr$minrain, format="f", big.mark=",", digits=2),
                       maxperchigh= tab_rr$maxperc1,
                       maxrainhigh= formatC(tab_rr$maxrain1, format="f", big.mark=",", digits=2),
                       maxperclow= tab_rr$maxperc2,
                       maxrainlow= formatC(tab_rr$maxrain2, format="f", big.mark=",", digits=2),
                       RRhigh = paste0(
                         formatC(tab_rr$RR1, format="f", big.mark=",", digits=2), " (",
                         formatC(tab_rr$RRlow1, format="f", big.mark=",", digits=2), "-",
                         formatC(tab_rr$RRhigh1, format="f", big.mark=",", digits=2), ")"),
                       RRlow = paste0(
                         formatC(tab_rr$RR2, format="f", big.mark=",", digits=2), " (",
                         formatC(tab_rr$RRlow2, format="f", big.mark=",", digits=2), "-",
                         formatC(tab_rr$RRhigh2, format="f", big.mark=",", digits=2), ")"),
                       tab_rr[1:19])

head(tab_format_rr)
# Save RR table
# write.csv(tab_format_rr, file="01rr_clim_tab01.csv")


