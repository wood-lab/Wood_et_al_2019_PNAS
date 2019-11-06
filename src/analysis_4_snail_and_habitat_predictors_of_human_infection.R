---
title: "Analysis 4: Identifying snail- and habitat-related predictors of human urogenital schistosomiasis burden"
author: "Isabel Jones and Chelsea Wood"
date: "updated 9 July 2019"
---
  
# Introduction
  
  # We want to know which snail- and habitat-variables best predict human schistosomiasis infection rates, here measured 
  # as re-infection to children over the course of one year following praziquantel treatment. 
  
  # WHICH SITES WERE INCLUDED?
  # For this analysis, which investigated the influence of snail, snail–habitat, and habitat variables on human infection burden, 
  # we wanted to include all site–year combinations that were unmanipulated; since human infection burdens were measured annually, 
  # year was the most temporally resolved unit of time that was possible to use in these two analyses. In year 1 (2016–2017), there 
  # were 12 unmanipulated villages (i.e., 4 villages were involved in the manipulative study that began in July of 2016). In year 2 
  # (2017-2018), an additional 2 villages were dedicated to the first manipulative study, and an additional 3 villages were dedicated 
  # to the second manipulative study, which began in June of 2017; however, 1 village that had been involved in the first manipulative 
  # study (Mbakhana) in year 1 was restored to its natural state by removal of all experimental manipulations in year 2, and we therefore 
  # had a total of 8 unmanipulated villages in year 2. See SI Appendix, Figure S3 for a schematic illustrating the decision-making process
  # for how villages and water-access sites were excluded in each analysis.

# Download required packages and import data

library(tidyverse)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(ggeffects)
library(plyr)
library(AICcmodavg)
library(sjPlot)
library(broom)
library(broom.mixed)
library(MuMIn) # for averaging model predictions
library(bbmle) # AICtab for negative binomial
library(Hmisc) # for CI in ggplot
library(rgdal)
library(spdep)
library(fields)
library(cowplot)


####### Download data - note that the individual-level data on which this analysis was conducted are not provided, to protect patient privacy.

villagedata <- read_csv("data/INDIVIDUAL-LEVEL DATA WITHHELD TO PROTECT PATIENT PRIVACY")

# keep only those lines of data from 2017 or 2018 AND where no manipulation (i.e., net or vegetation removal manipulation) was in place
data_subset <- villagedata %>% 
  filter((year==2017 & net!=1) | (year==2018 & (prawn|veg_removal)!=1))

# check village names
unique(data_subset$Village)

# make LakeYN into a factor
data_subset$LakeYN=as.factor(data_subset$LakeYN)

# keep only complete cases
data_subset_complete <- data_subset[complete.cases(data_subset), ]

# get the full list of included villages
data_subset_complete %>% summarise(unique(Village))

# get the list of villages included in 2017
data_subset_complete[data_subset_complete$year==2017,] %>% summarise(unique(Village))

# get the list of villages included in 2018
data_subset_complete[data_subset_complete$year==2018,] %>% summarise(unique(Village))

# check which villages are included in each year
data_subset %>% 
  filter(year==2017) %>% 
  distinct(Village)
data_subset %>% summarise(n_distinct(Village))

# check number of kids in study
data_subset %>% summarise(n_distinct(ID))


# GENERAL MODEL FORMULATIONS

# demographic variables and random effects:
# Class + sex + Pop_sc + LakeYN + (1|Village/ID)

# 1. 'null' model, with demographic predictors and random effects only: 
# Sh ~ Class + sex + Pop_sc + LakeYN + (1|Village/ID)

# 2. snail density only 
# Sh ~ BulinusDens_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID)

# 3. snail infection prevalence only
# Sh ~ ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID)

# 4. size and snail density separate
# Sh ~ BulinusDens_sc + TotalSize_enclosure_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID)

# 5. size and snail density combine for estimate of total snails 
# Sh ~ BulinusTotal_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID)

# 6. size and snail density and prevalence separate
# Sh ~ BulinusDens_sc + TotalSize_enclosure_sc + ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID)

# 7. size and snail density and prevalence combine to estimate total number of infected snails
# Sh ~ ShBulinus_Total_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID)

# 8. density and prevalence separate
# Sh ~ BulinusDens_sc + ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID)

# 9. density and prevalence combine to estimate infected snail density
# Sh ~ ShBulinusDens_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID)

# 10. size and percent other veg and percent mud separate
# Sh ~ TotalSize_enclosure_sc + PercOther_sc + PercMud_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID)

# 11. size and percent other veg and size and percent mud combine to estimate total area of other veg and total area of mud
# Sh ~ OtherVegTotal_sc + MudTotal_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID)

# 12. size and percent other and percent mud and average mass of veg sampled separate 
# Sh ~ TotalSize_enclosure_sc + VegMassAvg_sampled_sc + PercOther_sc + PercMud_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID)

# 13. size and percent other and average mass of veg sampled separate (+ percent mud) combine to estimate total mass of other vegetation 
# Sh ~ VegMassTotal_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID)

# 14. full orthogonal model 
# Sh ~ TotalSize_enclosure_sc + VegMassAvg_sampled_sc + BulinusDens_sc + ShPrev_Bulinus_sc + PercOther_sc + PercMud_sc +  Class + sex + Pop_sc + LakeYN + (1|Village/ID)



# LOGISTIC MODELS

########################################################
###########  S haematobium logistic models   ########### 
# every model: Class + sex + Pop_sc + LakeYN + (1|Village/ID)

# 1.
null.model <- glmer(Sh ~ Class + sex + Pop_sc + LakeYN + (1|Village/ID), data = data_subset_complete, 
                    family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(null.model)
data_subset_complete$model1_resids<-residuals(null.model)

# 2.
density.only <- glmer(Sh ~ BulinusDens_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), data = data_subset_complete, 
                      family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(density.only)
data_subset_complete$model2_resids<-residuals(density.only)

# 3.
prev.only <- glmer(Sh ~ ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), data = data_subset_complete, 
                   family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(prev.only)
data_subset_complete$model3_resids<-residuals(prev.only)

# 4. 
size.dens <- glmer(Sh ~ BulinusDens_sc + TotalSize_enclosure_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), data = data_subset_complete, 
                   family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(size.dens)
data_subset_complete$model4_resids<-residuals(size.dens)

# 5.
snail.total <- glmer(Sh ~ BulinusTotal_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), data = data_subset_complete, 
                     family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(snail.total)
data_subset_complete$model5_resids<-residuals(snail.total)

# 6.
size.dens.prev <- glmer(Sh ~ BulinusDens_sc + TotalSize_enclosure_sc + ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), data = data_subset_complete, 
                        family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(size.dens.prev)
data_subset_complete$model6_resids<-residuals(size.dens.prev)

# 7.
inf.snail.total <- glmer(Sh ~ ShBulinus_Total_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), data = data_subset_complete, 
                         family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(inf.snail.total)
data_subset_complete$model7_resids<-residuals(inf.snail.total)

# 8. BulinusDens_sc + ShPrev_Bulinus_sc = ShBulinusDens_sc
dens.prev <- glmer(Sh ~ BulinusDens_sc + ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), data = data_subset_complete, 
                   family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(dens.prev)
data_subset_complete$model8_resids<-residuals(dens.prev)

# 9. 
inf.snail.dens <- glmer(Sh ~ ShBulinusDens_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), data = data_subset_complete, 
                        family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(inf.snail.dens)
data_subset_complete$model9_resids<-residuals(inf.snail.dens)

# 10
size.perc.other <- glmer(Sh ~ TotalSize_enclosure_sc + PercOther_sc + PercMud_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), data = data_subset_complete, 
                         family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(size.perc.other)
data_subset_complete$model10_resids<-residuals(size.perc.other)

# 11. 
total.other <- glmer(Sh ~ OtherVegTotal_sc + MudTotal_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), data = data_subset_complete, 
                     family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(total.other)
data_subset_complete$model11_resids<-residuals(total.other)

# 12. TotalSize_enclosure_sc + PercOther_sc + PercMud_sc + VegMassAvg_sampled_sc = VegMassTotal_sc
size.perc.vegmass <- glmer(Sh ~ TotalSize_enclosure_sc + VegMassAvg_sampled_sc + PercOther_sc + PercMud_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), data = data_subset_complete, 
                           family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(size.perc.vegmass)
data_subset_complete$model12_resids<-residuals(size.perc.vegmass)

# 13.
veg.total.mass <- glmer(Sh ~ VegMassTotal_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), data = data_subset_complete, 
                        family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(veg.total.mass)
data_subset_complete$model13_resids<-residuals(veg.total.mass)

## full orthogonal model
# 14. TotalSize_enclosure_sc + PercOther_sc + VegMassAvg_sampled_sc + BulinusDens_sc + ShPrev_Bulinus_sc
full.orthog <- glmer(Sh ~ TotalSize_enclosure_sc + VegMassAvg_sampled_sc + BulinusDens_sc + ShPrev_Bulinus_sc + PercOther_sc + PercMud_sc +  Class + sex + Pop_sc + LakeYN + (1|Village/ID), data = data_subset_complete, 
                     family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(full.orthog)
data_subset_complete$model14_resids<-residuals(full.orthog)



### test for spatial auto-correlation in these models
# start by making a matrix of villages and residuals for each observation

residual_matrix<-cbind.data.frame(data_subset_complete$VillageCode,as.numeric(data_subset_complete$model1_resids),
                                  as.numeric(data_subset_complete$model2_resids),
                                  as.numeric(data_subset_complete$model3_resids),as.numeric(data_subset_complete$model4_resids),
                                  as.numeric(data_subset_complete$model5_resids),as.numeric(data_subset_complete$model6_resids),
                                  as.numeric(data_subset_complete$model7_resids),as.numeric(data_subset_complete$model8_resids),
                                  as.numeric(data_subset_complete$model9_resids),as.numeric(data_subset_complete$model10_resids),
                                  as.numeric(data_subset_complete$model11_resids),as.numeric(data_subset_complete$model12_resids),
                                  as.numeric(data_subset_complete$model13_resids),as.numeric(data_subset_complete$model14_resids))
residual_matrix<-as.data.frame(residual_matrix)
names(residual_matrix)<-c("village","model1","model2","model3","model4","model5","model6","model7","model8","model9",
                          "model10","model11","model12","model13","model14")

#Link village codes with village names

code_names<-c("Mbane","Malla Tack","Guidik","Syer","Merina Gewel","Diokoul","Malla","Diokhor",
              "Ndiol Maure","Lampsar","Mbarigot","Mbakhana","Ndiawdoune")
codes<-c("ME","MT","GK","ST","MG","DF","MA","DT","NM","LR","MO","MB","NE")
decode<-cbind(code_names,codes)
decode<-as.data.frame(decode)
names(decode)<-c("village_name","village")

decoded_resid_matrix<-merge(decode,residual_matrix,by=c("village"))
names(decoded_resid_matrix)<-c("village_code","village","model1","model2","model3","model4","model5","model6","model7",
                               "model8","model9","model10","model11","model12","model13","model14")
decoded_resid_matrix<-as.data.frame(decoded_resid_matrix)

aggregated<-aggregate(cbind(decoded_resid_matrix$model1,decoded_resid_matrix$model2,decoded_resid_matrix$model3,
                            decoded_resid_matrix$model4,decoded_resid_matrix$model5,decoded_resid_matrix$model6,
                            decoded_resid_matrix$model7,decoded_resid_matrix$model8,decoded_resid_matrix$model9,
                            decoded_resid_matrix$model10,decoded_resid_matrix$model11,decoded_resid_matrix$model12,
                            decoded_resid_matrix$model13,decoded_resid_matrix$model14),
                      by=list(decoded_resid_matrix$village),FUN=mean)
names(aggregated)<-c("village","model1","model2","model3","model4","model5","model6","model7","model8","model9",
                     "model10","model11","model12","model13","model14")

# based on recommendation here: https://stats.stackexchange.com/questions/40667/testing-for-spatial-autocorrelation-in-a-negative-binomial-regression-model

# need to make a listw using the nb2listw function
# start by importing all of the village-level polygons
village_polygons<-readOGR("data/polygons.shp")
str(village_polygons)
village_polygons$Name

# with help from this resource, we calculated the centroid of each polygon - not the average of the coordinates in 
# each dimension, but weights the component triangles of the polygon by area:
# https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf

coords<-coordinates(village_polygons)
IDs<-c("Diokhor","Diokoul","Guidik","Lampsar","Malla","Malla Tack","Mbakhana",
       "Mbane","Mbarigot","Merina Gewel","Ndiawdoune","Ndiol Maure","Syer")
village_nb<-tri2nb(coords,row.names=IDs)
class(village_nb)
length(village_nb)

# now we need to use the polygon centroids to create a listw object that will tell the Moran's i calculation how to
# weight the villages against one another. We want them to be weighted by distance.

dlist <- nbdists(village_nb, coords, longlat=T) 
dlist <- lapply(dlist, function(x) 1/x)

lw<-nb2listw(village_nb,dlist,zero.policy = T)
length(lw)
lw$weights

moran.mc(aggregated$model1,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model2,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model3,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model4,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model5,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model6,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model7,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model8,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model9,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model10,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model11,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model12,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model13,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model14,lw,na.action = na.omit,nsim=999)


# okay, there seems to be some spatial autocorrelation. What if we lumped sites that were both close together and
# similar in terms of their residuals? We would lose power, but could reduce auto-correlation. Let's proceed by lumping
# the closest sites.


coords<-coordinates(village_polygons)
IDs<-c("Diokhor","Diokoul","Guidik","Lampsar","Malla","Malla Tack","Mbakhana",
       "Mbane","Mbarigot","Merina Gewel","Ndiawdoune","Ndiol Maure","Syer")
spatialmatrix<-cbind.data.frame(coords)

dist<-rdist.earth(spatialmatrix,spatialmatrix,miles=FALSE)
distances<-cbind.data.frame(IDs,dist)
names(distances)<-c("","Diokhor","Diokoul","Guidik","Lampsar","Malla","Malla Tack","Mbakhana",
                    "Mbane","Mbarigot","Merina Gewel","Ndiawdoune","Ndiol Maure","Syer")

# start with Mbakhana and Mbarigot

data_subset_complete$mergevillage<-gsub("Mbakhana|Mbarigot","MbakhanaMbarigot",data_subset_complete$Village)
write.csv(data_subset_complete$mergevillage)


# Now that we've got Mbakhana and Mbarigot merged into a single village, we can proceed with analysis

########################################################
###########  S haematobium logistic models   ########### 
# every model: Class + sex + Pop_sc + LakeYN + (1|Village/ID)

# 1.
null.model <- glmer(Sh ~ Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
                    family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(null.model)
data_subset_complete$model1_resids<-residuals(null.model)

predictions<-predict(null.model,type="response")
data<-cbind(predictions,data_subset_complete$Sh)
data<-as.data.frame(data)

a_ii<-ggplot()+
  geom_jitter(data=data,aes(x=data$V2,y=data$predictions))+
  theme_minimal()+
  xlab('')+
  ylab('')+
  xlim(-0.5,1.5)+scale_x_continuous(breaks=c(0,1))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20,face="bold"))

# 2.
density.only <- glmer(Sh ~ BulinusDens_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
                      family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(density.only)
data_subset_complete$model2_resids<-residuals(density.only)

# 3.
prev.only <- glmer(Sh ~ ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
                   family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(prev.only)
data_subset_complete$model3_resids<-residuals(prev.only)

# 4. 
size.dens <- glmer(Sh ~ BulinusDens_sc + TotalSize_enclosure_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
                   family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(size.dens)
data_subset_complete$model4_resids<-residuals(size.dens)

# 5.
snail.total <- glmer(Sh ~ BulinusTotal_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
                     family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(snail.total)
data_subset_complete$model5_resids<-residuals(snail.total)

predictions<-predict(snail.total,type="response")
data<-cbind(predictions,data_subset_complete$Sh)
data<-as.data.frame(data)

a_iv<-ggplot()+
  geom_jitter(data=data,aes(x=data$V2,y=data$predictions))+
  theme_minimal()+
  xlab('')+
  ylab('')+
  xlim(-0.5,1.5)+scale_x_continuous(breaks=c(0,1))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20,face="bold"))

# 6.
size.dens.prev <- glmer(Sh ~ BulinusDens_sc + TotalSize_enclosure_sc + ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
                        family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(size.dens.prev)
data_subset_complete$model6_resids<-residuals(size.dens.prev)

# 7.
inf.snail.total <- glmer(Sh ~ ShBulinus_Total_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
                         family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(inf.snail.total)
data_subset_complete$model7_resids<-residuals(inf.snail.total)

predictions<-predict(inf.snail.total,type="response")
data<-cbind(predictions,data_subset_complete$Sh)
data<-as.data.frame(data)

a_vi<-ggplot()+
  geom_jitter(data=data,aes(x=data$V2,y=data$predictions))+
  theme_minimal()+
  xlab('')+
  ylab('')+
  xlim(-0.5,1.5)+scale_x_continuous(breaks=c(0,1))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20,face="bold"))


# 8. BulinusDens_sc + ShPrev_Bulinus_sc = ShBulinusDens_sc
dens.prev <- glmer(Sh ~ BulinusDens_sc + ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
                   family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(dens.prev)
data_subset_complete$model8_resids<-residuals(dens.prev)

# 9. 
inf.snail.dens <- glmer(Sh ~ ShBulinusDens_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
                        family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(inf.snail.dens)
data_subset_complete$model9_resids<-residuals(inf.snail.dens)

# 10
size.perc.other <- glmer(Sh ~ TotalSize_enclosure_sc + PercOther_sc + PercMud_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
                         family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(size.perc.other)
data_subset_complete$model10_resids<-residuals(size.perc.other)

predictions<-predict(size.perc.other,type="response")
data<-cbind(predictions,data_subset_complete$Sh)
data<-as.data.frame(data)

a_iii<-ggplot()+
  geom_jitter(data=data,aes(x=data$V2,y=data$predictions))+
  theme_minimal()+
  xlab('')+
  ylab('')+
  xlim(-0.5,1.5)+scale_x_continuous(breaks=c(0,1))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20,face="bold"))

# 11. 
total.other <- glmer(Sh ~ OtherVegTotal_sc + MudTotal_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
                     family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(total.other)
data_subset_complete$model11_resids<-residuals(total.other)

predictions<-predict(total.other,type="response")
data<-cbind(predictions,data_subset_complete$Sh)
data<-as.data.frame(data)

a_i<-ggplot()+
  geom_jitter(data=data,aes(x=data$V2,y=data$predictions))+
  theme_minimal()+
  xlab('')+
  ylab('')+
  xlim(-0.5,1.5)+scale_x_continuous(breaks=c(0,1))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20,face="bold"))

# 12. TotalSize_enclosure_sc + PercOther_sc + PercMud_sc + VegMassAvg_sampled_sc = VegMassTotal_sc
size.perc.vegmass <- glmer(Sh ~ TotalSize_enclosure_sc + VegMassAvg_sampled_sc + PercOther_sc + PercMud_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
                           family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(size.perc.vegmass)
data_subset_complete$model12_resids<-residuals(size.perc.vegmass)

# 13.
veg.total.mass <- glmer(Sh ~ VegMassTotal_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
                        family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(veg.total.mass)
data_subset_complete$model13_resids<-residuals(veg.total.mass)

predictions<-predict(veg.total.mass,type="response")
data<-cbind(predictions,data_subset_complete$Sh)
data<-as.data.frame(data)

a_v<-ggplot()+
  geom_jitter(data=data,aes(x=data$V2,y=data$predictions))+
  theme_minimal()+
  xlab('')+
  ylab('')+
  xlim(-0.5,1.5)+scale_x_continuous(breaks=c(0,1))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20,face="bold"))

## full orthogonal model
# 14. TotalSize_enclosure_sc + PercOther_sc + VegMassAvg_sampled_sc + BulinusDens_sc + ShPrev_Bulinus_sc
full.orthog <- glmer(Sh ~ TotalSize_enclosure_sc + VegMassAvg_sampled_sc + BulinusDens_sc + ShPrev_Bulinus_sc + PercOther_sc + PercMud_sc +  Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
                     family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(full.orthog)
data_subset_complete$model14_resids<-residuals(full.orthog)


#make a plot of observed versus predicted for all interpreted models

a<-ggdraw(plot=NULL,xlim=c(0,6),ylim=c(0,1))+
  draw_plot(a_i,x=0,y=0,width=1,height=1)+
  draw_plot(a_ii,x=1,y=0,width=1,height=1)+
  draw_plot(a_iii,x=2,y=0,width=1,height=1)+
  draw_plot(a_iv,x=3,y=0,width=1,height=1)+
  draw_plot(a_v,x=4,y=0,width=1,height=1)+
  draw_plot(a_vi,x=5,y=0,width=1,height=1)+
  draw_text("predicted",x=0.09,y=0.45,size=18,hjust=0,vjust=0,angle=90)+
  draw_text("(1)",x=0.35,y=0.97,size=24)+
  draw_text("(2)",x=1.35,y=0.97,size=24)+
  draw_text("(3)",x=2.35,y=0.97,size=24)+
  draw_text("(4)",x=3.35,y=0.97,size=24)+
  draw_text("(5)",x=4.35,y=0.97,size=24)+
  draw_text("(6)",x=5.35,y=0.97,size=24)


### test for spatial autocorrelation in these models
# start by making a matrix of villages and residuals for each observation

residual_matrix<-cbind.data.frame(data_subset_complete$mergevillage,as.numeric(data_subset_complete$model1_resids),
                                  as.numeric(data_subset_complete$model2_resids),
                                  as.numeric(data_subset_complete$model3_resids),as.numeric(data_subset_complete$model4_resids),
                                  as.numeric(data_subset_complete$model5_resids),as.numeric(data_subset_complete$model6_resids),
                                  as.numeric(data_subset_complete$model7_resids),as.numeric(data_subset_complete$model8_resids),
                                  as.numeric(data_subset_complete$model9_resids),as.numeric(data_subset_complete$model10_resids),
                                  as.numeric(data_subset_complete$model11_resids),as.numeric(data_subset_complete$model12_resids),
                                  as.numeric(data_subset_complete$model13_resids),as.numeric(data_subset_complete$model14_resids))
residual_matrix<-as.data.frame(residual_matrix)
names(residual_matrix)<-c("village","model1","model2","model3","model4","model5","model6","model7","model8","model9",
                          "model10","model11","model12","model13","model14")

aggregated<-aggregate(cbind(residual_matrix$model1,residual_matrix$model2,residual_matrix$model3,
                            residual_matrix$model4,residual_matrix$model5,residual_matrix$model6,
                            residual_matrix$model7,residual_matrix$model8,residual_matrix$model9,
                            residual_matrix$model10,residual_matrix$model11,residual_matrix$model12,
                            residual_matrix$model13,residual_matrix$model14),
                      by=list(residual_matrix$village),FUN=mean)
names(aggregated)<-c("village","model1","model2","model3","model4","model5","model6","model7","model8","model9",
                     "model10","model11","model12","model13","model14")

# based on recommendation here: https://stats.stackexchange.com/questions/40667/testing-for-spatial-autocorrelation-in-a-negative-binomial-regression-model

# need to make a listw using the nb2listw function
# start by importing all of the village-level polygons

village_polygons<-readOGR("data/sitesmerged-polygons.shp")
str(village_polygons)
village_polygons$Name

# With help from this resource, we calculated the centroid of each polygon - not the average of the coordinates in 
# each dimension, but weights the component triangles of the polygon by area:
# https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf

coords<-coordinates(village_polygons)
IDs<-c("Diokhor","Diokoul","Guidik","Lampsar","Malla","Malla Tack","MbakhanaMbarigot",
       "Mbane","Merina Gewel","Ndiawdoune","Ndiol Maure","Syer")
village_nb<-tri2nb(coords,row.names=IDs)
class(village_nb)
length(village_nb)

# Now we need to use the polygon centroids to create a listw object that will tell the Moran's i calculation how to
# weight the villages against one another. We want them to be weighted by distance.

dlist <- nbdists(village_nb, coords, longlat=T) 
dlist <- lapply(dlist, function(x) 1/x)

lw<-nb2listw(village_nb,dlist,style="W",zero.policy = T)
length(lw)
lw$weights

moran.mc(aggregated$model1,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model2,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model3,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model4,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model5,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model6,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model7,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model8,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model9,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model10,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model11,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model12,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model13,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated$model14,lw,na.action = na.omit,nsim=999)

# Okay. Just by merging the two closest sites (Mbakhana and Mbarigot), We have removed the problem with spatial
# auto-correlation.


# model selection (BIC)
Sh.mod.list <- c(null.model, density.only, prev.only, size.dens, snail.total, size.dens.prev, inf.snail.total, dens.prev, inf.snail.dens, size.perc.other, total.other, size.perc.vegmass, veg.total.mass, full.orthog)
Sh.mod.names <- c('Null', 'Sh ~ BulinusDens_sc', 'Sh ~ ShPrev_Bulinus_sc', 
                  'Sh ~ BulinusDens_sc + TotalSize_enclosure_sc',
                  'Sh ~ BulinusTotal_sc',
                  'Sh ~ BulinusDens_sc + TotalSize_enclosure_sc + ShPrev_Bulinus_sc',
                  'Sh ~ ShBulinus_Total_sc',
                  'Sh ~ BulinusDens_sc + ShPrev_Bulinus_sc',
                  'Sh ~ ShBulinusDens_sc',
                  'Sh ~ TotalSize_enclosure_sc + PercOther_sc + PercMud_sc',
                  'Sh ~ OtherVegTotal_sc + MudTotal_sc',
                  'Sh ~ TotalSize_enclosure_sc + VegMassAvg_sampled_sc + PercOther_sc + PercMud_sc',
                  'Sh ~ VegMassTotal_sc',
                  'Sh ~ TotalSize_enclosure_sc + VegMassAvg_sampled_sc + BulinusDens_sc + ShPrev_Bulinus_sc + PercOther_sc + PercMud_sc')
Sh.mod.names2 <- c('null.model', 'density.only', 'prev.only', 'size.dens', 'snail.total', 'size.dens.prev', 'inf.snail.total', 'dens.prev', 'inf.snail.dens', 'size.perc.other', 'total.other', 'size.perc.vegmass', 'veg.total.mass', 'full.orthog')

# BIC and MSE for all models is calculated in a separate script.
# For logistic/binomial models, "Analysis 4: Identifying snail- and habitat-related predictors of human urogenital schistosomiasis burden - Calculating BIC and MSE for logistic/binomial models"
# We used the BICs derived in those scripts to rank the models.

## top 6 models, within 10 deltaBIC
Sh.mod.subset <- c(total.other, null.model, size.perc.other, snail.total, veg.total.mass,inf.snail.total)
Sh.mod.subset_names <- c('total.other', 'null.model', 'size.perc.other', 'snail.total', 'veg.total.mass', 'inf.snail.total')


#get the marginal and conditional R2 for all models
tab_model(Sh.mod.list, show.r2 = TRUE, show.icc = FALSE, show.aic = FALSE,
          dv.labels = c('Null', 'Sh ~ BulinusDens_sc', 'Sh ~ ShPrev_Bulinus_sc', 
                        'Sh ~ BulinusDens_sc + TotalSize_enclosure_sc',
                        'Sh ~ BulinusTotal_sc',
                        'Sh ~ BulinusDens_sc + TotalSize_enclosure_sc + ShPrev_Bulinus_sc',
                        'Sh ~ ShBulinus_Total_sc',
                        'Sh ~ BulinusDens_sc + ShPrev_Bulinus_sc',
                        'Sh ~ ShBulinusDens_sc',
                        'Sh ~ TotalSize_enclosure_sc + PercOther_sc + PercMud_sc',
                        'Sh ~ OtherVegTotal_sc + MudTotal_sc',
                        'Sh ~ TotalSize_enclosure_sc + VegMassAvg_sampled_sc + PercOther_sc + PercMud_sc',
                        'Sh ~ VegMassTotal_sc',
                        'Sh ~ TotalSize_enclosure_sc + VegMassAvg_sampled_sc + BulinusDens_sc + ShPrev_Bulinus_sc + PercOther_sc + PercMud_sc'))





# top model coefficient plot (main text, Figure 3A) for logistic model

# get predictions for each model
M1.preds <- tidy(Sh.mod.subset[[1]], conf.int=TRUE) %>% 
  mutate(Model="Model 1",
         estimate = exp(estimate),
         conf.low = exp(conf.low),
         conf.high = exp(conf.high))
M2.preds <- tidy(Sh.mod.subset[[2]], conf.int=TRUE) %>% 
  mutate(Model="Model 2",
         estimate = exp(estimate),
         conf.low = exp(conf.low),
         conf.high = exp(conf.high))
M3.preds <- tidy(Sh.mod.subset[[3]], conf.int=TRUE) %>% 
  mutate(Model="Model 3",
         estimate = exp(estimate),
         conf.low = exp(conf.low),
         conf.high = exp(conf.high))
M4.preds <- tidy(Sh.mod.subset[[4]], conf.int=TRUE) %>% 
  mutate(Model="Model 4",
         estimate = exp(estimate),
         conf.low = exp(conf.low),
         conf.high = exp(conf.high))
M5.preds <- tidy(Sh.mod.subset[[5]], conf.int=TRUE) %>% 
  mutate(Model="Model 5",
         estimate = exp(estimate),
         conf.low = exp(conf.low),
         conf.high = exp(conf.high))
M6.preds <- tidy(Sh.mod.subset[[6]], conf.int=TRUE) %>% 
  mutate(Model="Model 6",
         estimate = exp(estimate),
         conf.low = exp(conf.low),
         conf.high = exp(conf.high))

Sh.odds = bind_rows(M1.preds,M2.preds,M3.preds,M4.preds,M5.preds,M6.preds) %>% 
  filter(term != '(Intercept)')

dodger = position_dodge(width = 0.9)
# elements like pointrange and position_dodge only work when the outcome
# is mapped to y, need to go through with OR set as y then flip at the end
ggplot(Sh.odds[Sh.odds$group=='fixed',], aes(y = estimate, x = term, color=term, shape = Model)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),position = dodger,size = 0.4) +
  scale_colour_manual(values=c("black","black","black","black","black","black","black","black","black","black","black","black")) +
  scale_x_discrete(limits=c("ShBulinus_Total_sc","BulinusTotal_sc","VegMassTotal_sc","PercMud_sc","PercOther_sc","MudTotal_sc",
                            "OtherVegTotal_sc","TotalSize_enclosure_sc","Class","Pop_sc","sexM","LakeYN1"),
                        labels=c("infected snail abundance",
                                                                "snail abundance","total mass of non-emergent vegetation",
                                                                "percent cover of mud",
                                                                "percent cover of non-emergent vegetation",
                                                                "area of mud",
                                                                "area of non-emergent vegetation","site area",
                                                                "school grade",
                                                                "village population","sex: male (vs. female)",
                                                                "location: lake (vs. river)")) +
  geom_hline(yintercept = 1.0, linetype = "dotted", size = 0.5) +
  scale_y_log10(breaks = c(0.25, 0.5, 1.0, 2.0, 5.0, 10, 20),minor_breaks = NULL) +
  labs(y = "odds ratio", x = "predictor") +
  theme(legend.position = "none") +
  coord_flip(ylim = c(0.25, 30)) +
  theme_classic() 

# display top model results 
tab_model(total.other, null.model, size.perc.other, snail.total, veg.total.mass,inf.snail.total,
          transform = "exp",
          pred.labels = c(),
          dv.labels = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6"),
          CSS = list(css.centeralign = 'text-align: left;'), 
          collapse.ci = TRUE,
          show.r2 = FALSE, show.icc = FALSE,
          show.aic = FALSE)



# NEGATIVE BINOMIAL MODELS 

#################       S haematobium negative binomial models     ##################  
#####################################################################################

# using glmmTMB with negative binomial with variance structure of (variance = µ(1 + µ/k): Hardin and Hilbe (2007)). 
# every model: PercMud_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID)

# 1.
Sh.nb.null.model <- glmmTMB(ShW ~ Class + sex + Pop_sc + LakeYN + (1|Village/ID), family=nbinom2, data = data_subset_complete)
summary(Sh.nb.null.model)
data_subset_complete$nb_model1_resids<-residuals(Sh.nb.null.model)

# 2.
Sh.nb.density.only <- glmmTMB(ShW ~ BulinusDens_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), family=nbinom2, data = data_subset_complete)
summary(Sh.nb.density.only)
data_subset_complete$nb_model2_resids<-residuals(Sh.nb.density.only)

# 3.
Sh.nb.prev.only <- glmmTMB(ShW ~ ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), family=nbinom2, data = data_subset_complete)
summary(Sh.nb.prev.only)
data_subset_complete$nb_model3_resids<-residuals(Sh.nb.prev.only)

# 4. 
Sh.nb.size.dens <- glmmTMB(ShW ~ BulinusDens_sc + TotalSize_enclosure_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), family=nbinom2, data = data_subset_complete)
summary(Sh.nb.size.dens)
data_subset_complete$nb_model4_resids<-residuals(Sh.nb.size.dens)

# 5.
Sh.nb.snail.total <- glmmTMB(ShW ~ BulinusTotal_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), 
                             family=nbinom2, data = data_subset_complete)
summary(Sh.nb.snail.total)
data_subset_complete$nb_model5_resids<-residuals(Sh.nb.snail.total)

# 6. 
Sh.nb.size.dens.prev <- glmmTMB(ShW ~ BulinusDens_sc + TotalSize_enclosure_sc + ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), 
                                family=nbinom2, data = data_subset_complete)
summary(Sh.nb.size.dens.prev)
data_subset_complete$nb_model6_resids<-residuals(Sh.nb.size.dens.prev)

# 7.
Sh.nb.inf.snail.total <- glmmTMB(ShW ~ ShBulinus_Total_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), 
                                 family=nbinom2, data = data_subset_complete,na.action=na.omit)
summary(Sh.nb.inf.snail.total)
data_subset_complete$nb_model7_resids<-residuals(Sh.nb.inf.snail.total)

# 8. BulinusDens_sc + ShPrev_Bulinus_sc = ShBulinusDens_sc
Sh.nb.dens.prev <- glmmTMB(ShW ~ BulinusDens_sc + ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), 
                           family=nbinom2, data = data_subset_complete)
summary(Sh.nb.dens.prev)
data_subset_complete$nb_model8_resids<-residuals(Sh.nb.dens.prev)

# 9.
Sh.nb.inf.snail.dens <- glmmTMB(ShW ~ ShBulinusDens_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), 
                                family=nbinom2, data = data_subset_complete)
summary(Sh.nb.inf.snail.dens)
data_subset_complete$nb_model9_resids<-residuals(Sh.nb.inf.snail.dens)

# 10.
Sh.nb.size.perc.other <- glmmTMB(ShW ~ TotalSize_enclosure_sc + PercOther_sc + PercMud_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), 
                                 family=nbinom2, data = data_subset_complete)
summary(Sh.nb.size.perc.other)
data_subset_complete$nb_model10_resids<-residuals(Sh.nb.size.perc.other)

# 11. 
Sh.nb.total.other <- glmmTMB(ShW ~ OtherVegTotal_sc + MudTotal_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), 
                             family=nbinom2, data = data_subset_complete)
summary(Sh.nb.total.other)
data_subset_complete$nb_model11_resids<-residuals(Sh.nb.total.other)

# 12. 
Sh.nb.size.perc.vegmass <- glmmTMB(ShW ~ TotalSize_enclosure_sc + PercOther_sc + PercMud_sc + VegMassAvg_sampled_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), 
                                   family=nbinom2, data = data_subset_complete)
summary(Sh.nb.size.perc.vegmass)
data_subset_complete$nb_model12_resids<-residuals(Sh.nb.size.perc.vegmass)

# 13. 
Sh.nb.veg.total.mass <- glmmTMB(ShW ~ VegMassTotal_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), 
                                family=nbinom2, 
                                data = data_subset_complete)
summary(Sh.nb.veg.total.mass)
data_subset_complete$nb_model13_resids<-residuals(Sh.nb.veg.total.mass)

## full orthogonal model
# 14. TotalSize_enclosure_sc + PercOther_sc + PercMud_sc + VegMassAvg_sampled_sc + BulinusDens_sc + ShPrev_Bulinus_sc
Sh.nb.full.orthog <- glmmTMB(ShW ~ TotalSize_enclosure_sc + PercOther_sc + PercMud_sc + VegMassAvg_sampled_sc + BulinusDens_sc + ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), 
                             family=nbinom2,
                             data = data_subset_complete)
summary(Sh.nb.full.orthog)
data_subset_complete$nb_model14_resids<-residuals(Sh.nb.full.orthog)


### test for spatial autocorrelation in these models
# start by making a matrix of villages and residuals for each observation

residual_matrix_nb<-cbind.data.frame(data_subset_complete$Village,as.numeric(data_subset_complete$nb_model1_resids),
                                     as.numeric(data_subset_complete$nb_model2_resids),
                                     as.numeric(data_subset_complete$nb_model3_resids),as.numeric(data_subset_complete$nb_model4_resids),
                                     as.numeric(data_subset_complete$nb_model5_resids),as.numeric(data_subset_complete$nb_model6_resids),
                                     as.numeric(data_subset_complete$nb_model7_resids),as.numeric(data_subset_complete$nb_model8_resids),
                                     as.numeric(data_subset_complete$nb_model9_resids),as.numeric(data_subset_complete$nb_model10_resids),
                                     as.numeric(data_subset_complete$nb_model11_resids),as.numeric(data_subset_complete$nb_model12_resids),
                                     as.numeric(data_subset_complete$nb_model13_resids),as.numeric(data_subset_complete$nb_model14_resids))
names(residual_matrix_nb)<-c("village","model1","model2","model3","model4","model5","model6","model7","model8","model9",
                             "model10","model11","model12","model13","model14")
residual_matrix_nb<-as.data.frame(residual_matrix_nb)

aggregated_nb<-aggregate(cbind(residual_matrix_nb$model1,residual_matrix_nb$model2,residual_matrix_nb$model3,
                               residual_matrix_nb$model4,residual_matrix_nb$model5,residual_matrix_nb$model6,
                               residual_matrix_nb$model7,residual_matrix_nb$model8,residual_matrix_nb$model9,
                               residual_matrix_nb$model10,residual_matrix_nb$model11,residual_matrix_nb$model12,
                               residual_matrix_nb$model13,residual_matrix_nb$model14),
                         by=list(residual_matrix_nb$village),function(x) mean(x,na.rm=TRUE))
names(aggregated_nb)<-c("village","model1","model2","model3","model4","model5","model6","model7","model8","model9",
                        "model10","model11","model12","model13","model14")


# based on recommendation here: https://stats.stackexchange.com/questions/40667/testing-for-spatial-autocorrelation-in-a-negative-binomial-regression-model

# need to make a listw using the nb2listw function
# start by importing all of the village-level polygons

village_polygons<-readOGR("data/polygons.shp")
str(village_polygons)
village_polygons$Name

# with help from this resource, we calculated the centroid of each polygon - not the average of the coordinates in 
# each dimension, but weights the component triangles of the polygon by area:
# https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf

coords<-coordinates(village_polygons)
IDs<-c("Diokhor","Diokoul","Guidik","Lampsar","Malla","Malla Tack","Mbakhana",
       "Mbane","Mbarigot","Merina Gewel","Ndiawdoune","Ndiol Maure","Syer")
village_nb<-tri2nb(coords,row.names=IDs)
class(village_nb)
length(village_nb)

# now we need to use the polygon centroids to create a listw object that will tell the Moran's i calculation how to
# weight the villages against one another. I want them to be weighted by distance.

dlist <- nbdists(village_nb, coords, longlat=T) 
dlist <- lapply(dlist, function(x) 1/x)

lw<-nb2listw(village_nb,dlist,zero.policy = T)
length(lw)
lw$weights

moran.mc(aggregated_nb$model1,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb$model2,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb$model3,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb$model4,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb$model5,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb$model6,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb$model7,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb$model8,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb$model9,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb$model10,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb$model11,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb$model12,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb$model13,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb$model14,lw,na.action = na.omit,nsim=999)


# okay, there isn't any significant spatial auto-correlation here, but it's close. To reduce autocorrelation and keep
# consistent with logistic analysis above, let's lump Mbakhana and Mbarigot agian.

# Now that we've got Mbakhana and Mbarigot merged into a single village, we can proceed with analysis
# negative binomial models with Mbakhana and Mbarigot lumped (mergevillage instead of village)

#################       S haematobium negative binomial models     #################       
# using glmmTMB with negative binomial with variance structure of (variance = µ(1 + µ/k): Hardin and Hilbe (2007)). 
# every model: PercMud_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID)

# 1.
Sh.nb.null.model <- glmmTMB(ShW ~ Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), family=nbinom2, data = data_subset_complete)
summary(Sh.nb.null.model)
data_subset_complete$nb_model1_resids<-residuals(Sh.nb.null.model)

# 2.
Sh.nb.density.only <- glmmTMB(ShW ~ BulinusDens_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), family=nbinom2, data = data_subset_complete)
summary(Sh.nb.density.only)
data_subset_complete$nb_model2_resids<-residuals(Sh.nb.density.only)

# 3.
Sh.nb.prev.only <- glmmTMB(ShW ~ ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), family=nbinom2, data = data_subset_complete)
summary(Sh.nb.prev.only)
data_subset_complete$nb_model3_resids<-residuals(Sh.nb.prev.only)
 
# 4. 
Sh.nb.size.dens <- glmmTMB(ShW ~ BulinusDens_sc + TotalSize_enclosure_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), family=nbinom2, data = data_subset_complete)
summary(Sh.nb.size.dens)
data_subset_complete$nb_model4_resids<-residuals(Sh.nb.size.dens)

# 5.
Sh.nb.snail.total <- glmmTMB(ShW ~ BulinusTotal_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), 
                             family=nbinom2, data = data_subset_complete)
summary(Sh.nb.snail.total)
data_subset_complete$nb_model5_resids<-residuals(Sh.nb.snail.total)

# 6. 
Sh.nb.size.dens.prev <- glmmTMB(ShW ~ BulinusDens_sc + TotalSize_enclosure_sc + ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), 
                                family=nbinom2, data = data_subset_complete)
summary(Sh.nb.size.dens.prev)
data_subset_complete$nb_model6_resids<-residuals(Sh.nb.size.dens.prev)

# 7.
Sh.nb.inf.snail.total <- glmmTMB(ShW ~ ShBulinus_Total_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), 
                                 family=nbinom2, data = data_subset_complete,na.action=na.omit)
summary(Sh.nb.inf.snail.total)
data_subset_complete$nb_model7_resids<-residuals(Sh.nb.inf.snail.total)

# 8. BulinusDens_sc + ShPrev_Bulinus_sc = ShBulinusDens_sc
Sh.nb.dens.prev <- glmmTMB(ShW ~ BulinusDens_sc + ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), 
                           family=nbinom2, data = data_subset_complete)
summary(Sh.nb.dens.prev)
data_subset_complete$nb_model8_resids<-residuals(Sh.nb.dens.prev)

# 9.
Sh.nb.inf.snail.dens <- glmmTMB(ShW ~ ShBulinusDens_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), 
                                family=nbinom2, data = data_subset_complete)
summary(Sh.nb.inf.snail.dens)
data_subset_complete$nb_model9_resids<-residuals(Sh.nb.inf.snail.dens)

# 10.
Sh.nb.size.perc.other <- glmmTMB(ShW ~ TotalSize_enclosure_sc + PercOther_sc + PercMud_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), 
                                 family=nbinom2, data = data_subset_complete)
summary(Sh.nb.size.perc.other)
data_subset_complete$nb_model10_resids<-residuals(Sh.nb.size.perc.other)

predictions<-predict(Sh.nb.size.perc.other,type="response")
data<-cbind(predictions,data_subset_complete$ShW)
data<-as.data.frame(data)

head(data)

b_ii<-ggplot()+
  geom_point(data=data,aes(x=data$V2,y=data$predictions))+
  theme_minimal()+
  geom_abline(slope=1,intercept=0)+
  xlab('')+
  ylab('')+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20,face="bold"))


# 11. 
Sh.nb.total.other <- glmmTMB(ShW ~ OtherVegTotal_sc + MudTotal_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), 
                             family=nbinom2, data = data_subset_complete)
summary(Sh.nb.total.other)
data_subset_complete$nb_model11_resids<-residuals(Sh.nb.total.other)

predictions<-predict(Sh.nb.total.other,type="response")
data<-cbind(predictions,data_subset_complete$ShW)
data<-as.data.frame(data)

head(data)

b_i<-ggplot()+
  geom_point(data=data,aes(x=data$V2,y=data$predictions))+
  theme_minimal()+
  geom_abline(slope=1,intercept=0)+
  xlab('')+
  ylab('')+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20,face="bold"))

# 12. 
Sh.nb.size.perc.vegmass <- glmmTMB(ShW ~ TotalSize_enclosure_sc + PercOther_sc + PercMud_sc + VegMassAvg_sampled_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), 
                                   family=nbinom2, data = data_subset_complete)
summary(Sh.nb.size.perc.vegmass)
data_subset_complete$nb_model12_resids<-residuals(Sh.nb.size.perc.vegmass)

predictions<-predict(Sh.nb.size.perc.vegmass,type="response")
data<-cbind(predictions,data_subset_complete$ShW)
data<-as.data.frame(data)

b_iii<-ggplot()+
  geom_point(data=data,aes(x=data$V2,y=data$predictions))+
  theme_minimal()+
  geom_abline(slope=1,intercept=0)+
  xlab('')+
  ylab('')+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20,face="bold"))

# 13. 
Sh.nb.veg.total.mass <- glmmTMB(ShW ~ VegMassTotal_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), 
                                family=nbinom2, 
                                data = data_subset_complete)
summary(Sh.nb.veg.total.mass)
data_subset_complete$nb_model13_resids<-residuals(Sh.nb.veg.total.mass)

predictions<-predict(Sh.nb.veg.total.mass,type="response")
data<-cbind(predictions,data_subset_complete$ShW)
data<-as.data.frame(data)

b_iv<-ggplot()+
  geom_point(data=data,aes(x=data$V2,y=data$predictions))+
  theme_minimal()+
  geom_abline(slope=1,intercept=0)+
  xlab('')+
  ylab('')+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20,face="bold"))

## full orthogonal model
# 14. TotalSize_enclosure_sc + PercOther_sc + PercMud_sc + VegMassAvg_sampled_sc + BulinusDens_sc + ShPrev_Bulinus_sc
Sh.nb.full.orthog <- glmmTMB(ShW ~ TotalSize_enclosure_sc + PercOther_sc + PercMud_sc + VegMassAvg_sampled_sc + BulinusDens_sc + ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), 
                             family=nbinom2,
                             data = data_subset_complete)
summary(Sh.nb.full.orthog)
data_subset_complete$nb_model14_resids<-residuals(Sh.nb.full.orthog)

predictions<-predict(Sh.nb.full.orthog,type="response")
data<-cbind(predictions,data_subset_complete$ShW)
data<-as.data.frame(data)

b_v<-ggplot()+
  geom_point(data=data,aes(x=data$V2,y=data$predictions))+
  theme_minimal()+
  geom_abline(slope=1,intercept=0)+
  xlab('')+
  ylab('')+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20,face="bold"))


# make a plot of observed versus predicted for all interpreted models (SI Appendix, Figure S7)

b<-ggdraw(plot=NULL,xlim=c(0,6),ylim=c(0,1))+
  draw_plot(b_i,x=0,y=0,width=1,height=1)+
  draw_plot(b_ii,x=1,y=0,width=1,height=1)+
  draw_plot(b_iii,x=2,y=0,width=1,height=1)+
  draw_plot(b_iv,x=3,y=0,width=1,height=1)+
  draw_plot(b_v,x=4,y=0,width=1,height=1)+
  draw_text("observed",x=3,y=0.01,size=18,hjust=0,vjust=0)+
  draw_text("predicted",x=0.09,y=0.45,size=18,hjust=0,vjust=0,angle=90)+
  draw_text("(1)",x=0.85,y=0.95,size=24)+
  draw_text("(2)",x=1.85,y=0.95,size=24)+
  draw_text("(3)",x=2.85,y=0.95,size=24)+
  draw_text("(4)",x=3.85,y=0.95,size=24)+
  draw_text("(5)",x=4.85,y=0.95,size=24)

final_prediction_plot<-ggdraw(plot=NULL,xlim=c(0,10),ylim=c(0,2.15))+
  draw_plot(a,x=0,y=1,width=10,height=1)+
  draw_plot(b,x=0,y=0,width=10,height=1)+
  draw_text("(a)",x=0.2,y=2.1,size=30)+
  draw_text("(b)",x=0.2,y=1.1,size=30)


### TEST FOR SPATIAL AUTOCORRELATION IN THESE MODELS
# Make a matrix of villages and residuals for each observation

residual_matrix_nb_merge<-cbind.data.frame(data_subset_complete$mergevillage,as.numeric(data_subset_complete$nb_model1_resids),
                                           as.numeric(data_subset_complete$nb_model2_resids),
                                           as.numeric(data_subset_complete$nb_model3_resids),as.numeric(data_subset_complete$nb_model4_resids),
                                           as.numeric(data_subset_complete$nb_model5_resids),as.numeric(data_subset_complete$nb_model6_resids),
                                           as.numeric(data_subset_complete$nb_model7_resids),as.numeric(data_subset_complete$nb_model8_resids),
                                           as.numeric(data_subset_complete$nb_model9_resids),as.numeric(data_subset_complete$nb_model10_resids),
                                           as.numeric(data_subset_complete$nb_model11_resids),as.numeric(data_subset_complete$nb_model12_resids),
                                           as.numeric(data_subset_complete$nb_model13_resids),as.numeric(data_subset_complete$nb_model14_resids))
residual_matrix_nb_merge<-as.data.frame(residual_matrix_nb_merge)
names(residual_matrix_nb_merge)<-c("village","model1","model2","model3","model4","model5","model6","model7","model8","model9",
                                   "model10","model11","model12","model13","model14")

aggregated_nb_merge<-aggregate(cbind(residual_matrix_nb_merge$model1,residual_matrix_nb_merge$model2,residual_matrix_nb_merge$model3,
                                     residual_matrix_nb_merge$model4,residual_matrix_nb_merge$model5,residual_matrix_nb_merge$model6,
                                     residual_matrix_nb_merge$model7,residual_matrix_nb_merge$model8,residual_matrix_nb_merge$model9,
                                     residual_matrix_nb_merge$model10,residual_matrix_nb_merge$model11,residual_matrix_nb_merge$model12,
                                     residual_matrix_nb_merge$model13,residual_matrix_nb_merge$model14),
                               by=list(residual_matrix_nb_merge$village),FUN=mean)
names(aggregated_nb_merge)<-c("village","model1","model2","model3","model4","model5","model6","model7","model8","model9",
                              "model10","model11","model12","model13","model14")

# based on recommendation here: https://stats.stackexchange.com/questions/40667/testing-for-spatial-autocorrelation-in-a-negative-binomial-regression-model

# need to make a listw using the nb2listw function
# start by importing all of the village-level polygons

village_polygons<-readOGR("data/sitesmerged-polygons.shp")
str(village_polygons)
village_polygons$Name

# with help from this resource, we calculated the centroid of each polygon - not the average of the coordinates in 
# each dimension, but weights the component triangles of the polygon by area:
# https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf

coords<-coordinates(village_polygons)
IDs<-c("Diokhor","Diokoul","Guidik","Lampsar","Malla","Malla Tack","MbakhanaMbarigot",
       "Mbane","Merina Gewel","Ndiawdoune","Ndiol Maure","Syer")
village_nb<-tri2nb(coords,row.names=IDs)
class(village_nb)
length(village_nb)

# now we need to use the polygon centroids to create a listw object that will tell the Moran's i calculation how to
#w eight the villages against one another. I want them to be weighted by distance.

dlist <- nbdists(village_nb, coords, longlat=T) 
dlist <- lapply(dlist, function(x) 1/x)

lw<-nb2listw(village_nb,dlist,zero.policy = T)
length(lw)
lw$weights

moran.mc(aggregated_nb_merge$model1,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb_merge$model2,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb_merge$model3,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb_merge$model4,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb_merge$model5,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb_merge$model6,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb_merge$model7,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb_merge$model8,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb_merge$model9,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb_merge$model10,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb_merge$model11,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb_merge$model12,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb_merge$model13,lw,na.action = na.omit,nsim=999)
moran.mc(aggregated_nb_merge$model14,lw,na.action = na.omit,nsim=999)



# Negative binomial model selection

## compare models by BIC
Sh.nb.mod.list <- c(Sh.nb.null.model, Sh.nb.density.only, Sh.nb.prev.only, Sh.nb.size.dens, Sh.nb.snail.total, Sh.nb.size.dens.prev, Sh.nb.inf.snail.total, Sh.nb.dens.prev, Sh.nb.inf.snail.dens, Sh.nb.size.perc.other, Sh.nb.total.other, Sh.nb.size.perc.vegmass, Sh.nb.veg.total.mass, Sh.nb.full.orthog)
Sh.nb.mod.names <- c('Sh.nb.null.model','Sh.nb.density.only','Sh.nb.prev.only','Sh.nb.size.dens', 'Sh.nb.snail.total', 'Sh.nb.size.dens.prev', 'Sh.nb.inf.snail.total', 'Sh.nb.dens.prev', 'Sh.nb.inf.snail.dens', 'Sh.nb.size.perc.other','Sh.nb.total.other', 'Sh.nb.size.perc.vegmass', 'Sh.nb.veg.total.mass', 'Sh.nb.full.orthog')

# BIC and MSE for all models is calculated in a separate script.
# For negative binomial models, "Analysis 4: Identifying snail- and habitat-related predictors of human urogenital schistosomiasis burden - Calculating BIC and MSE for negative binomial models"
# We used the BICs derived in those scripts to rank the models.

## models within 10 delta BIC of top model
Sh.nb.mod.subset <- c(Sh.nb.total.other, Sh.nb.size.perc.other, Sh.nb.size.perc.vegmass, Sh.nb.veg.total.mass, 
                      Sh.nb.full.orthog)

#get the marginal and conditional R2 for all models
tab_model(Sh.nb.null.model, Sh.nb.density.only, Sh.nb.prev.only, Sh.nb.size.dens, Sh.nb.snail.total, Sh.nb.size.dens.prev, Sh.nb.inf.snail.total, Sh.nb.dens.prev, Sh.nb.inf.snail.dens, Sh.nb.size.perc.other, Sh.nb.total.other, Sh.nb.size.perc.vegmass, Sh.nb.veg.total.mass, Sh.nb.full.orthog,
          show.r2 = TRUE, show.icc = FALSE, show.aic = FALSE, dv.labels = c('Sh.nb.null.model','Sh.nb.density.only','Sh.nb.prev.only','Sh.nb.size.dens',
                                                                            'Sh.nb.snail.total', 'Sh.nb.size.dens.prev', 'Sh.nb.inf.snail.total', 'Sh.nb.dens.prev',
                                                                            'Sh.nb.inf.snail.dens', 'Sh.nb.size.perc.other','Sh.nb.total.other', 'Sh.nb.size.perc.vegmass',
                                                                            'Sh.nb.veg.total.mass', 'Sh.nb.full.orthog'))



# top negative binomial model coefficient plot (main text, Figure 3B)

# top models = Sh.nb.total.other, Sh.nb.size.perc.other, Sh.nb.size.perc.vegmass, Sh.nb.veg.total.mass, Sh.nb.full.orthog
# get predictions for each model
# calculate CI for each estimate in model

M1.nb.ci <- as.data.frame(confint(Sh.nb.total.other)) #%>% # get confidence intervals from bblme package
M1.nb.preds <- tidy(Sh.nb.total.other)
M1.nb.preds <- bind_cols(M1.nb.preds,M1.nb.ci[1:9,]) %>% 
  dplyr::rename(conf.low = '2.5 %',
                conf.high = '97.5 %') %>% 
  mutate(Model="Model 1",
         estimate = exp(estimate),
         conf.low = exp(conf.low), 
         conf.high = exp(conf.high))

M2.nb.ci <- as.data.frame(confint(Sh.nb.size.perc.other)) #%>% # get confidence intervals from bblme package
M2.nb.preds <- tidy(Sh.nb.size.perc.other)
M2.nb.preds <- bind_cols(M2.nb.preds,M2.nb.ci[1:10,]) %>% 
  dplyr::rename(conf.low = '2.5 %',
                conf.high = '97.5 %') %>% 
  mutate(Model="Model 2",
         estimate = exp(estimate),
         conf.low = exp(conf.low), 
         conf.high = exp(conf.high))

M3.nb.ci <- as.data.frame(confint(Sh.nb.size.perc.vegmass)) #%>% # get confidence intervals from bblme package
M3.nb.preds <- tidy(Sh.nb.size.perc.vegmass)
M3.nb.preds <- bind_cols(M3.nb.preds,M3.nb.ci[1:11,]) %>% 
  dplyr::rename(conf.low = '2.5 %',
                conf.high = '97.5 %') %>% 
  mutate(Model="Model 3",
         estimate = exp(estimate),
         conf.low = exp(conf.low), 
         conf.high = exp(conf.high))

M4.nb.ci <- as.data.frame(confint(Sh.nb.veg.total.mass)) #%>% # get confidence intervals from bblme package
M4.nb.preds <- tidy(Sh.nb.veg.total.mass)
M4.nb.preds <- bind_cols(M4.nb.preds,M4.nb.ci[1:8,]) %>% 
  dplyr::rename(conf.low = '2.5 %',
                conf.high = '97.5 %') %>% 
  mutate(Model="Model 4",
         estimate = exp(estimate),
         conf.low = exp(conf.low), 
         conf.high = exp(conf.high))

M5.nb.ci <- as.data.frame(confint(Sh.nb.full.orthog)) #%>% # get confidence intervals from bblme package
M5.nb.preds <- tidy(Sh.nb.full.orthog)
M5.nb.preds <- bind_cols(M5.nb.preds,M5.nb.ci[1:13,]) %>% 
  dplyr::rename(conf.low = '2.5 %',
                conf.high = '97.5 %') %>% 
  mutate(Model="Model 5",
         estimate = exp(estimate),
         conf.low = exp(conf.low), 
         conf.high = exp(conf.high))

# bind to one df
ShW.odds = bind_rows(M1.nb.preds,M2.nb.preds,M3.nb.preds,M4.nb.preds,M5.nb.preds) %>% 
  filter(term != '(Intercept)')
ShW.odds$term

ShW.odds <- ShW.odds %>% filter(effect!="ran_pars")
dodger = position_dodge(width = 0.9)

# Elements like pointrange and position_dodge only work when the outcome
# is mapped to y, need to go through with OR set as y then flip at the end

ggplot(ShW.odds, aes(y = estimate, x = term, color=term, shape = Model)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),position = dodger,size = 0.4) +
  scale_colour_manual(values=c("black","black","black","black","black","black","black","black","black","black","black","black","black","black")) +
  scale_x_discrete(limits=c("ShPrev_Bulinus_sc","BulinusDens_sc","VegMassTotal_sc","PercMud_sc",
                            "PercOther_sc","MudTotal_sc","OtherVegTotal_sc","TotalSize_enclosure_sc","Class",
                            "Pop_sc","sexM","LakeYN1"),labels=c("snail prevalence","snail density","mass of non-emergent vegetation",
                                                                "percent cover of mud",
                                                                "percent cover of non-emergent vegetation",
                                                                "area of mud",
                                                                "area of non-emergent vegetation",
                                                                "site area",
                                                                "school grade", 
                                                                "village population","sex: male (vs. female)", 
                                                                "location: lake (vs. river)")) +
  geom_hline(yintercept = 1.0, linetype = "dotted", size = 0.5) +
  scale_y_log10(breaks = c(0.25, 0.5, 1.0, 2.0, 5.0, 10, 20),
                minor_breaks = NULL) +
  labs(y = "odds ratio", x = "predictor") +
  theme(legend.position = "none") +
  coord_flip(ylim = c(0.25, 30)) +
  theme_classic()

# show model results as dataframe 
tab_model(Sh.nb.total.other, Sh.nb.size.perc.other, Sh.nb.size.perc.vegmass, Sh.nb.veg.total.mass, Sh.nb.full.orthog,
          transform = "exp",
          pred.labels = c(),
          dv.labels = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5"),
          CSS = list(css.centeralign = 'text-align: left;'), 
          collapse.ci = TRUE,
          show.r2 = FALSE, show.icc = FALSE,
          show.aic = FALSE)
