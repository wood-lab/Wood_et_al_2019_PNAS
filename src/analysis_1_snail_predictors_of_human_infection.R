---
  title: "Analysis 1: Identifying snail-related predictors of human urogenital schistosomiasis burden"
author: "Isabel Jones and Chelsea Wood"
date: "updated 9 July 2019"
---
  
# Introduction

  # We want to know which snail variables best predict human schistosomiasis infection rates, here measured as reinfection to 
  # children over the course of one year following praziquantel treatment. 
  
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
library(visreg)
library(sjPlot)
library(ggeffects)
library(plyr)
library(broom) # for tidy
library(broom.mixed)
library(PerformanceAnalytics) # for correlation plots
library(AICcmodavg)
library(MuMIn) # for averaging model predictions
library(gridExtra)
library(bbmle) # AICtab for negative binomial
library(DHARMa) # for visualizing glmmTMB resids
library(Hmisc) # for CI in ggplot
library(fitdistrplus) # for checking nb fits and calculating mu


####### Download data - note that the individual-level data on which this analysis was conducted are not provided, to protect patient privacy.

villagedata <- read_csv("INDIVIDUAL-LEVEL DATASET NOT PROVIDED TO PROTECT PATIENT PRIVACY")
data_subset <- villagedata %>% 
  filter((year==2017 & net!=1) | (year==2018 & (prawn|veg_removal)!=1))
unique(data_subset$Village)
data_subset$LakeYN=as.factor(data_subset$LakeYN)
data_subset_complete <- data_subset[complete.cases(data_subset), ]
data_subset_complete %>% summarise(unique(Village))
data_subset_complete[data_subset_complete$year==2017,] %>% summarise(unique(Village))
data_subset_complete[data_subset_complete$year==2018,] %>% summarise(unique(Village))
# check which villages are included in each year
data_subset %>% 
  filter(year==2017) %>% 
  distinct(Village)
data_subset %>% summarise(n_distinct(Village))
# check number of kids in study
data_subset %>% summarise(n_distinct(ID))


# Upon testing for spatial autocorrelation, we discovered that there was mild but significant positive spatial autocorrelation in the residuals
# for the logistic models (i.e., models of re-infection probability) and mild but non-significant positive spatial autocorrelation in the 
# residuals for the negative binomial models (i.e., models of egg count). To address this issue, we grouped the nearest pair of villages, 
# Mbakhana and Mbarigot, into a single village. This eliminated the spatial autocorrelation issue. For details, see SI Appendix, Text S4.

data_subset_complete$mergevillage<-gsub("Mbakhana|Mbarigot","MbakhanaMbarigot",data_subset_complete$Village)


# Now that we've got Mbakhana and Mbarigot merged into a single village, we can proceed with analysis

# Predictors 

# demographic variables and random effects:
#  Class + sex + Pop_sc + LakeYN + (1|Village/ID)

# nested predictors --- orthogonal expanded and composite, to be tested in separate models
# 1/2. TotalSize_enclosure_sc + BulinusDens_sc = BulinusTotal_sc
# 3/4. TotalSize_enclosure_sc + BulinusDens_sc + ShPrev_Bulinus_sc = ShBulinus_Total_sc
# 5/6. BulinusDens_sc + ShPrev_Bulinus_sc = ShBulinusDens_sc
# 7/8. TotalSize_enclosure_sc + PercOther_sc + PercMud_sc = OtherVegTotal_sc + MudTotal_sc
# 9/10. TotalSize_enclosure_sc + PercOther_sc + PercMud_sc + VegMassAvgSampled_sc = VegMassTotal_sc

# full orthogonal model
# 11. TotalSize_enclosure_sc + PercOther_sc + PercMud_sc + VegMassAvgSampled_sc + BulinusDens_sc + ShPrev_Bulinus_sc

# Run the Sh logistic models (no interactions)

########################################################
###########  S haematobium logistic models   ########### 
# every model: Class + sex + Pop_sc + LakeYN + (1|Village/ID)
# 1. 
size.dens <- glmer(Sh ~ BulinusDens_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
                   family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(size.dens)
r.squaredGLMM(size.dens)

# 2.
snail.total <- glmer(Sh ~ BulinusTotal_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
                     family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(snail.total)
r.squaredGLMM(snail.total)

# 3.
size.dens.prev <- glmer(Sh ~ BulinusDens_sc + ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
                        family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(size.dens.prev)
r.squaredGLMM(size.dens.prev)

# 4.
inf.snail.total <- glmer(Sh ~ ShBulinus_Total_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
                         family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(inf.snail.total)
r.squaredGLMM(inf.snail.total)

# 5
prev <- glmer(Sh ~ ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
              family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(prev)
r.squaredGLMM(prev)

# 6. 
inf.snail.dens <- glmer(Sh ~ ShBulinusDens_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
                        family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(inf.snail.dens)
r.squaredGLMM(inf.snail.dens)

#7. NULL MODEL
null.model <- glmer(Sh ~ Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, 
                    family = 'binomial', control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(null.model)
r.squaredGLMM(null.model)


# SELECT THE TOP MODELS BY BIC


# Top model coefficient plot for Sh logistic model

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

Sh.odds$term


# set order of terms
# Sh.order <- rev(c("LakeYN1","TotalSize_enclosure_sc","PercOther_sc","OtherVegTotal_sc",
#                  "VegMassAvg_sampled_sc","StemsTotal_sc","BulinusDens_sc","ShPrev_Bulinus_sc","SexM","Class"))
# Sh.odds$term <- factor(Sh.odds$term, levels = Sh.order)
dodger = position_dodge(width = 0.9)
# Elements like pointrange and position_dodge only work when the outcome
#   is mapped to y, need to go through with OR set as y then flip at the end
ggplot(Sh.odds[Sh.odds$effect=='fixed',], aes(y = estimate, x = term, color=term, shape = Model)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),position = dodger,size = 0.4) +
  scale_colour_manual(values=c("black","black","black","black","black","black","black","black","black")) +
  scale_x_discrete(limits=c("ShBulinusDens_sc","ShPrev_Bulinus_sc","BulinusDens_sc","ShBulinus_Total_sc","BulinusTotal_sc","Class",
                            "Pop_sc","sexM","LakeYN1"),labels=c("infected snail density","snail prevalence","snail density","infected snail abundance",
                                                                "snail abundance","school grade", "village population","sex: male (vs. female)", 
                                                                "location: lake (vs. river)")) +
  geom_hline(yintercept = 1.0, linetype = "dotted", size = 0.5) +
  scale_y_log10(breaks = c(0.25, 0.5, 1.0, 2.0, 5.0, 10, 20),
                minor_breaks = NULL) +
  labs(y = "odds ratio", x = "predictor") +
  theme(legend.position = "none") +
  coord_flip(ylim = c(0.25, 40)) +
  theme_classic() 
ggsave("Eco Predictors/Figures for Chelsea paper/Sh Odds plot.png")

# show model results as dataframe 
tab_model(null.model, snail.total, inf.snail.total, size.dens, prev, inf.snail.dens, 
          transform = "exp",
          pred.labels = c(),
          dv.labels = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6"),
          CSS = list(css.centeralign = 'text-align: left;'), 
          collapse.ci = TRUE,
          show.r2 = TRUE, show.icc = FALSE,
          show.aic = TRUE)


---------------------------------------------------------------------

  
# Sh negative binomial models 

#################       S haematobium negative binomial models     #################       
# using glmmTMB with negative binomial with variance structure of (variance = µ(1 + µ/k): Hardin and Hilbe (2007)). 
# every model: PercMud_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID)

# 1. 
Sh.nb.size.dens <- glmmTMB(ShW ~ BulinusDens_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), family=nbinom2, data = data_subset_complete)
summary(Sh.nb.size.dens)
r2(Sh.nb.size.dens)

# 2.
Sh.nb.snail.total <- glmmTMB(ShW ~ BulinusTotal_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), 
                             family=nbinom2, data = data_subset_complete)
summary(Sh.nb.snail.total)
r2(Sh.nb.snail.total)

# 3. 
Sh.nb.dens.prev <- glmmTMB(ShW ~ BulinusDens_sc + ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), 
                           family=nbinom2, data = data_subset_complete)
summary(Sh.nb.dens.prev)
r2(Sh.nb.size.dens.prev)

# 4.
Sh.nb.inf.snail.total <- glmmTMB(ShW ~ ShBulinus_Total_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), 
                                 family=nbinom2, data = data_subset_complete,na.action=na.omit)
summary(Sh.nb.inf.snail.total)
r2(Sh.nb.inf.snail.total)

# 5. 
Sh.nb.prev <- glmmTMB(ShW ~ ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), 
                      family=nbinom2, data = data_subset_complete)
summary(Sh.nb.prev)

# 6.
Sh.nb.inf.snail.dens <- glmmTMB(ShW ~ ShBulinusDens_sc + Class + sex + Pop_sc + LakeYN + (1|Village/ID), 
                                family=nbinom2, data = data_subset_complete)
summary(Sh.nb.inf.snail.dens)

# 7.
null.model.nb <- glmmTMB(ShW ~ Class + sex + Pop_sc + LakeYN + (1|Village/ID), 
                         family=nbinom2, data = data_subset_complete)
summary(null.model.nb)


# USE BIC TO FIND TOP MODELS


# Top ShW NB model coefficient plot

```{r ShW plot top model odds ratios, echo=FALSE}
AICtab(Sh.nb.size.dens, Sh.nb.snail.total, Sh.nb.size.dens.prev, Sh.nb.inf.snail.total, Sh.nb.dens.prev, Sh.nb.inf.snail.dens, Sh.nb.size.perc.other, Sh.nb.total.other, Sh.nb.size.perc.vegmass, Sh.nb.veg.total.mass, Sh.nb.full.orthog, mnames = Sh.nb.mod.names)
# get predictions for each model
# calculate CI for each estimate in model
M1.nb.ci <- as.data.frame(confint(Sh.nb.snail.total)) #%>% # get confidence intervals from bblme package
M1.nb.preds <- tidy(Sh.nb.snail.total)
M1.nb.preds <- bind_cols(M1.nb.preds,M1.nb.ci[1:8,]) %>% 
  dplyr::rename(conf.low = '2.5 %',
                conf.high = '97.5 %') %>% 
  mutate(Model="Model 1",
         estimate = exp(estimate),
         conf.low = exp(conf.low), 
         conf.high = exp(conf.high))

M2.nb.ci <- as.data.frame(confint(Sh.nb.inf.snail.total)) #%>% # get confidence intervals from bblme package
M2.nb.preds <- tidy(Sh.nb.inf.snail.total)
M2.nb.preds <- bind_cols(M2.nb.preds,M2.nb.ci[1:8,]) %>% 
  dplyr::rename(conf.low = '2.5 %',
                conf.high = '97.5 %') %>% 
  mutate(Model="Model 2",
         estimate = exp(estimate),
         conf.low = exp(conf.low), 
         conf.high = exp(conf.high))


# bind to one df
ShW.odds = bind_rows(M1.nb.preds,M2.nb.preds) %>% 
  filter(term != '(Intercept)')
# TO-DO:
# set order of terms
# Sh.order <- rev(c("LakeYN1","TotalSize_enclosure_sc","PercOther_sc","OtherVegTotal_sc",
#                  "VegMassAvg_sampled_sc","PercMud_sc","BulinusDens_sc","ShPrev_Bulinus_sc","SexM","Class"))
# Sh.odds$term <- factor(Sh.odds$term, levels = Sh.order)
min(ShW.odds$estimate)
max(ShW.odds$estimate)
min(ShW.odds$conf.low)
max(ShW.odds$conf.high)
ShW.odds <- ShW.odds %>% filter(effect!="ran_pars")
dodger = position_dodge(width = 0.9)

ShW.odds$term
# Elements like pointrange and position_dodge only work when the outcome
#   is mapped to y, need to go through with OR set as y then flip at the end


ggplot(filter(ShW.odds,effect=="fixed"), aes(y = estimate, x = term, color=term, shape = Model)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),position = dodger,size = 0.4) +
  scale_colour_manual(values=c("black","black","black","black","black","black","black","black","black")) +
  scale_x_discrete(limits=c("ShBulinus_Total_sc","BulinusTotal_sc",
                            "Class","Pop_sc","sexM","LakeYN1"),labels=c("infected snail abundance","snail abundance",
                                                                        "school grade","village population","sex: male (vs. female)", 
                                                                        "location: lake (vs. river)")) +
  geom_hline(yintercept = 1.0, linetype = "dotted", size = 0.5) +
  scale_y_log10(breaks = c(0.25, 0.5, 1.0, 2.0, 5.0, 10, 20, 40),
                minor_breaks = NULL) +
  labs(y = "incidence rate ratio", x = "predictor") +
  theme(legend.position = "none") +
  coord_flip(ylim = c(0.25, 70)) +
  theme_classic() 



# show model results as dataframe 
tab_model(Sh.nb.snail.total, Sh.nb.inf.snail.total,
          transform = "exp",
          pred.labels = c(),
          dv.labels = c("Model 1", "Model 2"),
          CSS = list(css.centeralign = 'text-align: left;'), 
          collapse.ci = TRUE,
          show.r2 = TRUE, show.icc = FALSE,
          show.aic = TRUE)
```