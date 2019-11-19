---
title: "Analysis 4: Identifying snail- and habitat-related predictors of human urogenital schistosomiasis burden - Calculating BIC and MSE for negative binomial models"
author: "Grant Adams"
date: "updated 9 July 2019"
---
  
###############################################
# Load all the packages
library(plyr)
library(dplyr)
library(lme4)
library(tidyverse)
library(AICcmodavg)
library(ROCR)
library(glmmTMB)
library(foreach)
library(doParallel)

###############################################
# Set up parallel processing
numCores <- detectCores() - 1

###############################################
# Load all the data according to IZZY
villagedata <- read_csv("data/INDIVIDUAL-LEVEL DATA WITHHELD TO PROTECT PATIENT PRIVACY")
data_subset <- villagedata %>% 
  filter((year==2017 & net!=1) | (year==2018 & (prawn|veg_removal)!=1))
unique(data_subset$Village)
data_subset$LakeYN=as.factor(data_subset$LakeYN)
data_subset_complete <- data_subset[complete.cases(data_subset), ]

# Merge two sites: Mbakhana and Mbarigot. This helps to alleviate spatial autocorrelation.

data_subset_complete$Village
data_subset_complete$mergevillage<-gsub("Mbakhana|Mbarigot","MbakhanaMbarigot",data_subset_complete$Village)

#check that this worked properly
write.csv(data_subset_complete$mergevillage)


###############################################
# Fit base models
mod_list <- list()

# Fit
# 1.
mod_list[[1]] <- glmmTMB(ShW ~ Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), family=nbinom2, data = data_subset_complete)

# 2.
mod_list[[2]] <- glmmTMB(ShW ~ BulinusDens_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), family=nbinom2, data = data_subset_complete)

# 3.
mod_list[[3]] <- glmmTMB(ShW ~ ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), family=nbinom2, data = data_subset_complete)

# 4. 
mod_list[[4]] <- glmmTMB(ShW ~ BulinusDens_sc + TotalSize_enclosure_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), family=nbinom2, data = data_subset_complete)

# 5.
mod_list[[5]] <- glmmTMB(ShW ~ BulinusTotal_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), family=nbinom2, data = data_subset_complete)

# 6.
mod_list[[6]] <- glmmTMB(ShW ~ BulinusDens_sc + TotalSize_enclosure_sc + ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), family=nbinom2, data = data_subset_complete)

# 7.
mod_list[[7]] <- glmmTMB(ShW ~ ShBulinus_Total_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), family=nbinom2, data = data_subset_complete,na.action=na.omit)

# 8.
mod_list[[8]] <- glmmTMB(ShW ~ BulinusDens_sc + ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), family=nbinom2, data = data_subset_complete)

# 9. 
mod_list[[9]] <- glmmTMB(ShW ~ ShBulinusDens_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), family=nbinom2, data = data_subset_complete)

# 10
mod_list[[10]] <- glmmTMB(ShW ~ TotalSize_enclosure_sc + PercOther_sc + PercMud_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), family=nbinom2, data = data_subset_complete)

# 11. 
mod_list[[11]] <- glmmTMB(ShW ~ OtherVegTotal_sc + MudTotal_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), family=nbinom2, data = data_subset_complete)

# 12. 
mod_list[[12]] <- glmmTMB(ShW ~ TotalSize_enclosure_sc + PercOther_sc + PercMud_sc + VegMassAvg_sampled_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), family=nbinom2, data = data_subset_complete)

# 13.
mod_list[[13]] <- glmmTMB(ShW ~ VegMassTotal_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), family=nbinom2, data = data_subset_complete)

# 14. 
mod_list[[14]] <- glmmTMB(ShW ~ TotalSize_enclosure_sc + PercOther_sc + PercMud_sc + VegMassAvg_sampled_sc + BulinusDens_sc + ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), family=nbinom2, data = data_subset_complete)

###############################################
# Information criterion

# AIC
aic_vec <- sapply(mod_list, AIC)
daic_vec <- aic_vec - min(aic_vec)
aic_weights <- exp(-0.5 * daic_vec)
aic_weights <- aic_weights / sum(aic_weights)

# AICc
aicc_vec <- sapply(mod_list, AICcmodavg::AICc )
daicc_vec <- aicc_vec - min(aicc_vec)
aicc_weights <- exp(-0.5 * daicc_vec)
aicc_weights <- aicc_weights / sum(aicc_weights)

# BIC
bic_vec <- sapply(mod_list, BIC)
dbic_vec <- bic_vec - min(bic_vec)
bic_weights <- exp(-0.5 * dbic_vec)
bic_weights <- bic_weights / sum(bic_weights)

###############################################
# MSE Run set up

nmods <-  17 #  Models
nobs <- nrow(data_subset_complete) # number of obs 
# If running a subset to check for issues set nobs to 10 or something
# nobs <- 10

# Create MSE objects
mod_mat <- lapply(mod_list, function(x) model.matrix(x))


###############################################
# MSE Run parallel
registerDoParallel(numCores)  # use multicore, set to the number of our cores


stime <- system.time({
  mse_pred <- foreach(i = 1:nobs, .combine=rbind) %dopar% {
    library(glmmTMB)
    library(AICcmodavg)
    
    train_dat <- data_subset_complete[-i,] # Subset data
    test_dat <- data_subset_complete[i,]
    
    # List to store models
    mod_tmp <- list()
    
    # Refit models
    print(paste0("Fitting model subset ",i))
    mod_tmp <- lapply(mod_list, function(x) update(x, data = train_dat))
    
    # AIC of new models
    aic_vec_tmp <- sapply(mod_tmp, AIC)
    daic_vec_tmp <- aic_vec_tmp - min(aic_vec_tmp)
    aic_weights_tmp <- exp(-0.5 * daic_vec_tmp)
    aic_weights_tmp <- aic_weights_tmp / sum(aic_weights_tmp)
    
    # AICc of new models
    aicc_vec_tmp <- sapply(mod_tmp, AICcmodavg::AICc )
    daicc_vec_tmp <- aicc_vec_tmp - min(aicc_vec_tmp)
    aicc_weights_tmp <- exp(-0.5 * daicc_vec_tmp)
    aicc_weights_tmp <- aicc_weights_tmp / sum(aicc_weights_tmp)
    
    # BIC of new models
    bic_vec_tmp <- sapply(mod_tmp, BIC )
    dbic_vec_tmp <- bic_vec_tmp - min(bic_vec_tmp)
    bic_weights_tmp <- exp(-0.5 * dbic_vec_tmp)
    bic_weights_tmp <- bic_weights_tmp / sum(bic_weights_tmp)
    
    # Extract the model matrix from full model, but only the test row
    mod_mat_tmp <- lapply(mod_mat, function(x) x[i,])
    
    # Extract fixed effects
    fix_effects <- lapply(mod_tmp, function(x) fixef(x)$cond)
    
    # Multiply fixed effects and model matrix
    mat_mult <- Map('*',mod_mat_tmp,fix_effects)
    
    # Sum and put on natural scale
    pred_mat <- sapply(mat_mult, function(x) exp(sum(x)))
    
    # Model average predictions
    pred_mat <- c(pred_mat, sum(pred_mat[1:length(mod_list)] * aic_weights_tmp))  
    pred_mat <- c(pred_mat, sum(pred_mat[1:length(mod_list)] * aicc_weights_tmp))
    c(pred_mat, sum(pred_mat[1:length(mod_list)] * bic_weights_tmp))
    
  }
})


# When you're done, clean up the cluster
stopImplicitCluster()
pred_mat <- mse_pred


# Turn to matrices and name columns/rows
pred_mat <- as.data.frame(pred_mat); colnames(pred_mat) <- paste0("Model",1:ncol(pred_mat)); rownames(pred_mat) <- paste0("Observation_",1:nrow(pred_mat))

# Check things out
head(pred_mat)

# TRUE MSE
mse_vec <- colMeans((pred_mat - data_subset_complete$ShW)^2)
lmse_vec <- colMeans((log(pred_mat) - log(data_subset_complete$ShW))^2)

# Combine performance metrics
results <- data.frame(Model = 1:17, AIC =  c(aic_vec, NA, NA, NA), dAIC = c(daic_vec, NA, NA, NA), AICw = c(aic_weights, NA, NA, NA), AICc = c(aicc_vec, NA, NA, NA), dAICc = c(daicc_vec, NA, NA, NA), AICcw = c(aicc_weights, NA, NA, NA), BIC = c(bic_vec, NA, NA, NA), dBIC = c(dbic_vec, NA, NA, NA), BICw = c(bic_weights, NA, NA, NA), MSE = mse_vec, MSLE = lmse_vec, RMSE = sqrt(mse_vec), Formula = c(as.character(sapply(mod_list, formula)), "avgAIC", "avgAICc", "avgBIC") )

write.csv(results, file = "nbinom_mse_results.csv")

# Save values
pred_mat$TrueValue <- data_subset_complete$ShW
write.csv(pred_mat, file = "nbinom_mse_pred_mat.csv")
