---
title: "Analysis 4: Identifying snail- and habitat-related predictors of human urogenital schistosomiasis burden - Calculating BIC and MSE for logistic/binomial models"
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
villagedata <- read_csv("data/ecopredictors_data.csv")
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

# Model control
mod_ctl <- glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))

# Fit
# 1.
mod_list[[1]] <- glmer(Sh ~ Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, family = 'binomial', control=mod_ctl)

# 2.
mod_list[[2]] <- glmer(Sh ~ BulinusDens_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, family = 'binomial', control=mod_ctl)

# 3.
mod_list[[3]] <- glmer(Sh ~ ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, family = 'binomial', control=mod_ctl)

# 4. 
mod_list[[4]] <- glmer(Sh ~ BulinusDens_sc + TotalSize_enclosure_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, family = 'binomial', control=mod_ctl)

# 5.
mod_list[[5]] <- glmer(Sh ~ BulinusTotal_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, family = 'binomial', control=mod_ctl)

# 6.
mod_list[[6]] <- glmer(Sh ~ BulinusDens_sc + TotalSize_enclosure_sc + ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, family = 'binomial', control=mod_ctl)

# 7.
mod_list[[7]] <- glmer(Sh ~ ShBulinus_Total_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete,  family = 'binomial', control=mod_ctl)

# 8.
mod_list[[8]] <- glmer(Sh ~ BulinusDens_sc + ShPrev_Bulinus_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, family = 'binomial', control=mod_ctl)

# 9. 
mod_list[[9]] <- glmer(Sh ~ ShBulinusDens_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, family = 'binomial', control=mod_ctl)

# 10
mod_list[[10]] <- glmer(Sh ~ TotalSize_enclosure_sc + PercOther_sc + PercMud_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, family = 'binomial', control=mod_ctl)

# 11. 
mod_list[[11]] <- glmer(Sh ~ OtherVegTotal_sc + MudTotal_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, family = 'binomial', control=mod_ctl)

# 12. 
mod_list[[12]] <- glmer(Sh ~ TotalSize_enclosure_sc + VegMassAvg_sampled_sc + PercOther_sc + PercMud_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, family = 'binomial', control=mod_ctl)

# 13.
mod_list[[13]] <- glmer(Sh ~ VegMassTotal_sc + Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete,  family = 'binomial', control=mod_ctl)

# 14. 
mod_list[[14]] <- glmer(Sh ~ TotalSize_enclosure_sc + VegMassAvg_sampled_sc + BulinusDens_sc + ShPrev_Bulinus_sc + PercOther_sc + PercMud_sc +  Class + sex + Pop_sc + LakeYN + (1|mergevillage/ID), data = data_subset_complete, family = 'binomial', control=mod_ctl)

###############################################
# Information criterion

# AIC
aic_vec <- sapply(mod_list, AIC)
daic_vec <- aic_vec - min(aic_vec)
aic_weights <- exp(-0.5 * daic_vec)
aic_weights <- aic_weights / sum(aic_weights)

# AICc
aicc_vec <- sapply(mod_list, AICc)
daicc_vec <- aicc_vec - min(aicc_vec)
aicc_weights <- exp(-0.5 * daicc_vec)
aicc_weights <- aicc_weights / sum(aicc_weights)

# BIC
bic_vec <- sapply(mod_list, BIC)
dbic_vec <- bic_vec - min(bic_vec)
bic_weights <- exp(-0.5 * dbic_vec)
bic_weights <- bic_weights / sum(bic_weights)

###############################################
# MSE Run

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
    library(AICcmodavg)
    library(lme4)
    mod_ctl <- glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))
    
    # Subset data
    train_dat <- data_subset_complete[-i,]
    test_dat <- data_subset_complete[i,]
    
    # List to store models
    mod_tmp <- list()
    
    # Refit models
    print(paste0("Fitting model subset ",i))
    
    # Refit
    mod_tmp <- sapply(mod_list, function(x) update(x, data = train_dat))
    
    # AIC of new models
    aic_vec_tmp <- sapply(mod_tmp, AIC)
    daic_vec_tmp <- aic_vec_tmp - min(aic_vec_tmp)
    aic_weights_tmp <- exp(-0.5 * daic_vec_tmp)
    aic_weights_tmp <- aic_weights_tmp / sum(aic_weights_tmp)
    
    # AICc of new models
    aicc_vec_tmp <- sapply(mod_tmp, AICc )
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
    fix_effects <- lapply(mod_tmp, function(x) fixef(x))
    
    # Multiply fixed effects and model matrix
    mat_mult <- Map('*',mod_mat_tmp,fix_effects)
    
    # Sum and put on natural scale
    pred_mat <- sapply(mat_mult, function(x) 1/(1+exp(-sum(x))))
    
    # Model average predictions
    pred_mat <- c(pred_mat, sum(pred_mat[1:length(mod_list)] * aic_weights_tmp))  
    pred_mat <- c(pred_mat, sum(pred_mat[1:length(mod_list)] * aicc_weights_tmp))
    c(pred_mat, sum(pred_mat[1:length(mod_list)] * bic_weights_tmp))
  }
})


# When you're done, clean up the cluster
stopImplicitCluster()
pred_prob_mat <- mse_pred

# Get probs
expected_outcome_mat <- ifelse( mse_pred > 0.5, 1, 0)

# Get prediction accuracty
acc_mat <- mse_pred == data_subset_complete$Sh


# Turn to matrices and name columns/rows
pred_prob_mat <- as.data.frame(pred_prob_mat); colnames(pred_prob_mat) <- paste0("Model",1:ncol(pred_prob_mat)); rownames(pred_prob_mat) <- paste0("Observation_",1:nrow(pred_prob_mat))
expected_outcome_mat <- as.data.frame(expected_outcome_mat); colnames(expected_outcome_mat) <- paste0("Model",1:ncol(expected_outcome_mat)); rownames(expected_outcome_mat) <- paste0("Observation_",1:nrow(expected_outcome_mat))
acc_mat <- as.data.frame(acc_mat); colnames(acc_mat) <- paste0("Model",1:ncol(acc_mat)); rownames(acc_mat) <- paste0("Observation_",1:nrow(acc_mat))

# Check things out
head(pred_prob_mat) # Predicted probability of having schito varies across models
head(expected_outcome_mat) # Expected outcome is pretty much the same
head(acc_mat) # So is accuracy

# Calculate mean prediction accuracy
acc_vec <- colMeans(acc_mat, na.rm = TRUE); acc_vec # These are the results

# TRUE MSE
mse_vec <- colMeans((pred_prob_mat - data_subset_complete$Sh)^2)

# AUC
auc_pred <- sapply(pred_prob_mat, function(x) prediction(x, data_subset_complete$Sh))
auc_perf <- sapply(auc_pred, function(x) performance(x, measure = "tpr", x.measure = "fpr"))
auc_vals <- sapply(auc_pred, function(x) performance(x, measure = "auc")@y.values[[1]])

# Combine performance metrics
results <- data.frame(Model = 1:17, AIC =  c(aic_vec, NA, NA, NA), dAIC = c(daic_vec, NA, NA, NA), AICw = c(aic_weights, NA, NA, NA), AICc = c(aicc_vec, NA, NA, NA), dAICc = c(daicc_vec, NA, NA, NA), AICcw = c(aicc_weights, NA, NA, NA), BIC = c(bic_vec, NA, NA, NA), dBIC = c(dbic_vec, NA, NA, NA), BICw = c(bic_weights, NA, NA, NA), AUC = auc_vals, ACC = acc_vec, MSE = mse_vec, Formula = c(as.character(sapply(mod_list, formula)), "avgAIC", "avgAICc", "avgBIC") )

write.csv(results, file = "binom_mse_results.csv")

pred_prob_mat$TrueValue <- data_subset_complete$Sh
write.csv(pred_prob_mat, file = "pred_prob_mat.csv")