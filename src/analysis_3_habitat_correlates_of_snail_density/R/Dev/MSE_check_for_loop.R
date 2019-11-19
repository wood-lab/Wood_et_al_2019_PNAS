model = 2; incl_disease = 0; bulinus = T; BIC_sel = TRUE; fill_cont_mean = F; fill_cat_mode = F; remove_outlier = FALSE; method = "nlminb"; silent = TRUE; savedir = "Report/June_Check"


# Load functions
source("R/Support/1_read_dat_v3.R")
source("R/Support/2_fit_mesh_v2.R")
source("R/Support/3_build_data_v3.R")
source("R/Support/4_build_params_v3.R")
source("R/Support/5_build_map_v2.R")
source("R/Support/bad_params.R")
source("R/run_variogram.R")
source("R/Support/aic_selection_v2.R")
source("R/plot_coef_v2.R")

# Creat the directory and file names
if(!is.null(savedir)){
  suppressWarnings(dir.create(file.path(savedir)))
}

file_name_base <- paste(c("biom", "bulinus")[as.numeric(bulinus)+1],"model",  model, c("AIC", "BIC")[as.numeric(BIC_sel)+1], c("no_mean_cont_fill", "mean_cont_fill")[as.numeric(fill_cont_mean)+1], c("no_mode_cat_fill", "mode_cat_fill")[as.numeric(fill_cat_mode)+1], "remove_outlier", remove_outlier, sep = "-")


# Step 0 -- Load data, create mesh, and assign to list for TMB
data <- read_dat(exclude = T,  fill_cont_mean = fill_cont_mean, na_keep = T, quad = F, fill_cat_mode = fill_cat_mode, scale_cont = T, bulinus = bulinus, remove_outlier = remove_outlier) # read data in


# data <- fit_mesh(data, cluster = T, n_knots = 15) # Fit mesh
data_list <- build_data(data, model, debug = 0, bulinus = bulinus, parasite = 0) # Create for TMB


# Step 1 -- Set initial values for parameters
params <- build_params(model, data_list, incl_disease)
random = c("epsilon_mat", "omega_mat", "gamma_q", "x_q_cont_missing")


# Step 3 -- make and compile template file
setwd("src")
library(TMB)
# library(TMBdebug)
library(TMBhelper)
version = "temporal_model_spde_v3"
if(data_list$debug == 1){
  #dyn.unload(version)
  #file.remove(paste0(version,".dll"))
  #file.remove(paste0(version,".o"))
}
TMB::compile( paste0(version,".cpp") )
dyn.load(dynlib(version))
setwd("../")

# -----------------------------------------------------------------------------------------------------------------------------
# STEP 2
# -----------------------------------------------------------------------------------------------------------------------------
# Fit 1 -- Backward model selection
map <- build_map(model, data_list, incl_disease = incl_disease, params, space = 0, space_time = 0)



###############################################
# Load packages
library(plyr)
library(gtools)
library(TMB)
library(TMBhelper)
library(foreach)
library(doParallel)

source("R/Support/8_missing_cont_val_update_v2.R")
source("R/Support/7_build_x_q_array.R")
source("R/Optimize_optim.R")

###############################################
# Set up parallel processing
numCores <- detectCores() - 1

###############################################
# Set up data objects
mod_list <- list()
terms_list <- list()
Obj_list <- list()
MSE_list <- list()

ind = 1
remove_vec <- c()
model.list = list()
summ_stats = data.frame(matrix(NA, ncol = 9 , nrow = 1))
colnames(summ_stats) = c("Model", "AIC", "BIC","Log_Lik","N_Params","Param", "Max_Grad", "Converg", "MSE")
Pred <- data.frame(Model0 = rep(NA, length(data_list$c_q)))


###############################################
# Fit intercept model: No betas estimated
map$beta_c <- factor(rep(NA, length(map$beta_c)))
map$beta_p <- factor(rep(NA, length(map$beta_p)))
map$beta_log <- factor(rep(NA, length(map$beta_log)))
map$x_q_cont_missing <- factor(rep(NA, length(map$x_q_cont_missing)))

# Set up data
data_list$x_q_hat_est <- as.numeric(as.character(map$beta_p))
new_map <- map_update( map , params , data_list )
new_dat <- expand_xq( data_list, map )
params$x_q_cont_missing <- rep(0, sum(is.na(new_dat$x_q[which(new_dat$fit_ll == 1),1:(new_dat$n_x_q_cont)])))


# Fit base model
obj = TMB::MakeADFun( data = new_dat, parameters = params, DLL = version, map = new_map, random = random, silent = TRUE)
opt = suppressWarnings( Optimize( obj , n = log(nrow(new_dat$x_q)), loopnum = 3, getsd = FALSE))

#
# ###############################################
# # MSE Run parallel
nobs <- length(data_list$c_q)

# Get parameters
last_par_full = obj$env$parList()
last_par_full_base <- last_par_full
last_par = obj$env$last.par.best[-obj$env$random]

# Run
pred_mat <- matrix(999, nrow = 1, ncol = 4)
for(i in 1:nobs){


  data_list$fit_ll[i] <- 0

  # Refit model
  new_map <- map_update( map , params , data_list )
  new_dat <- expand_xq( data_list, map )
  new_params <- last_par_full
  new_params$gamma_q[i] <- 0
  new_params$x_q_cont_missing <- rep(0, sum(is.na(new_dat$x_q[which(new_dat$fit_ll == 1),1:(new_dat$n_x_q_cont)])))

  # Fit base model
  new_obj = TMB::MakeADFun( data = new_dat, parameters = new_params, DLL = version, map = new_map, random = random, silent = TRUE)
  new_opt = tryCatch( Optimize( new_obj , startpar = last_par, n = log(nrow(data_list$x_q)) - 1, loopnum = 1, getsd = FALSE), error = function(e) NULL)

  # Turn obs back on
  data_list$fit_ll[i] = 1

  # Return prediction
  pred_vec <- new_obj$report(new_obj$env$last.par.best)$pred_cq
  pred_mat <- rbind(pred_mat, c(new_opt$AIC, new_opt$AICc, new_opt$BIC, pred_vec[i]))
  print(i)
}

# When you're done, clean up the cluster
MSE_list[[ind]] <- pred_mat

# Save model objects
# mod_list[[ind]] <- opt
# Obj_list[[ind]] <- obj
BIC <- 2 * opt$objective + log(nrow(data_list$x_q)) * length(opt$par); if(length(BIC) == 0){ BIC = Inf }
AIC <- opt$AIC; if(is.null(AIC)){ AIC = Inf }
summ_stats <- rbind(summ_stats, c("Base", AIC, BIC, ifelse(is.null(opt$objective), NA, opt$objective) ,
                                  ifelse(is.null(length(opt$par)), length(obj$par), length(opt$par)) ,
                                  "Intercept",
                                  ifelse(is.null(opt$max_gradient), NA, opt$max_gradient),
                                  ifelse(is.null( opt$Convergence_check ), "Did not converge", opt$Convergence_check)),
                    mean((pred_mat[,4] - data_list$c_q)^2))
# Save parameter estimates
rep = sdreport(obj)
coef_save_p <- round(rep$value[which(names(rep$value) %in% c("beta_p"))], 3)
coef_save_c <- round(rep$value[which(names(rep$value) %in% c("beta_c"))], 3)
terms_list[[ind]] <- rbind( coef_save_p, coef_save_c )
colnames(terms_list[[ind]]) <- colnames(data_list$x_q)

# Set lowest AIC to starting AIC
if(BIC_sel == T){
  lowest_val <- BIC
}
if( BIC_sel == F){
  lowest_val <- AIC
}

message( "#########################" )
message( paste0("Fit model ", ind) )
message( "#########################" )

# Run the next models
ind <- ind + 1
rm(opt)
rm(obj)


#---------------------------------------------------------------------------
# FORWARD SELECTION
# Loop model formulations to see if it improves AIC
for(i in 1:ncol(data_list$x_q)){

  # Change the map to turn parameters on
  map[["beta_p"]] <- as.character( map[["beta_p"]] )
  map[["beta_c"]] <- as.character( map[["beta_c"]] )
  map[["beta_c"]][i] <- i
  map[["beta_p"]][i] <- i
  map[["beta_c"]] <- factor(map[["beta_c"]])
  map[["beta_p"]] <- factor(map[["beta_p"]])

  # Set up data
  data_list$x_q_hat_est <- as.numeric(as.character(map$beta_p))
  new_map <- map_update( map , params , data_list )
  new_dat <- expand_xq( data_list, map )
  params <- last_par_full_base
  params$x_q_cont_missing <- rep(0, sum(is.na(new_dat$x_q[which(new_dat$fit_ll == 1),1:(new_dat$n_x_q_cont)])))

  # Fit base model
  obj = TMB::MakeADFun( data = new_dat, parameters = params, DLL = version, map = new_map, random = random, silent = TRUE)
  opt = suppressWarnings( Optimize( obj , n = log(nrow(new_dat$x_q)), loopnum = 3, getsd = FALSE))

  ###############################################
  # MSE Run parallel
  registerDoParallel(numCores)  # use multicore, set to the number of our cores
  nobs <- length(data_list$c_q)

  # Get parameters
  last_par_full = obj$env$parList()
  last_par = obj$env$last.par.best[-obj$env$random]

  # Run
  mse_pred <- foreach(i = 1:nobs, .combine=rbind) %dopar% {

    # Load packages
    library(plyr)
    library(gtools)
    library(TMB)
    library(TMBhelper)
    library(foreach)
    library(doParallel)
    library(TMB)
    # library(TMBdebug)
    library(TMBhelper)

    # Load functions
    source("R/Support/8_missing_cont_val_update_v2.R")
    source("R/Support/7_build_x_q_array.R")
    source("R/Optimize_optim.R")
    source("R/Support/1_read_dat_v3.R")
    source("R/Support/2_fit_mesh_v2.R")
    source("R/Support/3_build_data_v3.R")
    source("R/Support/4_build_params_v3.R")
    source("R/Support/5_build_map_v2.R")
    source("R/Support/bad_params.R")
    source("R/run_variogram.R")
    source("R/Support/aic_selection_v2.R")
    source("R/plot_coef_v2.R")

    setwd("src")
    version = "temporal_model_spde_v3"
    if(data_list$debug == 1){
      #dyn.unload(version)
      #file.remove(paste0(version,".dll"))
      #file.remove(paste0(version,".o"))
    }
    TMB::compile( paste0(version,".cpp") )
    dyn.load(dynlib(version))
    setwd("../")


    data_list$fit_ll[i] <- 0

    # Refit model
    new_map <- map_update( map , params , data_list )
    new_dat <- expand_xq( data_list, map )
    new_params <- last_par_full
    new_params$gamma_q[i] <- 0
    new_params$x_q_cont_missing <- rep(0, sum(is.na(new_dat$x_q[which(new_dat$fit_ll == 1),1:(new_dat$n_x_q_cont)])))

    # Fit base model
    new_obj = TMB::MakeADFun( data = new_dat, parameters = new_params, DLL = version, map = new_map, random = random, silent = TRUE)
    new_opt = tryCatch( Optimize( new_obj , startpar = last_par, n = log(nrow(data_list$x_q)) - 1, loopnum = 1, getsd = FALSE), error = function(e) NULL)

    # Turn obs back on
    data_list$fit_ll[i] = 1

    # Return prediction
    pred_vec <- new_obj$report(new_obj$env$last.par.best)$pred_cq
    c(new_opt$AIC, new_opt$AICc, new_opt$BIC, pred_vec[i])
  }

  # When you're done, clean up the cluster
  stopImplicitCluster()
  pred_mat <- mse_pred
  MSE_list[[ind]] <- mse_pred

  # Save model objects
  # mod_list[[ind]] <- opt
  # Obj_list[[ind]] <- obj
  BIC <- 2 * opt$objective + log(nrow(data_list$x_q)) * length(opt$par); if(length(BIC) == 0){ BIC = Inf }
  AIC <- opt$AIC; if(is.null(AIC)){ AIC = Inf }
  params_in_model <- paste(colnames(data_list$x_q)[which(!is.na(as.character( map[["beta_p"]] )))], collapse = ", ")
  summ_stats <- rbind(summ_stats,
                      c("Base", AIC, BIC, ifelse(is.null(opt$objective), NA, opt$objective) ,
                        ifelse(is.null(length(opt$par)), length(obj$par), length(opt$par)) ,
                        params_in_model,
                        ifelse(is.null(opt$max_gradient), NA, opt$max_gradient),
                        ifelse(is.null( opt$Convergence_check ), "Did not converge", opt$Convergence_check)),
                      mean((pred_mat[,4] - data_list$c_q)^2))

  if(BIC_sel == T){
    val_sel <- BIC
  }
  if( BIC_sel == F){
    val_sel <- AIC
  }

  # If model did not converge converged
  if(is.null(opt$SD$value) == T){
    coef_save <- rep("Did not converge", length(terms_list[[ind - 1]]))
    terms_list[[ind]] <- coef_save

    # Run the next models (turn off parameters if they didnt lead to convergence)
    map[["beta_p"]] <- as.character( map[["beta_p"]] )
    map[["beta_c"]] <- as.character( map[["beta_c"]] )
    map[["beta_c"]][i] <- NA
    map[["beta_p"]][i] <- NA
    map[["beta_c"]] <- factor(map[["beta_c"]])
    map[["beta_p"]] <- factor(map[["beta_p"]])
  }

  # If model converged
  if(is.null(opt$SD$value) == F) {
    rep = sdreport(obj)
    coef_save_p <- round(rep$value[which(names(rep$value) %in% c("beta_p"))], 3)
    coef_save_c <- round(rep$value[which(names(rep$value) %in% c("beta_c"))], 3)
    terms_list[[ind]] <- rbind( coef_save_p, coef_save_c )
    colnames(terms_list[[ind]]) <- colnames(data_list$x_q)

    # Turn off parameters if did not improve fit
    if(lowest_val < val_sel){
      map[["beta_p"]] <- as.character( map[["beta_p"]] )
      map[["beta_c"]] <- as.character( map[["beta_c"]] )
      map[["beta_c"]][i] <- NA
      map[["beta_p"]][i] <- NA
      map[["beta_c"]] <- factor(map[["beta_c"]])
      map[["beta_p"]] <- factor(map[["beta_p"]])
    }

    if(lowest_val > val_sel){
      lowest_val <- val_sel
    }
  }

  # save(ind, file = paste(ind, ".RData"))
  rm(opt)
  rm(obj)

  message( "#########################" )
  message( paste0("Fit model ", ind) )
  message( "#########################" )

  ind <- ind + 1
}

summ_stats <- rbind(summ_stats, rep("Backwards selection", ncol(summ_stats)))
#---------------------------------------------------------------------------
# BACKWARD SELECTION
# Loop model formulations to see if it improves AIC
params_left <- which(!is.na(as.character(map[["beta_p"]]))) # which parameters remain from forward selection
for(i in params_left){

  # Change the map to turn parameters off
  map[["beta_p"]] <- as.character( map[["beta_p"]] )
  map[["beta_c"]] <- as.character( map[["beta_c"]] )
  map[["beta_c"]][i] <- NA
  map[["beta_p"]][i] <- NA
  map[["beta_c"]] <- factor(map[["beta_c"]])
  map[["beta_p"]] <- factor(map[["beta_p"]])

  # Set up data
  data_list$x_q_hat_est <- as.numeric(as.character(map$beta_p))
  new_map <- map_update( map , params , data_list )
  new_dat <- expand_xq( data_list, map )
  params <- last_par_full_base
  params$x_q_cont_missing <- rep(0, sum(is.na(new_dat$x_q[which(new_dat$fit_ll == 1),1:(new_dat$n_x_q_cont)])))

  # Fit base model
  obj = TMB::MakeADFun( data = new_dat, parameters = params, DLL = version, map = new_map, random = random, silent = TRUE)
  opt = suppressWarnings( Optimize( obj , n = log(nrow(new_dat$x_q)), loopnum = 3, getsd = FALSE))


  ###############################################
  # MSE Run parallel
  registerDoParallel(numCores)  # use multicore, set to the number of our cores
  nobs <- length(data_list$c_q)

  # Get parameters
  last_par_full = obj$env$parList()
  last_par = obj$env$last.par.best[-obj$env$random]

  # Run
  mse_pred <- foreach(i = 1:nobs, .combine=rbind) %dopar% {

    # Load packages
    library(plyr)
    library(gtools)
    library(TMB)
    library(TMBhelper)
    library(foreach)
    library(doParallel)
    library(TMB)
    # library(TMBdebug)
    library(TMBhelper)

    # Load functions
    source("R/Support/8_missing_cont_val_update_v2.R")
    source("R/Support/7_build_x_q_array.R")
    source("R/Optimize_optim.R")
    source("R/Support/1_read_dat_v3.R")
    source("R/Support/2_fit_mesh_v2.R")
    source("R/Support/3_build_data_v3.R")
    source("R/Support/4_build_params_v3.R")
    source("R/Support/5_build_map_v2.R")
    source("R/Support/bad_params.R")
    source("R/run_variogram.R")
    source("R/Support/aic_selection_v2.R")
    source("R/plot_coef_v2.R")


    setwd("src")
    version = "temporal_model_spde_v3"
    if(data_list$debug == 1){
      #dyn.unload(version)
      #file.remove(paste0(version,".dll"))
      #file.remove(paste0(version,".o"))
    }
    TMB::compile( paste0(version,".cpp") )
    dyn.load(dynlib(version))
    setwd("../")

    data_list$fit_ll[i] <- 0

    # Refit model
    new_map <- map_update( map , params , data_list )
    new_dat <- expand_xq( data_list, map )
    new_params <- last_par_full
    new_params$gamma_q[i] <- 0
    new_params$x_q_cont_missing <- rep(0, sum(is.na(new_dat$x_q[which(new_dat$fit_ll == 1),1:(new_dat$n_x_q_cont)])))

    # Fit base model
    new_obj = TMB::MakeADFun( data = new_dat, parameters = new_params, DLL = version, map = new_map, random = random, silent = TRUE)
    new_opt = tryCatch( Optimize( new_obj , startpar = last_par, n = log(nrow(data_list$x_q)) - 1, loopnum = 1, getsd = FALSE), error = function(e) NULL)

    # Turn obs back on
    data_list$fit_ll[i] = 1

    # Return prediction
    pred_vec <- new_obj$report(new_obj$env$last.par.best)$pred_cq
    c(new_opt$AIC, new_opt$AICc, new_opt$BIC, pred_vec[i])
  }

  # When you're done, clean up the cluster
  stopImplicitCluster()
  pred_mat <- mse_pred
  MSE_list[[ind]] <- mse_pred

  # Save model objects
  # mod_list[[ind]] <- opt
  # Obj_list[[ind]] <- obj
  BIC <- 2 * opt$objective + log(nrow(data_list$x_q)) * length(opt$par); if(length(BIC) == 0){ BIC = Inf }
  AIC <- opt$AIC; if(is.null(AIC)){ AIC = Inf }
  params_in_model <- paste(colnames(data_list$x_q)[which(!is.na(as.character( map[["beta_p"]] )))], collapse = ", ")
  summ_stats <- rbind(summ_stats, c("Base", AIC, BIC, ifelse(is.null(opt$objective), NA, opt$objective) ,
                                    ifelse(is.null(length(opt$par)), length(obj$par), length(opt$par)) ,
                                    params_in_model,
                                    ifelse(is.null(opt$max_gradient), NA, opt$max_gradient),
                                    ifelse(is.null( opt$Convergence_check ), "Did not converge", opt$Convergence_check)),
                      mean((pred_mat[,4] - data_list$c_q)^2))

  if(BIC_sel == T){
    val_sel <- BIC
  }
  if( BIC_sel == F){
    val_sel <- AIC
  }

  # If did not converge
  if(is.null(opt$SD$value) == T){
    coef_save <- rep("Did not converge", length(terms_list[[ind - 1]]))
    terms_list[[ind]] <- coef_save

    # Turn off parameters
    map[["beta_p"]] <- as.character( map[["beta_p"]] )
    map[["beta_c"]] <- as.character( map[["beta_c"]] )
    map[["beta_c"]][i] <- i
    map[["beta_p"]][i] <- i
    map[["beta_c"]] <- factor(map[["beta_c"]])
    map[["beta_p"]] <- factor(map[["beta_p"]])
  }

  # If it converged
  if(is.null(opt$SD$value) == F) {
    rep = sdreport(obj)
    coef_save_p <- round(rep$value[which(names(rep$value) %in% c("beta_p"))], 3)
    coef_save_c <- round(rep$value[which(names(rep$value) %in% c("beta_c"))], 3)
    terms_list[[ind]] <- rbind( coef_save_p, coef_save_c )
    colnames(terms_list[[ind]]) <- colnames(data_list$x_q)

    # Turn off parameters if it does not improve fit
    if(lowest_val < val_sel){
      map[["beta_p"]] <- as.character( map[["beta_p"]] )
      map[["beta_c"]] <- as.character( map[["beta_c"]] )
      map[["beta_c"]][i] <- i
      map[["beta_p"]][i] <- i
      map[["beta_c"]] <- factor(map[["beta_c"]])
      map[["beta_p"]] <- factor(map[["beta_p"]])
    }

    # If it improves fit, leave on and update selection criterion
    if(lowest_val > val_sel){
      lowest_val <- val_sel
    }
  }

  rm(opt)
  rm(obj)

  message( "#########################" )
  message( paste0("Fit model ", ind) )
  message( "#########################" )
  # save(ind, file = paste(ind, ".RData"))
  ind <- ind + 1
}
#
# FIT THE FINAL MODEL

# Set up data
data_list$x_q_hat_est <- as.numeric(as.character(map$beta_p))
new_map <- map_update( map , params , data_list )
new_dat <- expand_xq( data_list, map )
params <- last_par_full_base
params$x_q_cont_missing <- rep(0, sum(is.na(new_dat$x_q[which(new_dat$fit_ll == 1),1:(new_dat$n_x_q_cont)])))


# Fit base model
obj = TMB::MakeADFun( data = new_dat, parameters = params, DLL = version, map = new_map, random = random, silent = TRUE)
opt = suppressWarnings( Optimize( obj , n = log(nrow(new_dat$x_q)), loopnum = 3, getsd = FALSE))

###############################################
# MSE Run parallel
registerDoParallel(numCores)  # use multicore, set to the number of our cores
nobs <- length(data_list$c_q)

# Get parameters
last_par_full = obj$env$parList()
last_par = obj$env$last.par.best[-obj$env$random]

# Run
mse_pred <- foreach(i = 1:nobs, .combine=rbind) %dopar% {

  # Load packages
  library(plyr)
  library(gtools)
  library(TMB)
  library(TMBhelper)
  library(foreach)
  library(doParallel)
  library(TMB)
  # library(TMBdebug)
  library(TMBhelper)

  # Load functions
  source("R/Support/8_missing_cont_val_update_v2.R")
  source("R/Support/7_build_x_q_array.R")
  source("R/Optimize_optim.R")
  source("R/Support/1_read_dat_v3.R")
  source("R/Support/2_fit_mesh_v2.R")
  source("R/Support/3_build_data_v3.R")
  source("R/Support/4_build_params_v3.R")
  source("R/Support/5_build_map_v2.R")
  source("R/Support/bad_params.R")
  source("R/run_variogram.R")
  source("R/Support/aic_selection_v2.R")
  source("R/plot_coef_v2.R")

  setwd("src")
  version = "temporal_model_spde_v3"
  if(data_list$debug == 1){
    #dyn.unload(version)
    #file.remove(paste0(version,".dll"))
    #file.remove(paste0(version,".o"))
  }
  TMB::compile( paste0(version,".cpp") )
  dyn.load(dynlib(version))
  setwd("../")


  data_list$fit_ll[i] <- 0

  # Refit model
  new_map <- map_update( map , params , data_list )
  new_dat <- expand_xq( data_list, map )
  new_params <- last_par_full
  new_params$gamma_q[i] <- 0
  new_params$x_q_cont_missing <- rep(0, sum(is.na(new_dat$x_q[which(new_dat$fit_ll == 1),1:(new_dat$n_x_q_cont)])))

  # Fit base model
  new_obj = TMB::MakeADFun( data = new_dat, parameters = new_params, DLL = version, map = new_map, random = random, silent = TRUE)
  new_opt = tryCatch( Optimize( new_obj , startpar = last_par, n = log(nrow(data_list$x_q)) - 1, loopnum = 1, getsd = FALSE), error = function(e) NULL)

  # Turn obs back on
  data_list$fit_ll[i] = 1

  # Return prediction
  pred_vec <- new_obj$report(new_obj$env$last.par.best)$pred_cq
  c(new_opt$AIC, new_opt$AICc, new_opt$BIC, pred_vec[i])
}

# When you're done, clean up the cluster
stopImplicitCluster()
pred_mat <- mse_pred
MSE_list[[ind]] <- mse_pred

# Save model objects
mod_list[[ind]] <- opt
Obj_list[[ind]] <- obj
BIC <- 2 * opt$objective + log(nrow(data_list$x_q)) * length(opt$par); if(length(BIC) == 0){ BIC = Inf }
AIC <- opt$AIC; if(is.null(AIC)){ AIC = Inf }
params_in_model <- paste(colnames(data_list$x_q)[which(!is.na(as.character( map[["beta_p"]] )))], collapse = ", ")
summ_stats <- rbind(summ_stats, c("Base", AIC, BIC, ifelse(is.null(opt$objective), NA, opt$objective) ,
                                  ifelse(is.null(length(opt$par)), length(obj$par), length(opt$par)) ,
                                  params_in_model,
                                  ifelse(is.null(opt$max_gradient), NA, opt$max_gradient),
                                  ifelse(is.null( opt$Convergence_check ), "Did not converge", opt$Convergence_check)),
                    mean((pred_mat[,4] - data_list$c_q)^2))

# Save parameter estimates
if(is.null(opt$SD$value) == T){
  coef_save <- rep("Did not converge", length(terms_list[[ind - 1]]))
  terms_list[[ind]] <- coef_save
}
if(is.null(opt$SD$value) == F) {
  rep = sdreport(obj)
  coef_save_p <- round(rep$value[which(names(rep$value) == c("beta_p"))], 3)
  coef_save_c <- round(rep$value[which(names(rep$value) == c("beta_c"))], 3)
  terms_list[[ind]] <- rbind( coef_save_p, coef_save_c )
  colnames(terms_list[[ind]]) <- colnames(data_list$x_q)
}

message( "#########################" )
message( paste0("Fit final model") )
message( "#########################" )

# Return results
results_list = list(summ_stats = summ_stats, mod_list = mod_list, terms_list = terms_list, Obj_list = Obj_list, MSE_list = MSE_list)
save(results_list, file = "MSE_Run.RData")
return(results_list)
