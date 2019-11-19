# Function to do backward model selection using p-values, aic, or bic

aic_sel <- function( data_list = data_list, params = params, version = version, map = map, random = random , BIC_sel = T, method = "Nelder-Mead", silent = TRUE ){

  library(plyr)
  library(gtools)
  library(TMB)
  library(TMBhelper)
  source("R/Support/missing_cont_val_update_v2.R")
  source("R/Support/7_build_x_q_array.R")
  source("R/Optimize_optim.R")

  mod_list <- list()
  terms_list <- list()
  Obj_list <- list()

  ind = 1
  remove_vec <- c()
  model.list = list()
  summ_stats = data.frame(matrix(NA, ncol = 8 , nrow = 1))
  colnames(summ_stats) = c("Model", "AIC", "BIC","Log_Lik","N_Params","Param", "Max_Grad", "Converg")

  # Fit intercept model: No betas estimated
  map$beta_c <- factor(rep(NA, length(map$beta_c)))
  map$beta_p <- factor(rep(NA, length(map$beta_p)))
  map$beta_log <- factor(rep(NA, length(map$beta_log)))
  map$x_q_cont_missing <- factor(rep(NA, length(map$x_q_cont_missing)))

  # Fit
  data_list$x_q_hat_est <- as.numeric(as.character(map$beta_p))
  Obj = TMB::MakeADFun( data = expand_xq( data_list, map ), parameters = params, DLL = version, map = map_update( map , params , data_list ), random = random, silent = silent)
  Opt = tryCatch( Optimize( Obj , n = log(nrow(data_list$x_q)), loopnum = 6), error = function(e) NULL)


  # Save model objects
  # mod_list[[ind]] <- Opt
  # Obj_list[[ind]] <- Obj
  BIC <- 2 * Opt$objective + log(nrow(data_list$x_q)) * length(Opt$par); if(length(BIC) == 0){ BIC = Inf }
  AIC <- Opt$AIC; if(is.null(AIC)){ AIC = Inf }
  summ_stats <- rbind(summ_stats, c("Base", AIC, BIC, ifelse(is.null(Opt$objective), NA, Opt$objective) ,
                                    ifelse(is.null(length(Opt$par)), length(Obj$par), length(Opt$par)) ,
                                    "Intercept",
                                    ifelse(is.null(Opt$max_gradient), NA, Opt$max_gradient),
                                    ifelse(is.null( Opt$Convergence_check ), "Did not converge", Opt$Convergence_check)))
  # Save parameter estimates
  rep = sdreport(Obj)
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
  rm(Opt)
  rm(Obj)

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


    data_list$x_q_hat_est <- as.numeric(as.character(map$beta_p))

    # Fit model
    Obj = TMB::MakeADFun( data = expand_xq( data_list, map ), parameters = params, DLL = version, map = map_update( map , params , data_list ), random = random, silent = silent)
    Opt = tryCatch( Optimize( Obj , n = log(nrow(data_list$x_q)), loopnum = 6), error = function(e) NULL)

    # Save model objects
    # mod_list[[ind]] <- Opt
    # Obj_list[[ind]] <- Obj
    BIC <- 2 * Opt$objective + log(nrow(data_list$x_q)) * length(Opt$par); if(length(BIC) == 0){ BIC = Inf }
    AIC <- Opt$AIC; if(is.null(AIC)){ AIC = Inf }
    params_in_model <- paste(colnames(data_list$x_q)[which(!is.na(as.character( map[["beta_p"]] )))], collapse = ", ")
    summ_stats <- rbind(summ_stats, c("Base", AIC, BIC, ifelse(is.null(Opt$objective), NA, Opt$objective) ,
                                      ifelse(is.null(length(Opt$par)), length(Obj$par), length(Opt$par)) ,
                                      params_in_model,
                                      ifelse(is.null(Opt$max_gradient), NA, Opt$max_gradient),
                                      ifelse(is.null( Opt$Convergence_check ), "Did not converge", Opt$Convergence_check)))

    if(BIC_sel == T){
      val_sel <- BIC
    }
    if( BIC_sel == F){
      val_sel <- AIC
    }

    # If model did not converge converged
    if(is.null(Opt$SD$value) == T){
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
    if(is.null(Opt$SD$value) == F) {
      rep = sdreport(Obj)
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
    rm(Opt)
    rm(Obj)

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


    data_list$x_q_hat_est <- as.numeric(as.character(map$beta_p))

    # Fit model
    Obj = TMB::MakeADFun( data = expand_xq( data_list, map ), parameters = params, DLL = version, map = map_update( map , params , data_list ), random = random, silent = silent)
    Opt = tryCatch( Optimize( Obj , n = log(nrow(data_list$x_q)), loopnum = 6), error = function(e) NULL)

    # Save model objects
    # mod_list[[ind]] <- Opt
    # Obj_list[[ind]] <- Obj
    BIC <- 2 * Opt$objective + log(nrow(data_list$x_q)) * length(Opt$par); if(length(BIC) == 0){ BIC = Inf }
    AIC <- Opt$AIC; if(is.null(AIC)){ AIC = Inf }
    params_in_model <- paste(colnames(data_list$x_q)[which(!is.na(as.character( map[["beta_p"]] )))], collapse = ", ")
    summ_stats <- rbind(summ_stats, c("Base", AIC, BIC, ifelse(is.null(Opt$objective), NA, Opt$objective) ,
                                      ifelse(is.null(length(Opt$par)), length(Obj$par), length(Opt$par)) ,
                                      params_in_model,
                                      ifelse(is.null(Opt$max_gradient), NA, Opt$max_gradient),
                                      ifelse(is.null( Opt$Convergence_check ), "Did not converge", Opt$Convergence_check)))

    if(BIC_sel == T){
      val_sel <- BIC
    }
    if( BIC_sel == F){
      val_sel <- AIC
    }

    # If did not converge
    if(is.null(Opt$SD$value) == T){
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
    if(is.null(Opt$SD$value) == F) {
      rep = sdreport(Obj)
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

    rm(Opt)
    rm(Obj)

    message( "#########################" )
    message( paste0("Fit model ", ind) )
    message( "#########################" )
    # save(ind, file = paste(ind, ".RData"))
    ind <- ind + 1
  }

  # FIT THE FINAL MODEL
  data_list$x_q_hat_est <- as.numeric(as.character(map$beta_p))
  Obj = TMB::MakeADFun( data = expand_xq( data_list, map ), parameters = params, DLL = version, map = map_update( map , params , data_list ), random = random, silent = silent)
  Opt = tryCatch( Optimize( Obj , n = log(nrow(data_list$x_q)), loopnum = 6), error = function(e) NULL)

  # Save model objects
  mod_list[[ind]] <- Opt
  Obj_list[[ind]] <- Obj
  BIC <- 2 * Opt$objective + log(nrow(data_list$x_q)) * length(Opt$par); if(length(BIC) == 0){ BIC = Inf }
  AIC <- Opt$AIC; if(is.null(AIC)){ AIC = Inf }
  params_in_model <- paste(colnames(data_list$x_q)[which(!is.na(as.character( map[["beta_p"]] )))], collapse = ", ")
  summ_stats <- rbind(summ_stats, c("Base", AIC, BIC, ifelse(is.null(Opt$objective), NA, Opt$objective) ,
                                    ifelse(is.null(length(Opt$par)), length(Obj$par), length(Opt$par)) ,
                                    params_in_model,
                                    ifelse(is.null(Opt$max_gradient), NA, Opt$max_gradient),
                                    ifelse(is.null( Opt$Convergence_check ), "Did not converge", Opt$Convergence_check)))

  # Save parameter estimates
  if(is.null(Opt$SD$value) == T){
    coef_save <- rep("Did not converge", length(terms_list[[ind - 1]]))
    terms_list[[ind]] <- coef_save
  }
  if(is.null(Opt$SD$value) == F) {
    rep = sdreport(Obj)
    coef_save_p <- round(rep$value[which(names(rep$value) == c("beta_p"))], 3)
    coef_save_c <- round(rep$value[which(names(rep$value) == c("beta_c"))], 3)
    terms_list[[ind]] <- rbind( coef_save_p, coef_save_c )
    colnames(terms_list[[ind]]) <- colnames(data_list$x_q)
  }

  message( "#########################" )
  message( paste0("Fit final model") )
  message( "#########################" )

  # Return results
  results_list = list(summ_stats = summ_stats, mod_list = mod_list, terms_list = terms_list, Obj_list = Obj_list)
  return(results_list)
}
