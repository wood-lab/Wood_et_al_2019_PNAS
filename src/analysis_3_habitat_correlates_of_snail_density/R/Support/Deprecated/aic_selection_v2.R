# Function to do backward model selection using p-values, aic, or bic

aic_sel <- function( data_list = data_list, params = params, version = version, map = map, random = random , BIC_sel = T ){
  library(plyr)
  library(gtools)
  library(TMB)
  library(TMBhelper)
  source("R/Support/missing_cont_val_update_v2.R")
  source("R/Support/7_build_x_q_array.R")

  mod_list <- list()
  terms_list <- list()
  obj_list <- list()

  ind = 1
  remove_vec <- c()
  model.list = list()
  summ_stats = data.frame(matrix(NA, ncol = 6 , nrow = 1))
  colnames(summ_stats) = c("Model", "AIC", "BIC","Log_Lik","N_Params","Param")

  # Fit intercept model: No betas estimated
  map$beta_c <- factor(rep(NA, length(map$beta_c)))
  map$beta_p <- factor(rep(NA, length(map$beta_p)))
  map$beta_log <- factor(rep(NA, length(map$beta_log)))
  map$x_q_cont_missing <- factor(rep(NA, length(map$x_q_cont_missing)))

  data_list$x_q_hat_est <- as.numeric(as.character(map$beta_p))
  obj = TMB::MakeADFun( data = expand_xq( data_list, map ), parameters = params, DLL = version, map = map_update( map , params , data_list ), random = random)
  opt = tryCatch(optim(obj$par, obj$fn, obj$gr, control = list(maxit = 100000)), error = function(e) NULL)

  # Save model objects
  mod_list[[ind]] <- opt
  obj_list[[ind]] <- obj
  BIC <- 2 * opt$value + log(nrow(data_list$x_q)) * length(opt$par)
  AIC <- 2 * opt$value + 2 * length(opt$par); if(is.null(AIC)){ AIC = 1e6 }
  summ_stats <- rbind(summ_stats, (c("Base", AIC, BIC, opt$value,  length(opt$par), NA)))

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


  # Run the next models
  ind <- ind + 1

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

    # Turn on logistic parameters for missing data
    # if(i > data_list$n_x_q_cont){
    #   map[["beta_log"]] <- as.character( map[["beta_log"]] )
    #   map[["beta_log"]][i] <- i
    #   map[["beta_log"]] <- factor(map[["beta_log"]])
    # }

    data_list$x_q_hat_est <- as.numeric(as.character(map$beta_p))

    # Fit model
    obj = TMB::MakeADFun( data = expand_xq( data_list, map ), parameters = params, DLL = version, map = map_update( map , params , data_list ), random = random)
opt = tryCatch(optim(obj$par, obj$fn, obj$gr, control = list(maxit = 100000)), error = function(e) NULL)

    # Save model objects
    mod_list[[ind]] <- opt
    obj_list[[ind]] <- obj
    BIC <- 2 * opt$value + log(nrow(data_list$x_q)) * length(opt$par)
    AIC <- 2 * opt$value + 2 * length(opt$par); if(is.null(AIC)){ AIC = NA }
    params_in_model <- paste(colnames(data_list$x_q)[which(!is.na(as.character( map[["beta_p"]] )))], collapse = ", ")
    summ_stats <- rbind(summ_stats, (c("Base", AIC, BIC, opt$value,  length(opt$par), params_in_model)))

    if(BIC_sel == T){
      val_sel <- BIC
    }
    if( BIC_sel == F){
      val_sel <- AIC
    }

    # Save parameter estimates
    if(is.null(opt) == TRUE | opt$convergence > 0){
      coef_save <- rep("Did not converge", length(terms_list[[ind - 1]]))
      terms_list[[ind]] <- coef_save
    }
    if(is.null(opt) == FALSE) {
      rep = sdreport(obj)
      coef_save_p <- round(rep$value[which(names(rep$value) %in% c("beta_p"))], 3)
      coef_save_c <- round(rep$value[which(names(rep$value) %in% c("beta_c"))], 3)
      terms_list[[ind]] <- rbind( coef_save_p, coef_save_c )
      colnames(terms_list[[ind]]) <- colnames(data_list$x_q)
    }

    # Run the next models (turn of parameters if they dont do anything)
    if(is.null(opt) == TRUE | opt$convergence > 0){
      map[["beta_p"]] <- as.character( map[["beta_p"]] )
      map[["beta_c"]] <- as.character( map[["beta_c"]] )
      map[["beta_c"]][i] <- NA
      map[["beta_p"]][i] <- NA
      map[["beta_c"]] <- factor(map[["beta_c"]])
      map[["beta_p"]] <- factor(map[["beta_p"]])

      # Turn off logistic parameters for missing data
      # if(i > data_list$n_x_q_cont){
      #   map[["beta_log"]] <- as.character( map[["beta_log"]] )
      #   map[["beta_log"]][i] <- NA
      #   map[["beta_log"]] <- factor(map[["beta_log"]])
      # }
    }
    if(is.null(opt) == FALSE) {
      if(lowest_val < val_sel){
        map[["beta_p"]] <- as.character( map[["beta_p"]] )
        map[["beta_c"]] <- as.character( map[["beta_c"]] )
        map[["beta_c"]][i] <- NA
        map[["beta_p"]][i] <- NA
        map[["beta_c"]] <- factor(map[["beta_c"]])
        map[["beta_p"]] <- factor(map[["beta_p"]])

        # Turn on logistic parameters for missing data
        # if(i > data_list$n_x_q_cont){
        #   map[["beta_log"]] <- as.character( map[["beta_log"]] )
        #   map[["beta_log"]][i] <- NA
        #   map[["beta_log"]] <- factor(map[["beta_log"]])
        # }
      }
    }

    if(is.null(opt) == FALSE){
      if(lowest_val > val_sel){
        lowest_val <- val_sel
      }
    }
    # save(ind, file = paste(ind, ".RData"))
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

    # Turn on logistic parameters for missing data
    # if(i > data_list$n_x_q_cont){
    #   map[["beta_log"]] <- as.character( map[["beta_log"]] )
    #   map[["beta_log"]][i] <- NA
    #   map[["beta_log"]] <- factor(map[["beta_log"]])
    # }

    data_list$x_q_hat_est <- as.numeric(as.character(map$beta_p))

    # Fit model
    obj = TMB::MakeADFun( data = expand_xq( data_list, map ), parameters = params, DLL = version, map = map_update( map , params , data_list ), random = random)
    opt = tryCatch(optim(obj$par, obj$fn, obj$gr, control = list(maxit = 100000)), error = function(e) NULL)

    # Save model objects
    mod_list[[ind]] <- opt
    obj_list[[ind]] <- obj
    BIC <- 2 * opt$value + log(nrow(data_list$x_q)) * length(opt$par)
    AIC <- 2 * opt$value + 2 * length(opt$par); if(is.null(AIC)){ AIC = NA }
    params_in_model <- paste(colnames(data_list$x_q)[which(!is.na(as.character( map[["beta_p"]] )))], collapse = ", ")
    summ_stats <- rbind(summ_stats, (c("Base", AIC, BIC, opt$value,  length(opt$par), params_in_model)))

    if(BIC_sel == T){
      val_sel <- BIC
    }
    if( BIC_sel == F){
      val_sel <- AIC
    }

    # Save parameter estimates
    if(is.null(opt) == TRUE | opt$convergence > 0){
      coef_save <- rep("Did not converge", length(terms_list[[ind - 1]]))
      terms_list[[ind]] <- coef_save
    } else if(is.null(opt) == FALSE) {
      rep = sdreport(obj)
      coef_save_p <- round(rep$value[which(names(rep$value) %in% c("beta_p"))], 3)
      coef_save_c <- round(rep$value[which(names(rep$value) %in% c("beta_c"))], 3)
      terms_list[[ind]] <- rbind( coef_save_p, coef_save_c )
      colnames(terms_list[[ind]]) <- colnames(data_list$x_q)
    }

    # Run the next models: Turn parameters back on if the new model didn't converge
    if(is.null(opt) == TRUE | opt$convergence > 0){
      map[["beta_p"]] <- as.character( map[["beta_p"]] )
      map[["beta_c"]] <- as.character( map[["beta_c"]] )
      map[["beta_c"]][i] <- i
      map[["beta_p"]][i] <- i
      map[["beta_c"]] <- factor(map[["beta_c"]])
      map[["beta_p"]] <- factor(map[["beta_p"]])

      # Turn off logistic parameters for missing data
      # if(i > data_list$n_x_q_cont){
      #   map[["beta_log"]] <- as.character( map[["beta_log"]] )
      #   map[["beta_log"]][i] <- i
      #   map[["beta_log"]] <- factor(map[["beta_log"]])
      # }
    }
    if(is.null(opt) == FALSE) {
      if(lowest_val < val_sel){
        map[["beta_p"]] <- as.character( map[["beta_p"]] )
        map[["beta_c"]] <- as.character( map[["beta_c"]] )
        map[["beta_c"]][i] <- i
        map[["beta_p"]][i] <- i
        map[["beta_c"]] <- factor(map[["beta_c"]])
        map[["beta_p"]] <- factor(map[["beta_p"]])

        # Turn on logistic parameters for missing data
        # if(i > data_list$n_x_q_cont){
        #   map[["beta_log"]] <- as.character( map[["beta_log"]] )
        #   map[["beta_log"]][i] <- i
        #   map[["beta_log"]] <- factor(map[["beta_log"]])
        # }
      }
    }

    if(is.null(opt) == FALSE){
      if(lowest_val > val_sel){
        lowest_val <- val_sel
      }
    }
    # save(ind, file = paste(ind, ".RData"))
    ind <- ind + 1
  }

  # FIT THE FINAL MODEL
  data_list$x_q_hat_est <- as.numeric(as.character(map$beta_p))
  obj = TMB::MakeADFun( data = expand_xq( data_list, map ), parameters = params, DLL = version, map = map_update( map , params , data_list ), random = random)
  opt = tryCatch(optim(obj$par, obj$fn, obj$gr, control = list(maxit = 100000)), error = function(e) NULL)

  # Save model objects
  mod_list[[ind]] <- opt
  obj_list[[ind]] <- obj
  BIC <- 2 * opt$value + log(nrow(data_list$x_q)) * length(opt$par)
  AIC <- 2 * opt$value + 2 * length(opt$par); if(is.null(AIC)){ AIC = NA }
  params_in_model <- paste(colnames(data_list$x_q)[which(!is.na(as.character( map[["beta_p"]] )))], collapse = ", ")
  summ_stats <- rbind(summ_stats, (c("Base", AIC, BIC, opt$value,  length(opt$par), params_in_model)))

  # Save parameter estimates
  if(is.null(opt) == TRUE | opt$convergence > 0){
    coef_save <- rep("Did not converge", length(terms_list[[ind - 1]]))
    terms_list[[ind]] <- coef_save
  }
  if(is.null(opt) == FALSE) {
    rep = sdreport(obj)
    coef_save_p <- round(rep$value[which(names(rep$value) == c("beta_p"))], 3)
    coef_save_c <- round(rep$value[which(names(rep$value) == c("beta_c"))], 3)
    terms_list[[ind]] <- rbind( coef_save_p, coef_save_c )
    colnames(terms_list[[ind]]) <- colnames(data_list$x_q)
  }


  # Return results
  results_list = list(summ_stats = summ_stats, mod_list = mod_list, terms_list = terms_list, obj_list = obj_list, final_obj = obj, final_opt = opt)
  return(results_list)
}
