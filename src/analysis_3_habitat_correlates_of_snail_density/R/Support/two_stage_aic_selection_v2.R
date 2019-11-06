# Function to do forward and backward model selection using p-values, aic, or bic. Iteratively evaluates models on binomial and then positive density part.

aic_sel_two <- function( data_list = data_list, params = params, version = version, map = map, random = random , BIC_sel = T, bulinus, wd = getwd(), file_name = "NULL", two_stage = FALSE){
  source("R/Support/missing_cont_val_update_v2.R")
  source("R/Support/7_build_x_q_array.R")

  mod_list <- list()
  terms_list <- list()
  Obj_list <- list()

  method <- "BFGS"

  ind = 1
  remove_vec <- c()
  model.list = list()
  summ_stats = data.frame(matrix(NA, ncol = 7 , nrow = 1))
  colnames(summ_stats) = c("Model", "AIC", "BIC","Log_Lik","N_Params","P_params", "C_params")

  # Fit intercept model: No betas estimated
  map$beta_c <- factor(rep(NA, length(map$beta_c)))
  map$beta_p <- factor(rep(NA, length(map$beta_p)))
  map$beta_log <- factor(rep(NA, length(map$beta_log)))
  map$x_q_cont_missing <- factor(rep(NA, length(map$x_q_cont_missing)))

  data_list$x_q_hat_est <- as.numeric(as.character(map$beta_p))
  obj = TMB::MakeADFun( data = expand_xq( data_list, map ), parameters = params, DLL = version, map = map_update( map , params , data_list ), random = random, silent = TRUE)
  opt = tryCatch(optim(obj$par, obj$fn, obj$gr, control = list(maxit = 100000), method = method), error = function(e) NULL)

  # If it converged
  if(!is.null(opt)){
    # Refit
    # last_par = obj$env$parList(opt$par)
    # obj = TMB::MakeADFun(data = expand_xq( data_list, map ), parameters = last_par, DLL = version, map = map_update( map , params , data_list ), random = random, silent = TRUE)
    # opt = tryCatch(TMBhelper::Optimize( obj ), error = function(e) NULL)

    # Save model objects
    BIC <- 2 * opt$value + log(nrow(data_list$x_q)) * length(opt$par)
    AIC <- 2 * opt$value + 2 * length(opt$par); if(is.null(AIC)){ AIC = 1e6 }
    summ_stats <- rbind(summ_stats, (c("Base", AIC, BIC, opt$value,  length(opt$par), NA, NA)))

    # Set lowest AIC to starting AIC
    if(BIC_sel == T){
      lowest_val <- BIC
    }
    if( BIC_sel == F){
      lowest_val <- AIC
    }

    # Save parameter estimates
    rep = TMB::sdreport(obj)
    coef_save_p <- round(rep$value[which(names(rep$value) %in% c("beta_p"))], 3)
    coef_save_c <- round(rep$value[which(names(rep$value) %in% c("beta_c"))], 3)
    terms_list[[ind]] <- rbind( coef_save_p, coef_save_c )
    colnames(terms_list[[ind]]) <- colnames(data_list$x_q)
  }


  # If it did not converge
  if(is.null(opt)){

    # Save model objects
    BIC <- Inf
    AIC <- Inf
    summ_stats <- rbind(summ_stats, (c("Base", AIC, BIC, opt$value,  length(opt$par), NA, NA)))

    # Set lowest AIC to starting AIC
    lowest_val <- Inf

    # Save parameter estimates
    coef_save <- rep("Did not converge", ncol(data_list$x_q))
    terms_list[[ind]] <- coef_save
    colnames(terms_list[[ind]]) <- colnames(data_list$x_q)
  }

  # Run the next models
  ind <- ind + 1

  params_loop <- c("beta_p","beta_c")
  #---------------------------------------------------------------------------
  # FORWARD SELECTION
  # Loop model formulations to see if it improves AIC
  for(i in 1:ncol(data_list$x_q)){
    for(j in 1:2){

      # Change the map to turn binomial parameters on
      map[[ params_loop[j] ]] <- as.character( map[[ params_loop[j] ]] )
      map[[ params_loop[j] ]][i] <- i
      map[[ params_loop[j] ]] <- factor(map[[ params_loop[j] ]])

      # Turn on estimation of missing values
      data_list$x_q_hat_est <- as.numeric(as.character(map$beta_p))
      data_list$x_q_hat_est[which(is.na(data_list$x_q_hat_est))] <- as.numeric(as.character(map$beta_c))[which(is.na(data_list$x_q_hat_est))]

      # Fit model
      obj = TMB::MakeADFun( data = expand_xq( data_list, map ), parameters = params, DLL = version, map = map_update( map , params , data_list ), random = random, silent = TRUE)
      opt = tryCatch(optim(obj$par, obj$fn, obj$gr, control = list(maxit = 100000)), error = function(e) NULL)


      # Save parameter names
      p_params_in_model <- paste(colnames(data_list$x_q)[which(!is.na(as.character( map[["beta_p"]] )))], collapse = ", ")
      c_params_in_model <- paste(colnames(data_list$x_q)[which(!is.na(as.character( map[["beta_c"]] )))], collapse = ", ")

      # If the model did not converge
      if(is.null(opt) == T | opt$convergence > 0){
        # Save values
        coef_save <- rep("Did not converge", length(terms_list[[ind - 1]]))
        terms_list[[ind]] <- coef_save

        summ_stats <- rbind(summ_stats, (c(ind, opt$convergence, "Did not converge", "Did not converge",  length(opt$par), p_params_in_model, c_params_in_model)))

        # Run the next models (turn of parameters if they dont do anything)
        map[[params_loop[j]]] <- as.character( map[[params_loop[j]]] )
        map[[params_loop[j]]][i] <- NA
        map[[params_loop[j]]] <- factor(map[[params_loop[j]]])
        print(paste0("Model ", ind, "did not converge"))
      }


      # If model converges
      else if(is.null(opt) == F) {

        # # Refit
        # last_par = obj$env$parList(opt$par)
        # obj = TMB::MakeADFun(data = expand_xq( data_list, map ), parameters = last_par, DLL = version, map = map_update( map , params , data_list ), random = random, silent = TRUE)
        # opt = tryCatch(TMBhelper::Optimize( obj ), error = function(e) NULL)

        # Save output
        BIC <- 2 * opt$value + log(nrow(data_list$x_q)) * length(opt$par)
        AIC <- 2 * opt$value + 2 * length(opt$par); if(is.null(AIC)){ AIC = NA }
        summ_stats <- rbind(summ_stats, (c(ind, AIC, BIC, opt$value,  length(opt$par), p_params_in_model, c_params_in_model)))

        if(BIC_sel == T){
          val_sel <- BIC
        }
        if( BIC_sel == F){
          val_sel <- AIC
        }

        # Save parameter estimates
        rep = TMB::sdreport(obj)
        coef_save_p <- round(rep$value[which(names(rep$value) %in% c("beta_p"))], 3)
        coef_save_c <- round(rep$value[which(names(rep$value) %in% c("beta_c"))], 3)
        terms_list[[ind]] <- rbind( coef_save_p, coef_save_c )
        colnames(terms_list[[ind]]) <- colnames(data_list$x_q)

        # Turn off parameters if they did not improve model fit
        if(lowest_val < val_sel){
          map[[params_loop[j]]] <- as.character( map[[params_loop[j]]] )
          map[[params_loop[j]]][i] <- NA
          map[[params_loop[j]]] <- factor(map[[params_loop[j]]])
        }

        # Replace the lowest IC if parameter improves fit
        if(lowest_val > val_sel){
          lowest_val <- val_sel
        }
        print(paste0("Model ", ind, "converged"))
      }

      # Report
      # save(ind, file = paste0(wd,  "/ind_",ind, "-j_",j,"-i_",i, file_name, ".RData"))
      ind <- ind + 1
      rm(opt)
    }
  }



  #---------------------------------------------------------------------------
  # BACKWARD SELECTION
  summ_stats <- rbind(summ_stats, rep("Backwards selection", ncol(summ_stats)))
  # Loop model formulations to see if it improves AIC
  params_left <- rbind( as.character(map[["beta_p"]]), as.character(map[["beta_c"]]))  # which parameters remain from forward selection
  for(i in 1:ncol(params_left)){
    for(j in 1:nrow(params_left)){
      if(!is.na(params_left[ j, i ])){

        # Change the map to turn parameters off
        map[[ params_loop[j] ]] <- as.character( map[[ params_loop[j] ]] )
        map[[ params_loop[j] ]][i] <- NA
        map[[ params_loop[j] ]] <- factor(map[[ params_loop[j] ]])

        # Turn on estimation of missing values
        data_list$x_q_hat_est <- as.numeric(as.character(map$beta_p))
        data_list$x_q_hat_est[which(is.na(data_list$x_q_hat_est))] <- as.numeric(as.character(map$beta_c))[which(is.na(data_list$x_q_hat_est))]

        # Fit model
        obj = TMB::MakeADFun( data = expand_xq( data_list, map ), parameters = params, DLL = version, map = map_update( map , params , data_list ), random = random, silent = TRUE)
        opt <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 100000))

        # Save parameter names
        p_params_in_model <- paste(colnames(data_list$x_q)[which(!is.na(as.character( map[["beta_p"]] )))], collapse = ", ")
        c_params_in_model <- paste(colnames(data_list$x_q)[which(!is.na(as.character( map[["beta_c"]] )))], collapse = ", ")


        # If the model did not converge
        if(is.null(opt) == T | opt$convergence > 0){
          # Save values
          coef_save <- rep("Did not converge", length(terms_list[[ind - 1]]))
          terms_list[[ind]] <- coef_save

          summ_stats <- rbind(summ_stats, (c(ind, opt$convergence, "Did not converge", "Did not converge",  length(opt$par), p_params_in_model, c_params_in_model)))

          # Turn on parameters if removal led to non-convergence
          map[[params_loop[j]]] <- as.character( map[[params_loop[j]]] )
          map[[params_loop[j]]][i] <- i
          map[[params_loop[j]]] <- factor(map[[params_loop[j]]])
          print(paste0("Model ", ind, "did not converge"))
        }


        # If model converges
        else if(is.null(opt) == F)  {

          # # Refit
          # last_par = obj$env$parList(opt$par)
          # obj = TMB::MakeADFun(data = expand_xq( data_list, map ), parameters = last_par, DLL = version, map = map_update( map , params , data_list ), random = random, silent = TRUE)
          # opt = tryCatch(TMBhelper::Optimize( obj ), error = function(e) NULL)

          # Save output
          BIC <- 2 * opt$value + log(nrow(data_list$x_q)) * length(opt$par)
          AIC <- 2 * opt$value + 2 * length(opt$par); if(is.null(AIC)){ AIC = NA }
          summ_stats <- rbind(summ_stats, (c(ind, AIC, BIC, opt$value,  length(opt$par), p_params_in_model, c_params_in_model)))
          if(BIC_sel == T){
            val_sel <- BIC
          }
          if( BIC_sel == F){
            val_sel <- AIC
          }

          # Save parameter estimates
          rep = TMB::sdreport(obj)
          coef_save_p <- round(rep$value[which(names(rep$value) %in% c("beta_p"))], 3)
          coef_save_c <- round(rep$value[which(names(rep$value) %in% c("beta_c"))], 3)
          terms_list[[ind]] <- rbind( coef_save_p, coef_save_c )
          colnames(terms_list[[ind]]) <- colnames(data_list$x_q)

          # Turn on parameters if removal did not improve model fit
          if(lowest_val < val_sel){
            map[[params_loop[j]]] <- as.character( map[[params_loop[j]]] )
            map[[params_loop[j]]][i] <- i
            map[[params_loop[j]]] <- factor(map[[params_loop[j]]])
          }

          # Replace the lowest IC if parameter improves fit
          if(lowest_val > val_sel){
            lowest_val <- val_sel
          }

          print(paste0("Model ", ind, "converged"))
        }

        # Report
        # save(ind, file = paste0(c("biom", "bulinus")[as.numeric(bulinus)+1], "ind_",ind, "-j_",j,"-i_",i, ".RData"))

        rm(opt)
        ind <- ind + 1
      }
    }
  }

  # FIT THE FINAL MODEL
  data_list$x_q_hat_est <- as.numeric(as.character(map$beta_p))
  data_list$x_q_hat_est[which(is.na(data_list$x_q_hat_est))] <- as.numeric(as.character(map$beta_c))[which(is.na(data_list$x_q_hat_est))]
  obj = TMB::MakeADFun( data = expand_xq( data_list, map ), parameters = params, DLL = version, map = map_update( map , params , data_list ), random = random, silent = FALSE)
  opt <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 100000))

  # Save parameter names
  p_params_in_model <- paste(colnames(data_list$x_q)[which(!is.na(as.character( map[["beta_p"]] )))], collapse = ", ")
  c_params_in_model <- paste(colnames(data_list$x_q)[which(!is.na(as.character( map[["beta_c"]] )))], collapse = ", ")


  # If the model did not converge
  if(is.null(opt) == T | opt$convergence > 0){
    # Save values
    coef_save <- rep("Did not converge", length(terms_list[[ind - 1]]))
    terms_list[[ind]] <- coef_save

    summ_stats <- rbind(summ_stats, (c(ind, "Did not converge", "Did not converge", "Did not converge",  length(opt$par), p_params_in_model, c_params_in_model)))

    # Turn on parameters if removal led to non-convergence
    map[[params_loop[j]]] <- as.character( map[[params_loop[j]]] )
    map[[params_loop[j]]][i] <- i
    map[[params_loop[j]]] <- factor(map[[params_loop[j]]])
  }else if(is.null(opt) == F) {

    # # Refit
    # last_par = obj$env$parList(opt$par)
    # final_obj = TMB::MakeADFun(data = expand_xq( data_list, map ), parameters = last_par, DLL = version, map = map_update( map , params , data_list ), random = random, silent = TRUE)
    # final_opt = TMBhelper::Optimize( final_obj )

    # Save output
    BIC <- 2 * opt$value + log(nrow(data_list$x_q)) * length(opt$par)
    AIC <- 2 * opt$value + 2 * length(opt$par); if(is.null(AIC)){ AIC = NA }
    summ_stats <- rbind(summ_stats, (c(ind, AIC, BIC, opt$value,  length(opt$par), p_params_in_model, c_params_in_model)))

    # Save parameter estimates
    rep = TMB::sdreport(obj)
    coef_save_p <- round(rep$value[which(names(rep$value) %in% c("beta_p"))], 3)
    coef_save_c <- round(rep$value[which(names(rep$value) %in% c("beta_c"))], 3)
    terms_list[[ind]] <- rbind( coef_save_p, coef_save_c )
    colnames(terms_list[[ind]]) <- colnames(data_list$x_q)

  }


  # Return results
  results_list = list(summ_stats = summ_stats, mod_list = mod_list, terms_list = terms_list, Obj_list = Obj_list, final_obj = obj, final_opt = opt)
  return(results_list)
}
