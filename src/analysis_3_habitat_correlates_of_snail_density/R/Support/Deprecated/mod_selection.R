# Function to do backward model selection using p-values, aic, or bic

calc_wald <- function( rep ){

  rep <- rep

  omega <- rep$value
  se <- rep$sd

  test <- abs(omega/se)

  p_value <- 1 - pchisq(test, 1)

  return(p_value)
}

get_stars = function(p) {
  stars = findInterval(p, c(0, 0.001, 0.01, 0.05, 0.1))
  codes = c("***" , "**","*", ".", " ")
  codes[stars]
  codes[which(is.na(codes))] <- " "
}

update_map <- function( rep , map, significance){

  p_values <- calc_wald( rep )

  beta_p_p_values <- p_values[which(names(p_values) == "beta_p")]
  beta_p_c_values <- p_values[which(names(p_values) == "beta_c")]
  beta_p_y_values <- p_values[which(names(p_values) == "beta_y")]

  beta_p_remove <- which(beta_p_p_values > significance | is.nan(beta_p_p_values))
  beta_c_remove <- which(beta_p_c_values > significance | is.nan(beta_p_c_values))
  beta_y_remove <- which(beta_p_y_values > significance | is.nan(beta_p_y_values))

  map$beta_p[beta_p_remove] <- NA
  map$beta_c[beta_c_remove] <- NA
  map$beta_y[beta_y_remove] <- NA

  map$beta_p <- factor(map$beta_p)
  map$beta_c <- factor(map$beta_c)
  map$beta_y <- factor(map$beta_y)

  return(map)
}

mod_fun <- function( data = data_list, parameters = params, DLL = version, map = map, random = random, significance = 0.05 ){
  library(plyr)

  mod_list <- list()
  terms_list <- list()
  Obj_list <- list()

  ind = 1
  remove_vec <- c()
  model.list = list()
  summ_stats = data.frame(matrix(NA, ncol = 5 , nrow = 1))
  colnames(summ_stats) = c("Model", "AIC","Log_Lik","N_Params","Removed")

  # Fit first model
  Obj = TMBdebug::MakeADFun( data = data_list, parameters = params, DLL = version, map = map, random = random)
  Opt = Optimize( Obj )

  E <- eigen(Opt$h)

  # Save outputs
  mod_list[[ind]] <- Opt
  Obj_list[[ind]] <- Obj

  p_values <- calc_wald( Opt$SD )
  p_save <- p_values[which(names(p_values) %in% c("beta_p", "beta_c"))]

  terms_list[[ind]] <- rbind(params$beta_p, params$beta_c)
  coef_save <- paste0(round(Opt$SD$value[which(names(Opt$SD$value) %in% c("beta_p", "beta_c"))], 3), get_stars(p_save))
  terms_list[[ind]] <- replace(terms_list[[ind]], values = coef_save)
  colnames(terms_list[[ind]]) <- colnames(data_list$x_cq)

  summ_stats <- rbind(summ_stats, (c("Base", Opt$AIC, Opt$objective,  length(Opt$par), NA)))

  # Run the next models
  map <- update_map(Opt$SD, map , significance)


  ind <- ind + 1
  # Loop model formulations until all terms are signficant
  while( sum(p_values > significance, na.rm = T) > 0 ){

    # Fit model
    Obj = TMBdebug::MakeADFun( data = data_list, parameters = params, DLL = version, map = map, random = random)
    Opt = Optimize( Obj )

    # Save outputs
    mod_list[[ind]] <- Opt
    Obj_list[[ind]] <- Obj

    p_values <- calc_wald( Opt$SD )
    p_save <- p_values[which(names(p_values) %in% c("beta_p", "beta_c"))]

    terms_list[[ind]] <- rbind(params$beta_p, params$beta_c)
    coef_save <- paste0(round(Opt$SD$value[which(names(Opt$SD$value) %in% c("beta_p", "beta_c"))], 3), get_stars(p_save))
    terms_list[[ind]] <- replace(terms_list[[ind]], values = coef_save)
    colnames(terms_list[[ind]]) <- colnames(data_list$x_cq)

    summ_stats <- rbind(summ_stats, (c("Base", Opt$AIC, Opt$objective,  length(Opt$par), NA)))

    # Run the next models
    map <- update_map(Opt$SD, map , significance)

    ind <- ind + 1
  }

  # PART II - SEARCH AROUND THE BEST MODEL
  mod_list[[ind]] <- "BREAK"
  Obj_list[[ind]] <- "BREAK"
  terms_list[[ind]] <- "BREAK"
  ind <- ind + 1
  mod_list[[ind]] <- "LOCAL MODEL SEARCH"
  Obj_list[[ind]] <- "LOCAL MODEL SEARCH"
  terms_list[[ind]] <- "LOCAL MODEL SEARCH"
  ind <- ind + 1

  summ_stats <- rbind(summ_stats, rep("BREAK", ncol(summ_stats)), rep("BREAK", ncol(summ_stats)) )


  # Get parameters that are not estimated
  last_params <- Opt$SD$value[which(names(Opt$SD$value) %in% c("beta_p", "beta_c"))]
  not_estimated_params <- which(last_params == 0)


  # Loop through non-estimated parameters to see if they improve fit
  for(i in not_estimated_params){

    # Change map to estimate parameter
    param_type <- names(not_estimated_params[which(not_estimated_params == i)])
    if(param_type == "beta_p"){param_ind = i}
    if(param_type == "beta_c"){param_ind = i - ncol(data_list$x_cq)}
    map[[param_type]] <- as.character(map[[param_type]])
    map[[param_type]][param_ind] = param_ind
    map[[param_type]] <- factor(map[[param_type]])

    # Fit model
    Obj = TMBdebug::MakeADFun( data = data_list, parameters = params, DLL = version, map = map, random = random)
    Opt = Optimize( Obj )

    # Save outputs
    mod_list[[ind]] <- Opt
    Obj_list[[ind]] <- Obj

    p_values <- calc_wald( Opt$SD )
    p_save <- p_values[which(names(p_values) %in% c("beta_p", "beta_c"))]

    terms_list[[ind]] <- rbind(params$beta_p, params$beta_c)
    coef_save <- paste0(round(Opt$SD$value[which(names(Opt$SD$value) %in% c("beta_p", "beta_c"))], 3), get_stars(p_save))
    terms_list[[ind]] <- replace(terms_list[[ind]], values = coef_save)
    colnames(terms_list[[ind]]) <- colnames(data_list$x_cq)

    summ_stats <- rbind(summ_stats, (c("Base", Opt$AIC, Opt$objective,  length(Opt$par), NA)))

    # Run the next models
    map <- update_map(Opt$SD, map , significance)

    ind <- ind + 1

  }

  # FIT THE FINAL MODEL
  Obj = TMBdebug::MakeADFun( data = data_list, parameters = params, DLL = version, map = map, random = random)
  Opt = Optimize( Obj )

  # Save outputs
  mod_list[[ind]] <- Opt
  Obj_list[[ind]] <- Obj

  p_values <- calc_wald( Opt$SD )
  p_save <- p_values[which(names(p_values) %in% c("beta_p", "beta_c"))]

  terms_list[[ind]] <- rbind(params$beta_p, params$beta_c)
  coef_save <- paste0(round(Opt$SD$value[which(names(Opt$SD$value) %in% c("beta_p", "beta_c"))], 3), get_stars(p_save))
  terms_list[[ind]] <- replace(terms_list[[ind]], values = coef_save)
  colnames(terms_list[[ind]]) <- colnames(data_list$x_cq)

  summ_stats <- rbind(summ_stats, (c("Base", Opt$AIC, Opt$objective,  length(Opt$par), NA)))


  # Return results
  results_list = list(summ_stats = summ_stats, mod_list = mod_list, terms_list = terms_list, Obj_list = Obj_list)
  return(results_list)
}
