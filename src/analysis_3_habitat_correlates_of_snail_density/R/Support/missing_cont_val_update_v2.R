# This function updates the map to not imput values for missing data when the

map_update <- function( map , params , data_list ){

  beta_p <- as.character(map$beta_p)
  beta_c <- as.character(map$beta_c)

  # turn off logistic categorical variable parameters
  map[["beta_log"]] <- factor( c(rep(NA, data_list$n_x_q_cont), 1:(ncol(data_list$x_q) - data_list$n_x_q_cont)))
  # map[["beta_log"]][which(is.na(beta_p))] <- NA
  map[["beta_log"]] <- factor(map[["beta_log"]])

  map$x_q_cont_missing <- params$x_q_cont_missing # Get continuous variables
  map$x_q_cont_missing <- replace(map$x_q_cont_missing, values = 1:length(map$x_q_cont_missing)) # Which are missing 1 = missing, 0 = not missing
  map$x_q_cont_missing[!is.na(data_list$x_q[,1:(data_list$n_x_q_cont)])] <- NA
  map$x_q_cont_missing[,which(1:ncol(map$x_q_cont_missing) %in% which(is.na(beta_p)))] <- NA # Don't estimate data for the parameters that are turned off

  map$x_q_cont_missing <- factor(map$x_q_cont_missing)

  return(map)
}
