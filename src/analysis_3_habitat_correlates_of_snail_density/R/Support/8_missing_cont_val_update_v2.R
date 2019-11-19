# This function updates the map to not imput values for missing data when the

map_update <- function( map , params , data_list ){

  beta_p <- as.character(map$beta_p)
  beta_c <- as.character(map$beta_c)

  # turn off logistic categorical variable parameters
  map[["beta_log"]] <- factor( c(rep(NA, data_list$n_x_q_cont), 1:(ncol(data_list$x_q) - data_list$n_x_q_cont)))
  # map[["beta_log"]][which(is.na(beta_p))] <- NA
  map[["beta_log"]] <- factor(map[["beta_log"]])

  # Overdispersion
  map$gamma_q <- 1:length(params$gamma_q)
  if(model %in% c(0, 1, 3, 4, 5, 7, 8)){ # Turn off SD for overdispersion
    map$gamma_q <- replace(map$gamma_q, values = NA)
  }
  map$gamma_q[which(data_list$fit_ll == 0)] <- NA
  map$gamma_q <- factor(map$gamma_q)
  map$x_q_cont_missing <- NULL

  return(map)
}
