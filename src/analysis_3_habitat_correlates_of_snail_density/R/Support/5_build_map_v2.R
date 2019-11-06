# Function to build map arguments for TMB. Space and space-time are indicators wether: 0 - not fit; 1- only the presence part is fit, 2 - only the non-zero density is fit; 3- both presence and positive model compenents should be considered.

build_map <- function(model, data_list, incl_disease, params, space = 1, space_time = 0){

  map = list()

  map[["log_sigmas"]] = factor( c(1:6) )
  map[["beta0"]] = factor(c(1: (2 + data_list$n_p)))

  map[["beta_p"]] <- factor(1:length(params$beta_p))
  map[["beta_c"]] <- factor(1:length(params$beta_c))
  map[["beta_y"]] <- factor(1:length(params$beta_y))
  map[["beta_log"]] <- factor( rep(NA, data_list$n_x_q_cont), 1:(ncol(data_list$x_q) - data_list$n_x_q_cont))

  # -- 0. Not estimating probability of disease
  if(incl_disease == 0){ # Turn off random site effects for individual infection

    map[["log_sigmas"]][3:4] = factor( rep(NA, 2) )
    map[["beta0"]] = factor(c(1,2,rep(NA, data_list$n_p)))
    map[["beta_c_hat"]] = factor(rep(NA, data_list$n_p ))

    map[["epsilon_mat"]] = factor(matrix(c(1:(data_list$n_s*2), rep(NA, data_list$n_s * data_list$n_p)), nrow = 2 + data_list$n_p, ncol = data_list$n_s, byrow = T))

    map[["beta_y"]] = factor(matrix(NA, nrow = nrow(params$beta_y), ncol = ncol(params$beta_y)))

    # if(space == 0){
    #   map[["log_tau_O"]] = factor( rep(NA, 2) )
    # }
    #
    # if(space == 1){
    #   map[["log_tau_O"]] = factor(c(1, NA) )
    # }
    #
    # if(space == 2){
    #   map[["log_tau_O"]] = factor(c(NA, 1) )
    # }
    #
    # # Space time component
    # if(space_time == 0){ map[["log_tau_E"]] = factor( rep(NA, 2) ) }
    #
    # if(space_time == 1){ map[["log_tau_E"]] = factor( c(1, NA) ) }
    #
    # if(space_time == 2){ map[["log_tau_E"]] = factor( c(NA, 1) ) }
    #
    #
    # if(space_time == 0 & space == 0){ map[["log_kappa"]] = factor( rep(NA, 2) ) }
    #
    # if(space_time == 1 | space == 1){ map[["log_kappa"]] = factor( c(1, NA) ) }
    #
    # if(space_time == 2 | space == 2){ map[["log_kappa"]] = factor( c(NA, 1) ) }
  }


  # -- 1. Probability of infection estimated for ONE species
  if(data_list$n_p == 1 & incl_disease == 1){
    # Turn off random site SD and beta_y for 2nd parasite species
    map[["log_sigmas"]][4] = factor( NA )
    map[["beta_y"]] = factor(matrix( c(1:ncol(params$beta_y), rep(NA, ncol(params$beta_y))), nrow = nrow(params$beta_y), ncol = ncol(params$beta_y)))

    # if(space == F){
    #   map[["log_tau_O"]] = factor( rep(NA, 3) )
    # }
    #
    # if(space_time == F){
    #   map[["log_tau_E"]] = factor( rep(NA, 3) )
    # }
    #
    # if(space_time == F & space == F){
    #   map[["log_kappa"]] = factor( rep(NA, 3) )
    # }
  }



  if(model %in% c(1, 2, 3, 5, 6, 7)){ # Turn off CV for lognormal or p for tweedie
    map[["log_sigmas"]][6] = factor( NA )
  }

  if(model %in% c(0, 1, 3, 4, 5, 7, 8)){ # Turn off SD for overdispersion
    map[["log_sigmas"]][5] = factor( NA )
  }

  map[["log_sigmas"]] <- factor(as.character(map[["log_sigmas"]]))
  map[["beta0"]] = factor(as.character(map[["beta0"]]))

  # Spatial process
  if(incl_disease == 1){
    nparams <- 2 + data_list$n_p
  } else {
    nparams <- 2
  }

  # Missing continuous variables

    map$x_q_cont_missing <- params$x_q_cont_missing # Get continuous variables
    map$x_q_cont_missing <- replace(map$x_q_cont_missing, values = 1:length(map$x_q_cont_missing)) # Which are missing 1 = missing, 0 = not missing
    map$x_q_cont_missing[!is.na(data_list$x_q[,1:(data_list$n_x_q_cont)])] <- NA
    map$x_q_cont_missing <- factor(map$x_q_cont_missing)


  # Initialize map
  omega_s_count <- 1
  omega_st_count <- 1

  # map[["omega_s"]] = array( NA, dim = dim(params$omega_s))
  # map[["omega_st"]] = array( NA, dim = dim(params$omega_st))
  #
  # # Create data.frame of site and times where data are
  # data_check <- data.frame(s_q = data_list$s_q, t_q = data_list$t_q)
  # data_check <- as.data.frame(table(data_check$s_q, data_check$t_q)) # Count of quadrats in each field mission and site
  # colnames(data_check) <- c("site", "time", "freq")
  # data_check$site <- as.numeric(as.character(data_check$site)) + 1
  # data_check$time <- as.numeric(as.character(data_check$time)) + 1
  # data_check <- data_check[-which(data_check$freq == 0), ]
  #
  # # Fill space map
  # if(space %in% c(1,2,3)){
  #   space_params <- space
  #   if(space_params == 3){space_params = c(1,2)}
  #
  #   space_time_params <- space_time
  #   if(space_time_params == 3){space_params = c(1,2)}
  #
  #   # Map in parameters where data are
  #   for(p in space_params){ # Number of estimated paramters
  #     for(s in data_check$site){ # Number of sites
  #       for(k in 1:data_list$dim_M_s[s,1]){ # Number of vertices
  #
  #         map[["omega_s"]][k, s, p] <- omega_s_count
  #         omega_s_count <- omega_s_count + 1
  #
  #
  #       }
  #     }
  #   }
  # }
  #
  # # Fill space-time map
  # if(space_time %in% c(1,2,3)){
  #   for(ps in space_time_params){
  #     for(s in data_check$site){ # Number of sites
  #       for(k in 1:data_list$dim_M_s[s,1]){ # Number of vertices
  #         for(t in data_check$time){
  #           map[["omega_st"]][k, s, t, ps]  <- omega_st_count
  #           omega_st_count <- omega_st_count + 1
  #         }
  #       }
  #     }
  #   }
  # }
  #
  #
  # map[["omega_s"]] <- factor(map[["omega_s"]])
  # map[["omega_st"]] <- factor(map[["omega_st"]])

  return(map)
}
