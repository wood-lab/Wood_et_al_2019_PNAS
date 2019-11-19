# Function to build map arguments for TMB. Space and space-time are indicators wether: 0 - not fit; 1- only the presence part is fit, 2 - only the non-zero density is fit; 3- both presence and positive model compenents should be considered.

build_map <- function(model, data_list, incl_disease, params, space = 1, space_time = 0, village_on = FALSE){

  map = list()

  map[["log_sigmas"]] = factor( c(1:6) )
  map[["beta0"]] = factor(c(1:2))

  map[["beta_p"]] <- factor(1:length(params$beta_p))
  map[["beta_c"]] <- factor(1:length(params$beta_c))
  map[["beta_log"]] <- factor( rep(NA, data_list$n_x_q_cont), 1:(ncol(data_list$x_q) - data_list$n_x_q_cont))
  map$gamma_q <- 1:length(params$gamma_q)
map$x_q_cont_missing <- NULL

  if(model %in% c(1, 2, 3, 5, 6, 7)){ # Turn off CV for lognormal or p for tweedie
    map[["log_sigmas"]][6] = factor( NA )
  }

  if(model %in% c(0, 1, 3, 4, 5, 7, 8)){ # Turn off SD for overdispersion
    map[["log_sigmas"]][5] = factor( NA )
    map$gamma_q <- replace(map$gamma_q, values = NA)
  }
  map[["beta0"]] = factor(as.character(map[["beta0"]]))

  # Spatial process

  nparams <- 2

  # Overdispersion turn off if not estimating
  map$gamma_q[which(data_list$fit_ll == 0)] <- NA
  map$gamma_q <- factor(map$gamma_q)

  # Turn off random site effects where villages have 1 site
  site_dat <- data.frame(site = data_list$s_q, village = data_list$v_q)
  unique_sv <- unique(site_dat)
  unique_v <- table(unique_sv$village)
  single_v <- as.numeric(names(unique_v[which(unique_v == 1)]))
  site_turn_off <- unique_sv$site[which(unique_sv$village %in% single_v)] + 1 # add 1 because R indexing

  if(village_on){
    map$epsilon_mat <- params$epsilon_mat
    map$epsilon_mat <- replace(map$epsilon_mat, values = 1:(length(map$epsilon_mat)))
    map$epsilon_mat[,site_turn_off] <- NA
    map$epsilon_mat <- factor(map$epsilon_mat)
  }

  if(!village_on){
    map$omega_mat <- params$omega_mat
    map$omega_mat <- replace(map$omega_mat, values = NA)
    map$omega_mat <- factor(map$omega_mat)

    map$log_sigmas[3:4] <-  NA
  }



  map[["log_sigmas"]] <- factor(as.character(map[["log_sigmas"]]))

  return(map)
}
