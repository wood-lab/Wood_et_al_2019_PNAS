# Function to build parameter arguments for TMB

build_params <- function(model, data_list, incl_disease){

  params = list(

    # -- 2.0. Design matrix components
    "beta0" = rep(0, (2 + data_list$n_p)), # -- Regression intercepts
    "beta_p" = rep(0, ncol(data_list[["x_q"]])), # -- Regression coefficients for probability of presence
    "beta_c" = rep(0, ncol(data_list[["x_q"]])), # -- Regression coefficients for local density
    "beta_log" = rep(0, ncol(data_list[["x_q"]])),
    "beta_y" = matrix(0, nrow = data_list$n_p, ncol = dim(data_list[["x_i"]])[2]), # -- Regression coefficients for probability of disease
    "beta_c_hat" = rep(0, data_list$n_p),

    # -- 2.1. Random effects components
    "epsilon_mat" = matrix(0, nrow = 2 + data_list$n_p, ncol = data_list$n_s), # -- Random site effects
    "gamma_q" = rep(0, length(data_list[["c_q"]])), # -- Overdispersion parameters of log-normal Poisson

    # -- 2.2. Non-spatial variance components
    "log_sigmas" = rep(0, 6) # -- Variance components
  )


  if(incl_disease == 1){
    nparams <- 2 + data_list$n_p
  } else {
    nparams <- 2
  }

  params$gamma_q[which(data_list$c_q > 600)] <- 6

    params$x_q_cont_missing <- data_list$x_q[,1:(data_list$n_x_q_cont)]
    params$x_q_cont_missing <- replace(params$x_q_cont_missing, values = 0)


  # # -- 2.3. Spatial components
  # params[["log_tau_O"]] = rep(0, nparams)
  # params[["log_tau_E"]] = rep(0, nparams)
  # params[["log_kappa"]] = rep(0, nparams)
  #
  # # -- 2.4. Random spatio-temporal effects
  # params[["omega_s"]] = array(0, dim = c( max(data_list$dim_M_s), data_list$n_s, nparams))
  # params[["omega_st"]] = array(0, dim = c( max(data_list$dim_M_s), data_list$n_s, data_list$n_t, nparams))

  return(params)
}
