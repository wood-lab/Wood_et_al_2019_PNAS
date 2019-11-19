# Function to build parameter arguments for TMB

build_params <- function(model, data_list, incl_disease){

  params = list(

    # -- 2.0. Design matrix components
    "beta0" = rep(0, 2), # -- Regression intercepts
    "beta_p" = rep(0, ncol(data_list[["x_q"]])), # -- Regression coefficients for probability of presence
    "beta_c" = rep(0, ncol(data_list[["x_q"]])), # -- Regression coefficients for local density
    "beta_log" = rep(0, ncol(data_list[["x_q"]])),

    # -- 2.1. Random effects components
    "epsilon_mat" = matrix(0, nrow = 2, ncol = data_list$n_s), # -- Random site effects
    "omega_mat" = matrix(0, nrow = 2, ncol = data_list$n_v), # -- Random village effects

    "gamma_q" = rep(1, length(data_list[["c_q"]])), # -- Overdispersion parameters of log-normal Poisson

    # -- 2.2. Non-spatial variance components
    "log_sigmas" = rep(0, 6) # -- Variance components
  )

  nparams <- 2

  # Outliers, start high for convergence
  outliers <- which(data_list$c_q > 600)

  # Set to 6 if fitting, set to 0 if not
  if(length(outliers) > 0){
    params$gamma_q[outliers] <- sapply(data_list$fit_ll[outliers], function(x) ifelse( x == 1, 6, 0))
  }

  params$x_q_cont_missing <- rep(0, sum(is.na(data_list$x_q[which(data_list$fit_ll == 1),1:(data_list$n_x_q_cont)])))


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
