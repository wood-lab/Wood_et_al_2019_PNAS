# Function to fit simulated data from the wood-adams delta model

fit_sim<- function(model = 6, incl_disease = 1){

  # Step 0 -- Load simulator and simulate data. See file for specifications
  setwd("R")
  source("mod_sim.R")
  data <- sim_dat()
  setwd("../")

  # Assign data to list
  data_list <- list(
    # Settings
    model = model,
    incl_disease = 0,
    debug = 1,

    # -- 0.1 Local density data
    n_q = length(data$sim_quad$s_q),
    n_s = length(unique(data$sim_quad$s_q)),
    c_q = data$sim_quad$c_q,
    v_q = rep(1, length(data$sim_quad$s_q)),
    x_pq = as.matrix(rep(1, length(data$sim_quad$s_q))), # No covariates
    x_cq = as.matrix(rep(1, length(data$sim_quad$s_q))), # No Covariates

    # 0.2. -- Individual snail data
    n_i = dim(as.matrix(data$sim_y$y_i))[1],
    n_p = dim(as.matrix(data$sim_y$y_i))[2], # Number of parasites
    y_i = as.matrix(data$sim_y$y_i),
    q_i = data$sim_y$id_i - 1,
    x_i = as.matrix(rep(1, length(data$sim_y$q_i)))
  )

  if( min(data$sim_quad$s_q) == 1){ data_list[["s_q"]] = data$sim_quad$s_q - 1}
  if( min(data$sim_quad$s_q) == 0){ data_list[["s_q"]] = data$sim_quad$s_q}

  if( min(data$sim_y$s_i) == 1){ data_list[["s_i"]] = data$sim_y$s_i - 1}
  if( min(data$sim_y$s_i) == 0){ data_list[["s_i"]] = data$sim_y$s_i}

  data_list[["predTF_q"]] = rep(0, length(data_list$c_q))
  data_list[["predTF_i"]] = rep(0, nrow(data_list$y_i))
  data_list[["model"]] = model
  data_list[["incl_disease"]] = incl_disease


  # Step 1 -- Set initial values for parameters
  params = list(
    "beta0" = rep(1, (2 + data_list$n_p)), # -- Regression intercepts
    "beta_p" = rep(0, ncol(data_list[["x_pq"]])), # -- Regression coefficients for probability of presence
    "beta_c" = rep(0, ncol(data_list[["x_cq"]])), # -- Regression coefficients for local density
    "beta_y" = matrix(0, nrow = data_list$n_p, ncol = dim(data_list[["x_i"]])[2]), # -- Regression coefficients for probability of disease
    "beta_c_hat" = rep(0, data_list$n_p),
    "epsilon_mat" = matrix(0, nrow = 2 + data_list$n_p, ncol = data_list$n_s), # -- Random site effects
    "gamma_q" = rep(0, length(data_list[["c_q"]])), # -- Overdispersion parameters of log-normal Poisson
    "log_sigmas" = rep(0, 6) # -- Variance components
  )

  random = c("epsilon_mat", "gamma_q")


  # Step 2 -- Build map function.
  # -- 2.0. Turn off beta parameters
  map = list()
  map[["beta_p"]] = factor( NA )
  map[["beta_c"]] = factor( NA )
  map[["beta_y"]] = factor( NA )

  # -- 2.1. Make sure variance and intercepts is correct for model structure
  map[["log_sigmas"]] = factor( c(1:6) )
  map[["beta0"]] = factor(c(1: (2 + data_list$n_p)))

  if(incl_disease == 0){ # Turn off random site effects for individual infection
    map[["log_sigmas"]][4:5] = factor( rep(NA, 2) )
    map[["beta0"]] = factor(c(1,2,rep(NA, data_list$n_p)))
    map[["beta_c_hat"]] = factor(c(NA,NA))
    map[["epsilon_mat"]] = factor(matrix(c(1:(data_list$n_s*2), rep(NA, data_list$n_s * data_list$n_p)), nrow = 2 + data_list$n_p, ncol = data_list$n_s, byrow = T))
  }

  if(data_list$n_p == 1){ # Turn off random site effects for 2nd parasite species
    map[["log_sigmas"]][5] = factor( NA )
  }

  if(model %in% c(1, 2, 3, 5, 6, 7)){ # Turn off CV for lognormal
    map[["log_sigmas"]][3] = factor( NA )
  }

  if(model %in% c(0, 1, 3, 4, 5, 7)){ # Turn off SD for overdispersion
    map[["log_sigmas"]][6] = factor( NA )
  }

  map[["log_sigmas"]] <- factor(as.character(map[["log_sigmas"]]))
  map[["beta0"]] = factor(as.character(map[["beta0"]]))


  # Step 3 -- make and compile template file
  setwd("src")
  library(TMB)
  library(TMBdebug)
  compile( "spatio_temporal_model_test.cpp" )
  dyn.load( dynlib("spatio_temporal_model_test") )
  setwd("../")

  # Fit without RE and infection
  Obj = TMBdebug::MakeADFun( data = data_list, parameters = params, DLL="spatio_temporal_model_test", map = map)
  Obj$fn( Obj$par )
  Obj$gr( Obj$par )
  Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
  par = Obj$env$parList(Opt$par)
  rm(Opt)

  # Fit without infection
  Obj = TMBdebug::MakeADFun( data = data_list, parameters = par, DLL="spatio_temporal_model_test", map = map, random = random)
  Opt = TMBhelper::Optimize( obj=Obj, newtonsteps=1, getsd=FALSE, control = list(iter.max = 1e5))
  par = Obj$env$parList(Opt$par)
  rm(Opt)

  # Fit with infection
  if(incl_disease == 1){
    data_list$incl_disease = 1
    Obj = TMBdebug::MakeADFun( data = data_list, parameters = par, DLL="spatio_temporal_model_test", map = map, random = random)
    Opt = TMBhelper::Optimize( obj=Obj, newtonsteps=1, getsd=FALSE, control = list(iter.max = 1e5))
    par = Obj$env$parList(Opt$par)
  }
}
