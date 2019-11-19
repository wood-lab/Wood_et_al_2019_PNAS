rm(list = ls())
bulinus = T
model = 2
incl_disease = 0

# Load functions
source("R/Support/1_read_dat_v3.R")
source("R/Support/2_fit_mesh_v2.R")
source("R/Support/3_build_data_v3.R")
source("R/Support/4_build_params_v3.R")
source("R/Support/5_build_map_v2.R")
source("R/Support/bad_params.R")
source("R/run_variogram.R")
source("R/Support/aic_selection_v2.R")
source("R/plot_coef.R")


# Step 0 -- Load data, create mesh, and assign to list for TMB
data <- read_dat(exclude = T,  fill_cont_mean = F, na_keep = T, quad = F, fill_cat_mode = T, scale_cont = T, bulinus = bulinus) # read data in
#load("without_hmm.RData")



# data <- fit_mesh(data, cluster = T, n_knots = 15) # Fit mesh
data_list <- build_data(data, model, debug = 0, bulinus = bulinus, parasite = 0) # Create for TMB


# Step 1 -- Set initial values for parameters
params <- build_params(model, data_list, incl_disease)
random = c("epsilon_mat", "gamma_q", "x_q_cont_missing")#, "omega_s", "omega_st")


# Step 3 -- make and compile template file
setwd("src")
library(TMB)
library(TMBdebug)
library(TMBhelper)
version = "temporal_model_spde_v3"
if(data_list$debug == 1){
  #dyn.unload(version)
  #file.remove(paste0(version,".dll"))
  #file.remove(paste0(version,".o"))
}
TMB::compile( paste0(version,".cpp") )
dyn.load(dynlib(version))
setwd("../")

map <- build_map(model, data_list, incl_disease = incl_disease, params, space = 0, space_time = 0)

library(plyr)
source("R/Support/missing_cont_val_update_v2.R")
source("R/Support/7_build_x_q_array.R")

mod_list <- list()
terms_list <- list()
Obj_list <- list()

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
data_list <- expand_xq( data_list, map )
Obj = TMB::MakeADFun( data = data_list, parameters = params, DLL = version, map = map_update( map , params , data_list ), random = random)
Opt = Optimize( Obj )
rep2 <- Obj$report()
rep2$jnll_comp[,1007] # Should equal 231.40096
rep$jnll_comp[,1007]
#head(rep$linpred_pq_array[,,1])
#sum(data_list$x_q_array[,,1])
