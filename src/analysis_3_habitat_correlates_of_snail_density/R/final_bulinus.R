model = 2
incl_disease = 0
bulinus = T
BIC_sel = T
fill_cont_mean = F
fill_cat_mode = F
outlier = T

# Load packages
require(INLA)
require(RandomFields)
require(RANN)
require(sp)
require(raster)
require(rgdal)
require(TMB)
require(readxl)
require(lattice)
require(TMBdebug)
require(gstat)
require(TMBhelper)
require(gtools)
require(plyr)

# Load functions
source("R/Support/1_read_dat_v3.R")
source("R/Support/2_fit_mesh_v2.R")
source("R/Support/3_build_data_v3.R")
source("R/Support/4_build_params_v3.R")
source("R/Support/5_build_map_v2.R")
source("R/Support/bad_params.R")
source("R/Support/run_variogram.R")
source("R/Support/two_stage_aic_selection_v2.R")
source("R/Support/plot_coef_v2.R")
source("R/Support/missing_cont_val_update_v2.R")
source("R/Support/7_build_x_q_array.R")


# Step 0 -- Load data, create mesh, and assign to list for TMB
data <- read_dat(exclude = T,  fill_cont_mean = fill_cont_mean, na_keep = T, quad = F, fill_cat_mode = fill_cat_mode, scale_cont = T, bulinus = bulinus, outlier = outlier) # read data in


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

# -----------------------------------------------------------------------------------------------------------------------------
# STEP 2
# -----------------------------------------------------------------------------------------------------------------------------
# Fit 1 -- Backward model selection
map <- build_map(model, data_list, incl_disease = incl_disease, params, space = 0, space_time = 0)
mod_list <- list()
terms_list <- list()
Obj_list <- list()

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
opt = tryCatch(TMBhelper::Optimize( obj ), error = function(e) NULL)


load("Report/2018-10-11/parameters_for_bulinus-model-2-BIC-no_mean_cont_fill-no_mode_cat_fill-_outlier_-TRUE.Rdata")
betas_p <- rep$value[which(names(rep$value) %in% c("beta_p"))]
betas_c <- rep$value[which(names(rep$value) %in% c("beta_c"))]
betas <- rbind(betas_p, betas_c)

for(i in 1:2){
  beta_params <- c("beta_p", "beta_c")
  map[[beta_params[i]]] <- as.character(map[[beta_params[i]]])
  map[[beta_params[i]]] <- rep(NA, length(map[[beta_params[i]]]))
  keep <- which(betas[i,] != 0)
  map[[beta_params[i]]][keep] <- 1:length(keep)
  map[[beta_params[i]]] <- factor(map[[beta_params[i]]])
}

# Turn on estimation of missing values
data_list$x_q_hat_est <- as.numeric(as.character(map$beta_p))
data_list$x_q_hat_est[which(is.na(data_list$x_q_hat_est))] <- as.numeric(as.character(map$beta_c))[which(is.na(data_list$x_q_hat_est))]

# Fit model
obj = TMB::MakeADFun( data = expand_xq( data_list, map ), parameters = params, DLL = version, map = map_update( map , params , data_list ), random = random, silent = TRUE)
opt = tryCatch(TMBhelper::Optimize( obj ), error = function(e) NULL)
# Step 3 -- make and compile template file

obj = TMB::MakeADFun( data = expand_xq( data_list, map ), parameters = obj$env$parList(opt$par), DLL = version, map = map_update( map , params , data_list ), random = random, silent = TRUE)
opt = tryCatch(TMBhelper::Optimize( obj ), error = function(e) NULL)
# Step 3 -- make and compile template file

rep_final <- TMB::sdreport(obj)

pred <- obj$report()

plot(x = data_list$c_q, y = pred$pred_cq)
abline(0,1)
legend("topleft", paste("AIC = ", round(sum(pred$jnll_comp[1:8,]), 0),"; ZAPL"), bty = "n")
