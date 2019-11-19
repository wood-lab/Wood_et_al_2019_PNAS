model = 2; incl_disease = 0; bulinus = T; BIC_sel = TRUE; fill_cont_mean = T; fill_cat_mode = T; remove_outlier = FALSE; method = "nlminb"; silent = TRUE; savedir = "Report/June_Check"


# Load functions
library(glmmTMB)
source("R/Support/1_read_dat_v3.R")
source("R/Support/2_fit_mesh_v2.R")
source("R/Support/3_build_data_v3.R")
source("R/Support/4_build_params_v3.R")
source("R/Support/5_build_map_v2.R")
source("R/Support/bad_params.R")
source("R/run_variogram.R")
source("R/Support/aic_selection_v2.R")
source("R/plot_coef_v2.R")

# Creat the directory and file names
if(!is.null(savedir)){
  suppressWarnings(dir.create(file.path(savedir)))
}

file_name_base <- paste(c("biom", "bulinus")[as.numeric(bulinus)+1],"model",  model, c("AIC", "BIC")[as.numeric(BIC_sel)+1], c("no_mean_cont_fill", "mean_cont_fill")[as.numeric(fill_cont_mean)+1], c("no_mode_cat_fill", "mode_cat_fill")[as.numeric(fill_cat_mode)+1], "remove_outlier", remove_outlier, sep = "-")


# Step 0 -- Load data, create mesh, and assign to list for TMB
data <- read_dat(exclude = T,  fill_cont_mean = fill_cont_mean, na_keep = T, quad = F, fill_cat_mode = fill_cat_mode, scale_cont = T, bulinus = bulinus, remove_outlier = remove_outlier) # read data in
density_data <- data$density_data
density_data$snail_is_trunc_globo <- as.numeric(density_data$snail_is_trunc_globo)

# Model 1 - NBINOM 1
mod1 <- glmmTMB(snail_is_trunc_globo ~ cont_covar_floating_veg_g_cleaned + cat_covar_ludwigia + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_lake_river +  (1 | site),
                family = nbinom1, data = density_data)

mod1a <- glmmTMB(snail_is_trunc_globo ~ cont_covar_floating_veg_g_cleaned + cat_covar_ludwigia + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_lake_river +  (1 | site),
                 ziformula = ~ cont_covar_floating_veg_g_cleaned + cat_covar_ludwigia + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_lake_river +  (1 | site),
                 family = nbinom1, data = density_data)

# Model 2 - NBIMON 2
mod2 <- glmmTMB(snail_is_trunc_globo ~ cont_covar_floating_veg_g_cleaned + cat_covar_ludwigia + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_lake_river +  (1 | site),
                family = nbinom2, data = density_data)

mod2a <- glmmTMB(snail_is_trunc_globo ~ cont_covar_floating_veg_g_cleaned + cat_covar_ludwigia + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_lake_river +  (1 | site),
                 ziformula = ~ cont_covar_floating_veg_g_cleaned + cat_covar_ludwigia + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_lake_river +  (1 | site),
                 family = nbinom2, data = density_data)

# Model 3 - Poisson
mod3 <- glmmTMB(snail_is_trunc_globo ~ cont_covar_floating_veg_g_cleaned + cat_covar_ludwigia + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_lake_river +  (1 | site),
                family = poisson(link = "log"), data = density_data)

mod3a <- glmmTMB(snail_is_trunc_globo ~ cont_covar_floating_veg_g_cleaned + cat_covar_ludwigia + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_lake_river +  (1 | site),
                 ziformula = ~ cont_covar_floating_veg_g_cleaned + cat_covar_ludwigia + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_lake_river +  (1 | site),
                 family = poisson(link = "log"), data = density_data)

# Model 4 - Lognormal
mod4 <- glmmTMB(snail_is_trunc_globo ~ cont_covar_floating_veg_g_cleaned + cat_covar_ludwigia + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_lake_river +  (1 | site), data = density_data)

mod4a <- glmmTMB(snail_is_trunc_globo ~ cont_covar_floating_veg_g_cleaned + cat_covar_ludwigia + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_lake_river +  (1 | site),
                 ziformula = ~ cont_covar_floating_veg_g_cleaned + cat_covar_ludwigia + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_lake_river +  (1 | site), data = density_data)


# Model 4 - tweedie
mod5 <- glmmTMB(snail_is_trunc_globo ~ cont_covar_floating_veg_g_cleaned + cat_covar_ludwigia + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_lake_river +  (1 | site), family = tweedie(link = "log"),  data = density_data)


#######################################################################################
# delta-PL

# Load functions
source("R/Support/1_read_dat_v3.R")
source("R/Support/2_fit_mesh_v2.R")
source("R/Support/3_build_data_v3.R")
source("R/Support/4_build_params_v3.R")
source("R/Support/5_build_map_v2.R")
source("R/Support/bad_params.R")
source("R/run_variogram.R")
source("R/Support/aic_selection_v2.R")
source("R/plot_coef_v2.R")

# Creat the directory and file names
if(!is.null(savedir)){
  suppressWarnings(dir.create(file.path(savedir)))
}

file_name_base <- paste(c("biom", "bulinus")[as.numeric(bulinus)+1],"model",  model, c("AIC", "BIC")[as.numeric(BIC_sel)+1], c("no_mean_cont_fill", "mean_cont_fill")[as.numeric(fill_cont_mean)+1], c("no_mode_cat_fill", "mode_cat_fill")[as.numeric(fill_cat_mode)+1], "remove_outlier", remove_outlier, sep = "-")


# Step 0 -- Load data, create mesh, and assign to list for TMB
data <- read_dat(exclude = T,  fill_cont_mean = fill_cont_mean, na_keep = T, quad = F, fill_cat_mode = fill_cat_mode, scale_cont = T, bulinus = bulinus, remove_outlier = remove_outlier) # read data in


# data <- fit_mesh(data, cluster = T, n_knots = 15) # Fit mesh
data_list <- build_data(data, model, debug = 0, bulinus = bulinus, parasite = 0) # Create for TMB


# Step 1 -- Set initial values for parameters
params <- build_params(model, data_list, incl_disease)
random = c("epsilon_mat", "omega_mat", "gamma_q", "x_q_cont_missing")


# Step 3 -- make and compile template file
setwd("src")
library(TMB)
# library(TMBdebug)
library(TMBhelper)
version = "temporal_model_spde_v4"
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



###############################################
# Load packages
library(plyr)
library(gtools)
library(TMB)
library(TMBhelper)
library(foreach)
library(doParallel)

source("R/Support/8_missing_cont_val_update_v2.R")
source("R/Support/7_build_x_q_array.R")
source("R/Optimize_optim.R")

###############################################
# Fit intercept model: No betas estimated
map$beta_c <- factor(c(NA, NA, NA, NA, NA, NA, NA, 1, NA, NA, NA, NA, NA, 2,NA, NA, 3,  4, NA, NA, NA, NA, 5))
map$beta_p <- factor(c(NA, NA, NA, NA, NA, NA, NA, 1, NA, NA, NA, NA, NA, 2,NA, NA, 3,  4, NA, NA, NA, NA, 5))
map$beta_log <- factor(rep(NA, length(map$beta_log)))
map$x_q_cont_missing <- factor(rep(NA, length(map$x_q_cont_missing)))

# Set up data
data_list$x_q_hat_est <- as.numeric(as.character(map$beta_p))
new_map <- map_update( map , params , data_list )
new_dat <- expand_xq( data_list, map )
params$x_q_cont_missing <- rep(0, sum(is.na(new_dat$x_q[which(new_dat$fit_ll == 1),1:(new_dat$n_x_q_cont)])))


# Fit base model
obj = TMB::MakeADFun( data = new_dat, parameters = params, DLL = version, map = new_map, random = random, silent = TRUE)
opt = suppressWarnings( Optimize( obj , n = log(nrow(new_dat$x_q)), loopnum = 3, getsd = FALSE))

pred_vec <- obj$report(obj$env$last.par.best)$pred_cq
jnll_comp <- obj$report(obj$env$last.par.best)$jnll_comp

bic_temp <-2 * sum(jnll_comp[1:7,]) + log(nrow(data_list$x_q)) * length(opt$par)


#######################################################################################
# Plot it
mod_list <- list(mod1, mod1a, mod2, mod2a, mod3, mod3a, mod4, mod4a, mod5)
mod_names <- c("Nbinom1", "zi-Nbinom1", "Nbinom2", "zi-Nbinom2", "Poisson", "zi-Poisson", "Log-normal", "zi-Lognormal", "tweedie")

par(mfrow = c(2,5))

for(i in 1:length(mod_list)){
  plot(density_data$snail_is_trunc_globo, predict(mod_list[[i]], type = "response"), ylab = "Predicted", xlab = "Observed", ylim = c(0, 150), xlim = c(0, 150))
  abline(0, 1)
  legend("topleft", c(mod_names[i], paste0("BIC = ", round(BIC(mod_list[[i]]), 0))), bty = "n")
}

plot(density_data$snail_is_trunc_globo, pred_vec, ylab = "Predicted", xlab = "Observed", ylim = c(0, 150), xlim = c(0, 150))
abline(0, 1)
legend("topleft", c("delta-PL", paste0("BIC = ", round(bic_temp, 0))), bty = "n")





