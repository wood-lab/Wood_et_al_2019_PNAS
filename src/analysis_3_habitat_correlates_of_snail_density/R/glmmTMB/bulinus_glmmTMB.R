# Load functions
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
library(glmmTMB)


# Step 0 -- Load data, create mesh, and assign to list for TMB
data <- read_dat(exclude = T,  fill_cont_mean = T, na_keep = T, quad = F, fill_cat_mode = T, scale_cont = T, bulinus = F) # read data in
density_data <- data$density_data
density_data$snail_is_trunc_globo <- as.numeric(density_data$snail_is_trunc_globo)

density_data$cat_covar_ceratophyllum

par(mfrow = c(2,3))
# ZINB1
mod1 <- glmmTMB(snail_is_trunc_globo ~ (1 | site) + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_ludwigia + cont_covar_floating_veg_g_cleaned + cat_covar_lake_river,
                family = nbinom1,
                ziformula = ~ (1 | site) + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_ludwigia + cont_covar_floating_veg_g_cleaned + cat_covar_lake_river, data = density_data)

glmm_pred <- predict(object = mod1, newdata = density_data, type = "response")
plot(density_data$snail_is_trunc_globo, glmm_pred, xlim = c(0, 1000), ylim = c(0, 1000))
abline(0,1)
legend("topleft", paste("AIC = ", round(AIC(mod1) , 0),"; ZINB1"), bty = "n")

# ZINB2
mod2 <- glmmTMB(snail_is_trunc_globo ~ (1 | site) + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_ludwigia + cont_covar_floating_veg_g_cleaned + cat_covar_lake_river,
                family = nbinom2,
                ziformula = ~ (1 | site) + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_ludwigia + cont_covar_floating_veg_g_cleaned + cat_covar_lake_river, data = density_data)

glmm_pred <- predict(object = mod2, newdata = density_data, type = "response")
plot(density_data$snail_is_trunc_globo, glmm_pred, xlim = c(0, 1000), ylim = c(0, 1000))
abline(0,1)
legend("topleft", paste("AIC = ", round(AIC(mod2), 0),"; ZINB2"), bty = "n")

# ZIP
mod3 <- glmmTMB(snail_is_trunc_globo ~ (1 | site) + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_ludwigia + cont_covar_floating_veg_g_cleaned + cat_covar_lake_river,
                family = poisson(link = "log"),
                ziformula = ~ (1 | site) + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_ludwigia + cont_covar_floating_veg_g_cleaned + cat_covar_lake_river, data = density_data)

glmm_pred <- predict(object = mod3, newdata = density_data, type = "response")
plot(density_data$snail_is_trunc_globo, glmm_pred, xlim = c(0, 1000), ylim = c(0, 1000))
abline(0,1)
legend("topleft", paste("AIC = ", round(AIC(mod3) , 0),"; ZIP"), bty = "n")

# NB1
mod4 <- glmmTMB(snail_is_trunc_globo ~ (1 | site) + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_ludwigia + cont_covar_floating_veg_g_cleaned + cat_covar_lake_river,
                family = poisson(link = "log"), data = density_data)

glmm_pred <- predict(object = mod4, newdata = density_data, type = "response")
plot(density_data$snail_is_trunc_globo, glmm_pred, xlim = c(0, 1000), ylim = c(0, 1000))
abline(0,1)
legend("topleft", paste("AIC = ", round(AIC(mod4), 0),"; NB1"), bty = "n")

# NB2
mod5 <- glmmTMB(snail_is_trunc_globo ~ (1 | site) + cat_covar_ceratophyllum + cat_covar_potamogeton + cat_covar_ludwigia + cont_covar_floating_veg_g_cleaned + cat_covar_lake_river,
                family = nbinom1, data = density_data)

glmm_pred <- predict(object = mod5, newdata = density_data, type = "response")
plot(density_data$snail_is_trunc_globo, glmm_pred, xlim = c(0, 1000), ylim = c(0, 1000))
abline(0,1)
legend("topleft", paste("AIC = ", round(AIC(mod5), 0),"; NB2"), bty = "n")

# ZIPL
plot(x = data_list$c_q, y = pred$pred_cq, xlim = c(0, 1000), ylim = c(0, 1000))
abline(0,1)
legend("topleft", paste("AIC = ", round(sum(pred$jnll_comp[1:8,]), 0),"; ZAPL"), bty = "n")

