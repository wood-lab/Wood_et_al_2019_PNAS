# Load functions
source("R/Support/1_read_dat_v3.R")
source("R/Support/2_fit_mesh_v2.R")
source("R/Support/3_build_data.R")
source("R/Support/4_build_params.R")
source("R/Support/5_build_map.R")
source("R/Support/bad_params.R")
source("R/run_variogram.R")
source("R/Support/mod_selection.R")
source("R/Support/aic_selection.R")
source("R/plot_coef.R")
library(glmmTMB)


# Step 0 -- Load data, create mesh, and assign to list for TMB
data <- read_dat(exclude = T,  fill_cont_mean = T, na_keep = T, quad = F, fill_cat_mode = T, scale_cont = T, bulinus = F) # read data in
density_data <- data$density_data
density_data$snail_is_biomph <- as.numeric(density_data$snail_is_biomph)

# Model 1
mod1 <- glmmTMB(snail_is_biomph ~ (1 | site) + 1,
                     family = nbinom1,
                     ziformula = ~ (1 | site) + 1, data = density_data)

pred <- predict(object = mod1, newdata = density_data, type = "response")
plot(density_data$snail_is_biomph, pred)
abline(0,1)


# Model 2
mod2 <- glmmTMB(snail_is_biomph ~ cont_covar_floating_veg_g_cleaned + (1 | site),
                family = nbinom2,
                ziformula = ~ cont_covar_floating_veg_g_cleaned + (1 | site), data = density_data)

pred <- predict(object = mod2, newdata = density_data, type = "response")
plot(density_data$snail_is_biomph, pred)
abline(0,1)


# Model 4
mod4 <- glmmTMB(snail_is_biomph ~ cont_covar_floating_veg_g_cleaned + (1 | site),
                family = genpois(link = "log"),
                ziformula = ~ cont_covar_floating_veg_g_cleaned + (1 | site), data = density_data)

pred <- predict(object = mod4, newdata = density_data, type = "response")
plot(density_data$snail_is_biomph, pred)
abline(0,1)


# Model 5
mod5 <- glmmTMB(snail_is_biomph ~ data_list$x_cq + (1 | site) + (1|field_mission),
                family = tweedie(link = "log"), data = density_data)

pred <- predict(object = mod5, newdata = density_data, type = "response")
plot(density_data$snail_is_biomph, pred)
abline(0,1)

