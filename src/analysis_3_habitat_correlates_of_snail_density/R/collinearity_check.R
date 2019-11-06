model = 2; incl_disease = 0; bulinus = T; BIC_sel = T; fill_cont_mean = T; fill_cat_mode = T; outlier = T; two_stage = T


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
library(corrplot)

# Load functions
source("src/analysis_3_habitat_correlates_of_snail_density/R/Support/1_read_dat_v3.R")
source("src/analysis_3_habitat_correlates_of_snail_density/R/Support/2_fit_mesh_v2.R")
source("src/analysis_3_habitat_correlates_of_snail_density/R/Support/3_build_data_v3.R")
source("src/analysis_3_habitat_correlates_of_snail_density/R/Support/4_build_params_v3.R")
source("src/analysis_3_habitat_correlates_of_snail_density/R/Support/5_build_map_v2.R")
source("src/analysis_3_habitat_correlates_of_snail_density/R/Support/bad_params.R")
source("src/analysis_3_habitat_correlates_of_snail_density/R/Support/run_variogram.R")
source("src/analysis_3_habitat_correlates_of_snail_density/R/Support/two_stage_aic_selection_v2_optimx.R")
source("src/analysis_3_habitat_correlates_of_snail_density/R/Support/aic_selection_v3.R")
source("src/analysis_3_habitat_correlates_of_snail_density/R/Support/plot_coef_v2.R")


# Step 0 -- Load data, create mesh, and assign to list for TMB
data <- read_dat(exclude = T,  fill_cont_mean = fill_cont_mean, na_keep = T, quad = F, fill_cat_mode = fill_cat_mode, scale_cont = T, bulinus = bulinus, outlier = outlier) # read data in

data_list <- build_data(data, model, debug = 0, bulinus = bulinus, parasite = 0) # Create for TMB

res <- cor(data_list$x_q, use = "pairwise.complete.obs")
corrplot(res)


design <- data_list$x_q

e <- eigen(t(design) %*% design)
e$val

condition_check <-  sqrt(e$val[1] / e$val)

library(faraway)
collinear <- cbind(condition_check, vif(design))
write.csv(vif(design), file = "VIF.csv")


x <- runif(100, 0, 1)

vif(cbind(x,x))
