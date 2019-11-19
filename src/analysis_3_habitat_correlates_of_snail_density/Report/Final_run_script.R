---
  title: "Analysis 3: Identifying habitat-related predictors of snail abundance"
author: "Grant Adams"
date: "updated June 2019"
---

# Load dependencies
library( TMB )
library( TMBhelper )
library( lattice )
library( readxl )
library( RandomFields )
library( RANN )
library( sp )
library( rgdal )
library( gtools )
library( sp )
library( gstat )
library( FSA )
library( Hmisc )
library( plyr )

# Final model runs
# Load function
source("R/fit_mod.R")

# Run with outlier
fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = TRUE, fill_cont_mean = F, fill_cat_mode = F, remove_outlier = FALSE, method = "nlminb", silent = TRUE, savedir = "Report/Final_run_outlier") # Bulinus delta-pl
fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = TRUE, fill_cont_mean = F, fill_cat_mode = T, remove_outlier = FALSE, method = "nlminb", silent = TRUE, savedir = "Report/Final_run_outlier") # Bulinus delta-pl
fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = TRUE, fill_cont_mean = T, fill_cat_mode = F, remove_outlier = FALSE, method = "nlminb", silent = TRUE, savedir = "Report/Final_run_outlier") # Bulinus delta-pl
fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = TRUE, fill_cont_mean = T, fill_cat_mode = T, remove_outlier = FALSE, method = "nlminb", silent = TRUE, savedir = "Report/Final_run_outlier") # Bulinus delta-pl


# Run without outlier
fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = TRUE, fill_cont_mean = F, fill_cat_mode = F, remove_outlier = TRUE, method = "nlminb", silent = TRUE, savedir = "Report/Final_run_no_outlier") # Bulinus delta-pl
fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = TRUE, fill_cont_mean = F, fill_cat_mode = T, remove_outlier = TRUE, method = "nlminb", silent = TRUE, savedir = "Report/Final_run_no_outlier") # Bulinus delta-pl
fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = TRUE, fill_cont_mean = T, fill_cat_mode = F, remove_outlier = TRUE, method = "nlminb", silent = TRUE, savedir = "Report/Final_run_no_outlier") # Bulinus delta-pl
fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = TRUE, fill_cont_mean = T, fill_cat_mode = T, remove_outlier = TRUE, method = "nlminb", silent = TRUE, savedir = "Report/Final_run_no_outlier") # Bulinus delta-pl
