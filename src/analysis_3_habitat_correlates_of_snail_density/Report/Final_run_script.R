
# Install dependencies
devtools::install_github("kaskr/TMB_contrib_R/TMBdebug")

# Final model runs
# Load functions
source("R/fit_mod.R")

# With outlier
fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = TRUE, fill_cont_mean = F, fill_cat_mode = F, remove_outlier = FALSE, method = "nlminb", silent = TRUE, savedir = "Report/Final_run") # Bulinus delta-pl
fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = TRUE, fill_cont_mean = F, fill_cat_mode = T, remove_outlier = FALSE, method = "nlminb", silent = TRUE, savedir = "Report/Final_run") # Bulinus delta-pl
fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = TRUE, fill_cont_mean = T, fill_cat_mode = F, remove_outlier = FALSE, method = "nlminb", silent = TRUE, savedir = "Report/Final_run") # Bulinus delta-pl
fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = TRUE, fill_cont_mean = T, fill_cat_mode = T, remove_outlier = FALSE, method = "nlminb", silent = TRUE, savedir = "Report/Final_run") # Bulinus delta-pl



# Without outlier
fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = TRUE, fill_cont_mean = F, fill_cat_mode = F, remove_outlier = TRUE, method = "nlminb", silent = TRUE, savedir = "Report/Final_run_no_outlier") # Bulinus delta-pl
fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = TRUE, fill_cont_mean = F, fill_cat_mode = T, remove_outlier = TRUE, method = "nlminb", silent = TRUE, savedir = "Report/Final_run_no_outlier") # Bulinus delta-pl
fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = TRUE, fill_cont_mean = T, fill_cat_mode = F, remove_outlier = TRUE, method = "nlminb", silent = TRUE, savedir = "Report/Final_run_no_outlier") # Bulinus delta-pl
fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = TRUE, fill_cont_mean = T, fill_cat_mode = T, remove_outlier = TRUE, method = "nlminb", silent = TRUE, savedir = "Report/Final_run_no_outlier") # Bulinus delta-pl
