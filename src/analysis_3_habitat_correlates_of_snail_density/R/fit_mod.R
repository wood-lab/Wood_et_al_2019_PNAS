# Function to fit data from the wood-adams delta model

fit_mod <- function(model = 2, incl_disease = 0, bulinus = T, BIC_sel = T, fill_cont_mean = F, fill_cat_mode = F, remove_outlier = F, method = "Nelder-Mead", silent = TRUE, savedir = "Report", save_file = FALSE){

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

  file_name_base <- paste(c("biom", "bulinus")[as.numeric(bulinus)+1],"model",  model, c("AIC", "BIC")[as.numeric(BIC_sel)+1], c("no_mean_cont_fill", "mean_cont_fill")[as.numeric(fill_cont_mean)+1], c("no_mode_cat_fill", "mode_cat_fill")[as.numeric(fill_cat_mode)+1], "remove_outlier",remove_outlier, sep = "-")


  # Step 0 -- Load data, create mesh, and assign to list for TMB
  data <- read_dat(exclude = T,  fill_cont_mean = fill_cont_mean, na_keep = T, quad = F, fill_cat_mode = fill_cat_mode, scale_cont = T, bulinus = bulinus, remove_outlier = remove_outlier) # read data in


  # data <- fit_mesh(data, cluster = T, n_knots = 15) # Fit mesh
  data_list <- build_data(data, model, debug = 0, bulinus = bulinus, parasite = 0) # Create for TMB


  # Step 1 -- Set initial values for parameters
  params <- build_params(model, data_list, incl_disease)
  random = c("epsilon_mat", "gamma_q", "x_q_cont_missing")#, "omega_s", "omega_st")


  # Step 3 -- make and compile template file
  setwd("src")
  library(TMB)
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
  mod_sel <- aic_sel( data = data_list, params = params, version = version, map = map, random = random, BIC_sel = BIC_sel, method = method, silent = silent)


  # Extract final objects
  Final_Opt <- mod_sel$mod_list[[length(mod_sel$mod_list)]]
  Final_Obj <- mod_sel$Obj_list[[length(mod_sel$Obj_list)]]
  Final_Terms <- mod_sel$terms_list[[length(mod_sel$terms_list)]]
  mod_sel_results <- mod_sel$summ_stats


  # -----------------------------------------------------------------------------------------------------------------------------
  # STEP 3: Report results
  # -----------------------------------------------------------------------------------------------------------------------------
  pred <- Final_Obj$report(Final_Obj$env$last.par.best)
  rep <- sdreport(Final_Obj)
  dat <- list( mod_sel = mod_sel, data_list = data_list, data = data, params = params, map = map, bulinus = bulinus, model = model, fill_cont_mean = fill_cont_mean, fill_cat_mode = fill_cat_mode, BIC_sel = BIC_sel, file_name_base = file_name_base, savedir = savedir)

  if(save_file){
    save(dat , file = paste0(savedir,"/model_selection", file_name_base, ".RData"))
  }


  # Plot partial response curves:
  plot_coef(  Final_Opt, Final_Obj, Final_Terms, data_list, data, bulinus, model , fill_cont_mean, fill_cat_mode, BIC_sel, pred, file_name_base, savedir)

  # Variograms
  # run_variogram( data, data_list, pred, bulinus, model )


  # Plot predicted vs observed
  filename <- paste(savedir, "/observed_vs_predicted_count_", file_name_base, ".png", sep = "")
  png( file = filename , width=8, height=6, res=200, units='in')
  # par(mfrow = c(1,2))
  plot(x = data_list$c_q, y = pred$pred_cq, xlab = "Observed count", ylab = "Predicted count", main = "ZIPL")
  abline(0,1)
  #plot(x = data_list$y_i, y = pred$zero_prob_yi, xlab = "Observed infection", ylab = "Predicted infection", main = NA)
  dev.off()


  # Model selection results
  filename <- paste(savedir,"/model_selection_", file_name_base, ".csv", sep = "")
  write.csv(mod_sel_results, file = filename)


  # Save parameters
  filename <- paste(savedir, "/parameters_for_", file_name_base, ".RData", sep = "")
  save(rep, file = filename)

  param_names <-  names(rep$value)
  param_names[which(param_names %in% c("beta_p", "beta_c"))] <- paste(param_names[which(param_names %in% c("beta_p", "beta_c"))], colnames(data_list$x_q), sep = "-")
  params_sdreport <- data.frame( names = param_names, values = rep$value, sd = rep$sd)

  filename <- paste(savedir,"/parameters_for_", file_name_base, ".csv", sep = "")
  write.csv(params_sdreport, file = filename)

  # Site random effects
  intercepts <- Final_Opt$SD$par.fixed[which(names(Final_Opt$SD$par.fixed) == "beta0")]
  epsilon <- params$epsilon_mat[1:2,]
  epsilon_est <- rep$par.random[which(names(rep$par.random) == "epsilon_mat")]
  epsilon <- replace(epsilon, values = epsilon_est)
  site_id <- as.data.frame(unique( cbind( data$density_data$site, data_list$s_q ) ))
  colnames(site_id) <- c("site", "s_q")
  site_id$s_q <- as.numeric(as.character(site_id$s_q))
  site_id <- site_id[order(site_id$s_q),]
  site_id$epsilon_s_p <- epsilon[1,]
  site_id$epsilon_s_c <- epsilon[2,]
  site_id$partial_response <- (1/(1+exp(- (intercepts[1] + site_id$epsilon_s_p)))) * exp(intercepts[2] + site_id$epsilon_s_c)
  site_id <- rbind(as.matrix(site_id), c("Global", NA, 0, 0, (1/(1+exp(- (intercepts[1])))) * exp(intercepts[2])))


  filename <- paste(savedir, "/site_parameters_for_", file_name_base, ".csv", sep = "")
  write.csv(site_id, file = filename)
}

#
# # Run models - With outlier
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = T, fill_cont_mean = F, fill_cat_mode = F, remove_outlier = F, method = "Nelder-Mead", silent = TRUE, savedir = "Report/Nov_bulinus_with_outlier-Nelder-Mead") # Bulinus delta-pl
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = T, fill_cont_mean = T, fill_cat_mode = T, remove_outlier = F, method = "Nelder-Mead", silent = TRUE, savedir = "Report/Nov_bulinus_with_outlier-Nelder-Mead") # Bulinus delta-pl
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = T, fill_cont_mean = F, fill_cat_mode = T, remove_outlier = F, method = "Nelder-Mead", silent = TRUE, savedir = "Report/Nov_bulinus_with_outlier-Nelder-Mead") # Bulinus delta-pl
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = T, fill_cont_mean = T, fill_cat_mode = F, remove_outlier = F, method = "Nelder-Mead", silent = TRUE, savedir = "Report/Nov_bulinus_with_outlier-Nelder-Mead") # Bulinus delta-pl
#
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = T, fill_cont_mean = F, fill_cat_mode = F, remove_outlier = F, method = "L-BFGS-B", silent = TRUE, savedir = "Report/Nov_bulinus_with_outlier-L-BFGS-B") # Bulinus delta-pl
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = T, fill_cont_mean = T, fill_cat_mode = T, remove_outlier = F, method = "L-BFGS-B", silent = TRUE, savedir = "Report/Nov_bulinus_with_outlier-L-BFGS-B") # Bulinus delta-pl
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = T, fill_cont_mean = F, fill_cat_mode = T, remove_outlier = F, method = "L-BFGS-B", silent = TRUE, savedir = "Report/Nov_bulinus_with_outlier-L-BFGS-B") # Bulinus delta-pl
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = T, fill_cont_mean = T, fill_cat_mode = F, remove_outlier = F, method = "L-BFGS-B", silent = TRUE, savedir = "Report/Nov_bulinus_with_outlier-L-BFGS-B") # Bulinus delta-pl
#
# # Run models - No outlier
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = T, fill_cont_mean = F, fill_cat_mode = F, remove_outlier = T, method = "Nelder-Mead", silent = TRUE, savedir = "Report/Nov_bulinus_no_outlier-Nelder-Mead") # Bulinus delta-pl
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = T, fill_cont_mean = T, fill_cat_mode = T, remove_outlier = T, method = "Nelder-Mead", silent = TRUE, savedir = "Report/Nov_bulinus_no_outlier-Nelder-Mead") # Bulinus delta-pl
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = T, fill_cont_mean = F, fill_cat_mode = T, remove_outlier = T, method = "Nelder-Mead", silent = TRUE, savedir = "Report/Nov_bulinus_no_outlier-Nelder-Mead") # Bulinus delta-pl
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = T, fill_cont_mean = T, fill_cat_mode = F, remove_outlier = T, method = "Nelder-Mead", silent = TRUE, savedir = "Report/Nov_bulinus_no_outlier-Nelder-Mead") # Bulinus delta-pl
#
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = T, fill_cont_mean = F, fill_cat_mode = F, remove_outlier = T, method = "L-BFGS-B", silent = TRUE, savedir = "Report/Nov_bulinus_no_outlier-L-BFGS-B") # Bulinus delta-pl
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = T, fill_cont_mean = T, fill_cat_mode = T, remove_outlier = T, method = "L-BFGS-B", silent = TRUE, savedir = "Report/Nov_bulinus_no_outlier-L-BFGS-B") # Bulinus delta-pl
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = T, fill_cont_mean = F, fill_cat_mode = T, remove_outlier = T, method = "L-BFGS-B", silent = TRUE, savedir = "Report/Nov_bulinus_no_outlier-L-BFGS-B") # Bulinus delta-pl
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = T, fill_cont_mean = T, fill_cat_mode = F, remove_outlier = T, method = "L-BFGS-B", silent = TRUE, savedir = "Report/Nov_bulinus_no_outlier-L-BFGS-B") # Bulinus delta-pl
#
#
#
#
# # Run models - AIC
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = FALSE, fill_cont_mean = F, fill_cat_mode = F, remove_outlier = F, method = "nlminb", silent = TRUE, savedir = "Report/Nov_bulinus_without_outlier-nlminb-AIC") # Bulinus delta-pl
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = FALSE, fill_cont_mean = T, fill_cat_mode = T, remove_outlier = F, method = "nlminb", silent = TRUE, savedir = "Report/Nov_bulinus_without_outlier-nlminb-AIC") # Bulinus delta-pl
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = FALSE, fill_cont_mean = F, fill_cat_mode = T, remove_outlier = F, method = "nlminb", silent = TRUE, savedir = "Report/Nov_bulinus_without_outlier-nlminb-AIC") # Bulinus delta-pl
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = FALSE, fill_cont_mean = T, fill_cat_mode = F, remove_outlier = F, method = "nlminb", silent = TRUE, savedir = "Report/Nov_bulinus_without_outlier-nlminb-AIC") # Bulinus delta-pl
#
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = FALSE, fill_cont_mean = F, fill_cat_mode = F, remove_outlier = TRUE, method = "nlminb", silent = TRUE, savedir = "Report/Nov_bulinus_with_outlier-nlminb-AIC") # Bulinus delta-pl
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = FALSE, fill_cont_mean = T, fill_cat_mode = T, remove_outlier = TRUE, method = "nlminb", silent = TRUE, savedir = "Report/Nov_bulinus_with_outlier-nlminb-AIC") # Bulinus delta-pl
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = FALSE, fill_cont_mean = F, fill_cat_mode = T, remove_outlier = TRUE, method = "nlminb", silent = TRUE, savedir = "Report/Nov_bulinus_with_outlier-nlminb-AIC") # Bulinus delta-pl
# fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = FALSE, fill_cont_mean = T, fill_cat_mode = F, remove_outlier = TRUE, method = "nlminb", silent = TRUE, savedir = "Report/Nov_bulinus_with_outlier-nlminb-AIC") # Bulinus delta-pl
#
#
# # Final model
# dpl <- fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = TRUE, fill_cont_mean = F, fill_cat_mode = F, remove_outlier = TRUE, method = "nlminb", silent = TRUE, savedir = "Report/March_DPL_final_no_outlier") # Bulinus delta-pl
# zipl <- fit_mod(model = 6, incl_disease = 0, bulinus = T, BIC_sel = TRUE, fill_cont_mean = F, fill_cat_mode = F, remove_outlier = TRUE, method = "nlminb", silent = TRUE, savedir = "Report/March_ZIPL_final_no_outlier") # Bulinus z-pl
#
# # fit_mod(model = 2, incl_disease = 0, bulinus = T, BIC_sel = T, fill_cont_mean = F, fill_cat_mode = F) # Bulinus delta-pl - WORKSb
# # # fit_mod(model = 6, incl_disease = 0, bulinus = T, BIC_sel = T, fill_cont_mean = F, fill_cat_mode = F) # Bulinus ZI-PL - WORKS
# # fit_mod(model = 0, incl_disease = 0, bulinus = F, BIC_sel = T, fill_cont_mean = F, fill_cat_mode = F) # Biomp delta-lognormal - DOES NOT CONVERGE
# # fit_mod(model = 1, incl_disease = 0, bulinus = F, BIC_sel = T, fill_cont_mean = F, fill_cat_mode = F) # Biomp delta-poisson - DOES NOT CONVERGE
# # fit_mod(model = 2, incl_disease = 0, bulinus = F, BIC_sel = T, fill_cont_mean = F, fill_cat_mode = F) # Biomp delta-lognormal poisson - WORKS
# # #  fit_mod(model = 4, incl_disease = 0, bulinus = F, BIC_sel = T, fill_cont_mean = T, fill_cat_mode = F) # Biomp ZI-lognormal - DOES NOT WORK
# # # fit_mod(model = 5, incl_disease = 0, bulinus = F, BIC_sel = T, fill_cont_mean = T, fill_cat_mode = F) # Biomp ZI-poisson - DOES NOT WORK
# # fit_mod(model = 6, incl_disease = 0, bulinus = F, BIC_sel = T, fill_cont_mean = T, fill_cat_mode = F) # Biomp ZI-PL - WORKS
# # fit_mod(model = 8, incl_disease = 0, bulinus = F, BIC_sel = T, fill_cont_mean = T, fill_cat_mode = F) # Biomp tweedie - WORKS
