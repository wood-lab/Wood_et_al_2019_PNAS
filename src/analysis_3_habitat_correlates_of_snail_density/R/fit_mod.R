# Function to fit data from the wood-adams delta model

fit_mod <- function(model = 2, incl_disease = 0, bulinus = T, BIC_sel = T, fill_cont_mean = 0, fill_cat_mode = 0, outlier = T, two_stage = T){

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
  source("R/Support/two_stage_aic_selection_v2_optimx.R")
  source("R/Support/aic_selection_v3.R")
  source("R/Support/plot_coef_v2.R")


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


  # Create files and directories
  filename_base <- paste(c("biom", "bulinus")[as.numeric(bulinus)+1],"-model",  model, c("-AIC", "-BIC")[as.numeric(BIC_sel)+1], "-mean_cont_fill_", fill_cont_mean, "-mode_cat_fill_", fill_cat_mode,"-outlier_", outlier, "-two_stage_",two_stage, sep = "")
  mainDir = paste0(getwd(), "/Report")
  subDir = paste0(c("Biom", "Bulinus")[as.numeric(bulinus)+1], "_", c("AIC", "BIC")[as.numeric(BIC_sel)+1], "_", Sys.Date())
  report_dir <- paste0(mainDir, "/",subDir)

  if (file.exists(paste(mainDir, subDir, "/", sep = "/", collapse = "/"))) {
    cat("subDir exists in mainDir and is a directory")
  } else if (file.exists(paste(mainDir, subDir, sep = "/", collapse = "/"))) {
    cat("subDir exists in mainDir but is a file")
    # you will probably want to handle this separately
  } else {
    cat("subDir does not exist in mainDir - creating")
    dir.create(file.path(mainDir, subDir))
    dir.create(file.path(report_dir, "Figures"))
    dir.create(file.path(report_dir, "Figures/Variograms"))
  }

  # -----------------------------------------------------------------------------------------------------------------------------
  # STEP 2
  # -----------------------------------------------------------------------------------------------------------------------------
  map <- build_map(model, data_list, incl_disease = incl_disease, params, space = 0, space_time = 0)

  # Fit 1 -- Backward model selection

    mod_sel <- aic_sel_two( data = data_list, params = params, version = version, map = map, random = random, BIC_sel = BIC_sel, bulinus = bulinus, wd = report_dir, file_name = filename_base, two_stage = two_stage)


  # Extract final objects
  Final_Opt <- mod_sel$final_opt
  Final_Obj <- mod_sel$final_obj
  Final_Terms <- mod_sel$terms_list[[length(mod_sel$terms_list)]]
  mod_sel_results <- mod_sel$summ_stats


  # -----------------------------------------------------------------------------------------------------------------------------
  # STEP 3: Report results
  # -----------------------------------------------------------------------------------------------------------------------------
  pred <- Final_Obj$report()
  rep <- TMB::sdreport(Final_Obj)
  dat <- list(Final_Terms = Final_Terms, Final_Opt = Final_Opt, Final_Obj = Final_Obj, Report = rep, Predicted = pred, Data_List = data_list, Raw_Data = data,  Model_Selection = mod_sel)


  # Plot partial response curves:
  plot_coef(  rep , Final_Obj, Final_Terms, data_list, data, bulinus, model , fill_cont_mean, fill_cat_mode, BIC_sel, file_name = filename_base, wd = report_dir )

  # Variograms
  # run_variogram( data, data_list, pred, bulinus, model, wd = report_dir, file_name = filename_base )


  # Plot predicted vs observed
  filename <- paste(report_dir, "/observed_vs_predicted_count_", filename_base, ".png", sep = "")
  png( file = filename , width=8, height=6, res=200, units='in')
  # par(mfrow = c(1,2))
  plot(x = data_list$c_q, y = pred$pred_cq, xlab = "Observed count", ylab = "Predicted count", main = "ZIPL", xlim = c(0,150), ylim = c(0,150))
  abline(0,1)
  #plot(x = data_list$y_i, y = pred$zero_prob_yi, xlab = "Observed infection", ylab = "Predicted infection", main = NA)
  dev.off()


  # Model selection results
  filename <- paste(report_dir, "/model_selection_", filename_base, ".csv", sep = "")
  write.csv(mod_sel_results, file = filename)


  # Save parameters
  filename <- paste(report_dir, "/parameters_for_", filename_base, ".Rdata", sep = "")
  save(rep, file = filename)

  #
  filename <- paste(report_dir, "/model_objects_for_", filename_base, ".Rdata", sep = "")
  save(dat, file = filename)

  # Site random effects
  intercepts <- rep$value[which(names(rep$value) == "beta0")][1:2]
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


  filename <- paste(report_dir,"/site_parameters_for_", filename_base, ".csv", sep = "")
  write.csv(site_id, file = filename)

  return( dat )
}

