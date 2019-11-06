# Function to read the data used in the analysis. Also conducts preliminary data analysis on covariates. Returns a list of the FULL data to be used in TMB. Covariates include "lake_river" where lake = 1, river = 0. Exlude will exlude treatments where a fence of prawn were included. fill_cont_mean will fill missing values with the mean if a continuous. na_keep is a switch to remove covariates that have na. Quad adds quadratic variables to the model

# Filling missing variables is as follows: 0 - leave as NA, 1 - fill in as mean/mode across all variables, 2 - fill in mean/mode of site specific variable

Mode <- function(x, na.rm = T) {
  ux <- unique(x)

  if(na.rm == T){
    if(length(which(is.na(ux))) > 0){
      ux <- ux[-which(is.na(ux))]
    }
  }
  ux[which.max(tabulate(match(x, ux)))]
}

read_dat <- function(exclude = T, fill_cont_mean =0, na_keep = T, quad = F, scale_cont = T, fill_cat_mode = 0, bulinus = T, outlier = T){

  source("R/support/fixNAs.R")

  # Step 0 -- Load data and packages
  density_data <- read.csv("data/aggregated_data_with_treatments_indicated.csv", na.strings=c(""))
  individual_data <- read.csv("data/merged_data.csv", na.strings=c(""))
  lake_river <- data.frame(readxl::read_xlsx("data/Info on lake_river for Grant.xlsx", sheet = 1))
  site_area <- data.frame(readxl::read_xlsx("data/PrawnBoundaryAreas_clwupdated.xlsx"))

  # Emergent vs floating interaction
  density_data$emergent_vs_floating <- as.numeric(as.character(density_data$emergent_stems_cleaned)) * as.numeric(as.character(density_data$floating_veg_g_cleaned))

  # Step 1 -- Convert factors to characters
  density_data[] <- lapply(density_data, as.character)
  individual_data[] <- lapply(individual_data, as.character)
  site_area[] <- lapply(site_area, as.character)
  colnames(density_data)[colnames(density_data) == "sampling_site"] <- "site"
  colnames(site_area)[colnames(site_area) == "FM"] <- "field_mission"
  colnames(site_area)[colnames(site_area) == "site"] <- "site_short"


  # Step 1 -- Remove incomplete spatial data
  # density_data$long <- as.numeric(density_data$long)
  # density_data$lat <- as.numeric(density_data$lat)
  # density_data <- density_data[complete.cases(density_data[,c("lat", "long")]),]

  # Remove outlier
  if(outlier == FALSE){
    density_data <- density_data[which(as.numeric(as.character(density_data$snail_is_trunc_globo)) < 700), ]
  }


  # Step 2 -- Set snail height to 4 if >4 or <4
  individual_data$lab_snail_height_mm[which(individual_data$lab_snail_height_mm %in% c("< 4", ">4"))] <- 4
  individual_data$lab_snail_height_mm <- scale(as.numeric(as.character(individual_data$lab_snail_height_mm)))
  individual_data <- individual_data[complete.cases(individual_data[,"lab_snail_height_mm"]),]


  # Step 3 -- Merge lake_river identifiers with other data
  lake_river$lake_river <- ifelse( lake_river$lake_river == "lake", 1, 0) # Make it a dummy variable where lake = 1, river = 0
  individual_data <- merge(individual_data, lake_river, by = "site")
  density_data <- merge(density_data, lake_river, by = "site")

  # Merge Site area
  density_data$site_short <- gsub(" ", "", density_data$site)
  individual_data$site_short <- gsub(" ", "", individual_data$site)
  individual_data <- merge(individual_data, site_area, by = c("site_short" , "field_mission"))
  density_data <- merge(density_data, site_area, by = c("site_short" , "field_mission"))


  # Step 4 -- Schisto data treatment
  individual_data$cerc_molecular_ID  <- as.character(individual_data$cerc_molecular_ID )
  individual_data$cerc_molecular_ID <- ifelse( is.na(individual_data$cerc_molecular_ID), "none", individual_data$cerc_molecular_ID)

  # -- 4.1. Make a column indicating when you have molecular confirmation of hybrid identity
  individual_data$molec_conf_hybrid <- ifelse(individual_data$cerc_molecular_ID=="Schistosoma haematobium/bovis hybrid"|individual_data$cerc_molecular_ID=="Schistosoma haematobium/bovis hybrid + Plagiorchis + Echinostoma"|individual_data$cerc_molecular_ID=="hybrid", 1, 0)

  # -- 4.2. Make a column indicating when you have molecular confirmation of S.haem identity
  individual_data$molec_conf_haem <- ifelse(individual_data$cerc_molecular_ID=="Schistosoma haematobium"|individual_data$cerc_molecular_ID=="Some CO1 S haematobium + paramphistome"|individual_data$cerc_molecular_ID=="Aptorchis + Paraamphistome (probably metacerc) + Schistosoma haematobium", 1, 0)

  # -- 4.3. Make a column indicating when you have molecular confirmation of S.mansoni identity
  individual_data$molec_conf_mansoni <- ifelse(individual_data$cerc_molecular_ID=="Schistosoma mansoni"|individual_data$cerc_molecular_ID=="Schistsoma mansoni",1,0)


  # Step 5 -- Names of data we want (double check these with chelsea)
  continuous_density_covariate_names <- c("emergent_stems_cleaned", "floating_veg_g_cleaned", "emergent_vs_floating", "secchi_depth_cm", "pH", "nitrate_mg.L_NO3", "nitrite_mg.L_NO2", "phosphate_mg.L_PO4", "depth_cm_cleaned_agg") # No longer include = "oxygen_mg.L_o2", "salinity_brix", "flowmeter_rotations", "temp_C", "area_prawn_m2_UTM28N
  categorical_density_covariate_names <- c("lake_river", "typha", "phrag", "cyperus", "poa_grass", "ludwigia", "nymphaea", "salvinia", "ceratophyllum", "potamogeton", "pistia", "azolla", "ipomea", "lemna")

  density_names <- c("snail_is_trunc_globo", "snail_is_biomph", "point_ID", "lat", "long", "unique_quad", "site", "sampling_date", "field_mission", "sampling_start_time_cleaned", "include_in_which_analyses")

  # Add quadratic relationship
  if(quad == T){
    continuous_quad <- sapply(density_data[,continuous_density_covariate_names], function(x) as.numeric(as.character(x))) # make numeric
    continuous_quad = apply(continuous_quad, c(1,2), function(x) x^2)

    colnames(continuous_quad) <- paste(continuous_density_covariate_names, "squared", sep ="_")

    continuous_density_covariate_names <- c(continuous_density_covariate_names, colnames(continuous_quad) )
    density_data <- cbind(density_data, continuous_quad)
  }

  density_names <- c(continuous_density_covariate_names, categorical_density_covariate_names, density_names)


  individual_covariate_names <- c("lab_snail_height_mm")
  individual_names <- c("snail_visual_ID", "point_ID", "lat", "long", "site", "field_mission", "sampling_date", "sampling_start_time_cleaned")
  individual_response_names <- c("molec_conf_haem", "molec_conf_hybrid", "molec_conf_mansoni") # was the snail infected with schisto? 1 = yes, 0 = no, na = if snail was not diagnosable - e.g., if it was dead and decomposed - or if the snail was not a host
  individual_names <- c(individual_covariate_names, individual_names, individual_response_names)


  # Step 6 -- Subset the individual data for bulinus or biomphalaria
  species_i = individual_data[,"snail_visual_ID"]
  if(bulinus == T){
    bulinus_snails <- which( species_i %in% c("trunc_globo", "globosus", "truncatus"))
    individual_data <- individual_data[bulinus_snails,]
  }
  # Biomphalaria
  if(bulinus == F){
    biomph_snails <- which( species_i %in% c("biomph"))
    individual_data <- individual_data[biomph_snails,]
  }


  # Step 7 -- Extract data
  density_data <- density_data[,which(colnames(density_data) %in% density_names)]
  individual_data <- individual_data[,which(colnames(individual_data) %in% individual_names)]

  # Change column names for indexing later
  colnames(density_data)[which(colnames(density_data) %in% c(continuous_density_covariate_names))] <- paste0("cont_covar_", colnames(density_data)[which(colnames(density_data) %in% c(continuous_density_covariate_names))])
  colnames(density_data)[which(colnames(density_data) %in% c(categorical_density_covariate_names))] <- paste0("cat_covar_", colnames(density_data)[which(colnames(density_data) %in% c(categorical_density_covariate_names))])

  continuous_density_covariate_names <-  paste0("cont_covar_", c(continuous_density_covariate_names))
  categorical_density_covariate_names <-  paste0("cat_covar_", c(categorical_density_covariate_names))


  # -- 8.1. Individual quadrate index
  quad_ids <- c('point_ID','site', "sampling_date", "field_mission")


  # -- 8.2. Exluded data we want to
  exclude_rows <- c()
  if( exclude == T ){
    exclude_rows <- which(density_data$include_in_which_analyses == "exclude")
  }

  # Turn variables into numeric
  for(i in 1:length(categorical_density_covariate_names)){
    density_data[,categorical_density_covariate_names[i]] <- as.numeric(as.character(density_data[,categorical_density_covariate_names[i]]))
  }

  for(i in 1:length(continuous_density_covariate_names)){
    density_data[,continuous_density_covariate_names[i]] <- as.numeric(as.character(density_data[,continuous_density_covariate_names[i]]))
  }

  # Data summary
  names <- c(continuous_density_covariate_names, categorical_density_covariate_names)
  percent_na <- c(colSums(is.na(density_data[-exclude_rows,continuous_density_covariate_names])), colSums(is.na(density_data[-exclude_rows ,categorical_density_covariate_names]))) / nrow(density_data)
  mean_mode <- c(apply(density_data[-exclude_rows,continuous_density_covariate_names],2, mean, na.rm = T), apply(density_data[-exclude_rows, categorical_density_covariate_names],2, Mode, na.rm = T))
  sd_vals <- c(apply(density_data[-exclude_rows,continuous_density_covariate_names],2, sd, na.rm = T), rep(NA, length(categorical_density_covariate_names)))

  dat_summ <- data.frame( names = names, mean_mode = mean_mode, sd_vals = sd_vals, percent_na = percent_na)
  write.csv(dat_summ, file = paste("Report/data_summary_for_", c("biom", "bulinus")[as.numeric(bulinus)+1], c("_no_exclude", "_exclude")[as.numeric(exclude)+1], ".csv", sep = ""))


  ####################################################################################
  # Remove NAs
  if(na_keep == F){
    density_data <- density_data[complete.cases(density_data[,continuous_density_covariate_names]),]

    density_data <- density_data[complete.cases(density_data[,categorical_density_covariate_names]),]
  }


  # Fill with mode across all data?
  if(fill_cat_mode == 1){
    for(i in 1:length(categorical_density_covariate_names)){

      density_data[which(is.na(density_data[,categorical_density_covariate_names[i]])),categorical_density_covariate_names[i]] <- Mode(density_data[which(!is.na(density_data[-exclude_rows, categorical_density_covariate_names[i]])),categorical_density_covariate_names[i]] )
    }
  }


  # Step - 9 Scale continuous covariates
  for(i in 1:length(continuous_density_covariate_names)){
    # Scale
    if(scale_cont == T){
      density_data[,continuous_density_covariate_names[i]] <- scale(density_data[,continuous_density_covariate_names[i]], center = mean(density_data[-exclude_rows,continuous_density_covariate_names[i]], na.rm = T), scale = sd(density_data[-exclude_rows,continuous_density_covariate_names[i]], na.rm = T))
    }

    #individual_data[,continuous_density_covariate_names[i]] <- scale(individual_data[,continuous_density_covariate_names[i]])
  }

  # Fill with median across sites (i.e. 0)
  if(fill_cont_mean == 1){
    for(i in 1:length(continuous_density_covariate_names)){
      density_data[which(is.na(density_data[,continuous_density_covariate_names[i]])),continuous_density_covariate_names[i]] <- median(density_data[-exclude_rows, continuous_density_covariate_names[i]], na.rm = T)
    }
  }


  # Fill with site specific mode?
  if(fill_cat_mode == 2){
    density_data <- fixNAs(density_data, categorical_density_covariate_names, exclude_rows)
  }

  if(fill_cont_mean == 2){
    density_data <- fixNAs(density_data, continuous_density_covariate_names, exclude_rows)
  }


  # Remove rows
  density_data <- density_data[-exclude_rows,]


  # Merge environmental data
  individual_data <- merge(individual_data,  density_data[, c(continuous_density_covariate_names, categorical_density_covariate_names, quad_ids)], by = quad_ids)


  data = list( density_data = density_data, individual_data = individual_data, dat_summ = dat_summ)

  return(data)
}
