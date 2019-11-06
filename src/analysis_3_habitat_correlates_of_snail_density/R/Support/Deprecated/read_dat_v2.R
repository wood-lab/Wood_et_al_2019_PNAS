# Function to read the data used in the analysis. Also conducts preliminary data analysis on covariates. Returns a list of the FULL data to be used in TMB. Covariates include "lake_river" where lake = 1, river = 0. Allows for specifying which parasite you want to use if bulinus (0 = both, 1 = haem, 2 = hybrid). Exlude will exlude treatments where a fence of prawn were included.

read_dat <- function(pda = F, bulinus = T, parasite = 0, exclude = T){

  # Step 0 -- Load data and packages
  library(lattice)
  library(readxl)

  density_data <- read.csv("data/aggregated_data_with_treatments_indicated.csv", na.strings=c(""))
  individual_data <- read.csv("data/Wood_et_al_individual_snail_data_5_12_2018.csv", na.strings=c(""))
  lake_river <- data.frame(readxl::read_xlsx("data/Info on lake_river for Grant.xlsx", sheet = 1))


  # Step 0 -- Convert factors to characters
  density_data[] <- lapply(density_data, as.character)
  individual_data[] <- lapply(individual_data, as.character)


  # Step 1 -- Remove na rows of density data
  complete_case_col <- c("lat", "long", "floating_veg_g_cleaned", "emergent_stems_cleaned")

  density_data$floating_veg_g_cleaned <- scale(as.numeric(density_data$floating_veg_g_cleaned))
  density_data$emergent_stems_cleaned <- scale(as.numeric(density_data$emergent_stems_cleaned))

  # Make into the mean
  density_data$emergent_stems_cleaned[which(is.na(density_data$emergent_stems_cleaned))] <- 0
  density_data$floating_veg_g_cleaned[which(is.na(density_data$floating_veg_g_cleaned))] <- 0

  density_data$long <- as.numeric(density_data$long)
  density_data$lat <- as.numeric(density_data$lat)
  density_data <- density_data[complete.cases(density_data[,complete_case_col]),]


  # Step 2 -- Remove na rows of individual data
  individual_data$floating_veg_g_cleaned <- scale(as.numeric(as.character(individual_data$floating_veg_g_cleaned)))
  individual_data$emergent_stems_cleaned <- scale(as.numeric(as.character(individual_data$emergent_stems_cleaned)))

  # -- Set snail height to 4 if >4 or <4
  individual_data$lab_snail_height_mm[which(individual_data$lab_snail_height_mm %in% c("< 4", ">4"))] <- 4
  individual_data$lab_snail_height_mm <- scale(as.numeric(as.character(individual_data$lab_snail_height_mm)))
  complete_case_col <- c(complete_case_col, "lab_snail_height_mm")
  individual_data <- individual_data[complete.cases(individual_data[,complete_case_col]),]


  # Step 3 -- Merge lake_river identifiers with other data
  lake_river$lake_river <- ifelse( lake_river$lake_river == "lake", 1, 0) # Make it a dummy variable where lake = 1, river = 0
  individual_data <- merge(individual_data, lake_river, by = "site")
  colnames(lake_river) = c("sampling_site", "lake_river", "cluster")
  density_data <- merge(density_data, lake_river, by = "sampling_site")


  # Step 4 -- Schisto data treatment
  individual_data$cerc_molecular_ID  <- as.character(individual_data$cerc_molecular_ID )
  individual_data$cerc_molecular_ID <- ifelse( is.na(individual_data$cerc_molecular_ID), "none", individual_data$cerc_molecular_ID)

  # -- 4.1. Make a column indicating when you have molecular confirmation of hybrid identity
  individual_data$molec_conf_hybrid <- ifelse(individual_data$cerc_molecular_ID=="Schistosoma haematobium/bovis hybrid"|individual_data$cerc_molecular_ID=="Schistosoma haematobium/bovis hybrid + Plagiorchis + Echinostoma"|individual_data$cerc_molecular_ID=="hybrid", 1, 0)

  # -- 4.2. Make a column indicating when you have molecular confirmation of S.haem identity
  individual_data$molec_conf_haem <- ifelse(individual_data$cerc_molecular_ID=="Schistosoma haematobium"|individual_data$cerc_molecular_ID=="Some CO1 S haematobium + paramphistome"|individual_data$cerc_molecular_ID=="Aptorchis + Paraamphistome (probably metacerc) + Schistosoma haematobium",1,0)

  # -- 4.3. Make a column indicating when you have molecular confirmation of S.mansoni identity
  individual_data$molec_conf_mansoni <- ifelse(individual_data$cerc_molecular_ID=="Schistosoma mansoni"|individual_data$cerc_molecular_ID=="Schistsoma mansoni",1,0)


  # Step 5 -- Names of data we want (double check these with chelsea)
  density_covariate_names <- c("emergent_stems_cleaned", "floating_veg_g_cleaned", "lake_river") # c("flowmeter_rotations", "secchi_depth_cm", "temp_C", "salinity_brix", "pH", "nitrate_mg.L_NO3", "nitrite_mg.L_NO2", "phosphate_mg.L_PO4", "oxygen_mg.L_o2", "depth_cm_cleaned_agg", "typha", "phrag", "cyperus", "poa_grass", "ludwigia", "nymphaea", "salvinia", "ceratophyllum", "potamogeton", "pistia", "azolla", "ipomea", "lemna", "emergent_stems_cleaned", "floating_veg_g_cleaned", "lake_river")
  density_response_names <- c("snail_is_trunc_globo", "snail_is_biomph")
  density_spatio_names <- c("point_ID", "lat", "long", "unique_quad", "sampling_site") # Lat, long are in WGS84
  density_temporal_names <- c("sampling_date", "field_mission", "sampling_start_time_cleaned", "include_in_which_analyses")
  density_names <- c(density_covariate_names, density_temporal_names, density_response_names, density_spatio_names, density_temporal_names)

  individual_covariate_names <- c("lab_snail_height_mm")
  individual_species <- c("snail_visual_ID")
  individual_response_names <- c("molec_conf_haem", "molec_conf_hybrid", "molec_conf_mansoni") # was the snail infected with schisto? 1 = yes, 0 = no, na = if snail was not diagnosable - e.g., if it was dead and decomposed - or if the snail was not a host
  individual_spatio_names <- c("point_ID", "lat", "long", "site")
  individual_temporal_names <- c("field_mission", "sampling_date", "sampling_start_time_cleaned")
  individual_names <- c(individual_covariate_names, density_covariate_names, individual_species, individual_response_names, individual_spatio_names, individual_temporal_names)


  # Step 6 -- Subset the individual data for bulinus or biomphalaria
  species_i = individual_data[,which(colnames(individual_data) %in% individual_species)]
  if(bulinus == T){
    bulinus_snails <- which( species_i %in% c("trunc_globo", "globosus", "truncatus"))
    individual_data <- individual_data[bulinus_snails,]
  }
  # Biomphalaria
  else {
    biomph_snails <- which( species_i %in% c("biomph"))
    individual_data <- individual_data[biomph_snails,]
  }


  # Step 7 -- Extract data
  density_data <- density_data[,which(colnames(density_data) %in% density_names)]
  individual_data <- individual_data[,which(colnames(individual_data) %in% individual_names)]


  # Step 8 -- Build indices
  # -- 8.1. Individual quadrate index
  quad_ids <- c('point_ID','sampling_site', "sampling_date", "field_mission")
  quad_df <- unique(density_data[,quad_ids])
  names(quad_df)[2] <- "site"
  quad_ids[2] <- "site"
  quad_df$quad_no = (1:nrow(quad_df)) - 1
  density_data$quad_no <- quad_df$quad_no
  individual_data <- merge(individual_data,  quad_df, by = quad_ids)

  # -- 8.2. Exluded data we want to
  if( exclude == T ){
    exclude_quads <- density_data$quad_no[which(density_data$include_in_which_analyses == "exclude")]
    density_data <- density_data[-which(density_data$quad_no %in% exclude_quads),]
    individual_data <- individual_data[-which(individual_data$quad_no %in% exclude_quads ),]
  }

  # -- 8.3. Site and village indices
  site_df <- data.frame(site = as.character(unique(density_data$sampling_site)),
                        site_no = (1:length(unique(density_data$sampling_site))) - 1
  )
  village_df <- data.frame(village = unique( gsub( '[[:digit:]]+', "", as.character( density_data$sampling_site ) ) ),
                           village_no = (1:length(unique(gsub( '[[:digit:]]+', "", as.character( density_data$sampling_site ) ) ) ) ) - 1
  )


  # Step 9 -- Assign data to list for use in TMB
  data <- list(
    n_s = length( unique( density_data$sampling_site ) ), # Number of sites
    n_t = length( unique( density_data$field_mission ) ),
    q_q = quad_df$quad_no, # quadrate ID
    t_q = as.numeric(density_data$field_mission) - 1,
    s_q = site_df[match( as.character( density_data$sampling_site ), as.character(site_df$site) ) , 2], # Site ID
    v_q = village_df[match( gsub( '[[:digit:]]+', "", as.character( density_data$sampling_site ) ), as.character(village_df$village) ) , 2], # Get village identifiers
    X_pq = as.matrix(density_data[,which(colnames(density_data) %in% density_covariate_names)]), # Covariates for presence/absence
    X_cq = as.matrix(density_data[,which(colnames(density_data) %in% density_covariate_names)]), # covariates for local density

    q_i = individual_data$quad_no, # unique quadrate ID
    s_i = site_df[match( as.character( individual_data$site ), as.character(site_df$site) ) , 2], # Site ID
    v_i = village_df[match( gsub( '[[:digit:]]+', "", as.character( individual_data$site ) ), as.character(village_df$village) ) , 2], # Get village identifiers
    x_i = as.matrix(individual_data[ , which(colnames(individual_data) %in% c(individual_covariate_names, density_covariate_names))]) # Covariates predicting snail infection
  )


  # Step 10 -- Extract the response variables of desire
  # -- 10.1. Bulinus
  if(bulinus == T){
    if(parasite == 0){ # Both parasite species
      data[["y_i"]] <- cbind(individual_data$molec_conf_haem, individual_data$molec_conf_hybrid)
    }
    if(parasite == 1){ # Haemotobium
      data[["y_i"]] <- cbind(individual_data$molec_conf_haem, individual_data$molec_conf_hybrid)
    }
    if(parasite == 2){ # Hybrid
      data[["y_i"]] <- cbind(individual_data$molec_conf_haem, individual_data$molec_conf_hybrid)
    }
    data[["c_q"]] <- as.numeric(density_data$snail_is_trunc_globo)

    data[["x_i"]] <- array(unlist(list(data[["x_i"]], data[["x_i"]])), dim = c(nrow(data[["x_i"]]), ncol(data[["x_i"]]), 2))
  }
  # -- 10.2. Biomphalaria
  if(bulinus == F){
    data[["c_q"]] <- as.numeric(density_data$snail_is_biomph)
    data[["y_i"]] <- as.matrix(individual_data$molec_conf_mansoni)
  }

  data[["density_data"]] <- density_data
  data[["site_df"]] <- site_df

  return(data)
}
