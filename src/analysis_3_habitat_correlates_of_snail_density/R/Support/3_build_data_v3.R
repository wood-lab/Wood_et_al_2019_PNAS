# Function to build data list for TMB. Allows for specifying which parasite you want to use if bulinus (0 = both, 1 = haem, 2 = hybrid). Bulinus selects what snail species

build_data <- function(data, model, debug = 1, bulinus = T, parasite = 0){
  library(gtools)
  density_data <- data$density_data
  individual_data <- data$individual_data


  # STEP 1 -- BUILD INDICES
  density_data$quad_no <- c(1:nrow(density_data))-1
  individual_data <- merge(individual_data,  density_data[, c('point_ID','site', "sampling_date", "field_mission", "quad_no")], by = c('point_ID','site', "sampling_date", "field_mission"))

  site_df <- data.frame(site = as.character(unique(density_data$site)),
                        site_no = (1:length(unique(density_data$site))) - 1
  )

  village_df <- data.frame(village = unique( gsub( '[[:digit:]]+', "", as.character( density_data$site ) ) ),
                           village_no = (1:length(unique(gsub( '[[:digit:]]+', "", as.character( density_data$site ) ) ) ) ) - 1
  )

  field_mission_df <- data.frame(field_mission = as.character(unique(density_data$field_mission)),
                                 field_mission_no = (1:length(unique(density_data$field_mission))) - 1
  )


  # STEP 2 -- BUILD DATA LIST
  data_list <- list(
    # Settings
    model = model,
    incl_disease = 0,
    debug = debug,

    # -- 1.0. Indices
    n_s = length( unique( density_data$site ) ),
    n_t = length( unique( density_data$field_mission ) )
  )

  # -- 1.1 Local density data
  if( bulinus == T){ data_list$c_q = as.numeric(as.character(density_data$snail_is_trunc_globo))  }
  if( bulinus == F){ data_list$c_q = as.numeric( as.character( density_data$snail_is_biomph ) ) }
  data_list$q_q = density_data$quad_no
  data_list$t_q = site_df[match( as.character( density_data$field_mission ), as.character(field_mission_df$field_mission) ) , 2]
  data_list$s_q = site_df[match( as.character( density_data$site ), as.character(site_df$site) ) , 2]
  data_list$v_q = village_df[match( gsub( '[[:digit:]]+', "", as.character( density_data$site ) ), as.character(village_df$village) ) , 2]


  data_list$n_v = length( unique( data_list$v_q ))

  # Index to fit or not
  data_list$fit_ll = rep(1, length(density_data$quad_no)) # 1 = TRUE/FIT, 0 = FALSE/DON'T FIT

  # Get design matrix components
  x_i_cat <- density_data[,grep("cat_covar_", colnames(density_data))]
  x_i_cont <- density_data[,grep("cont_covar_", colnames(density_data))]

  data_list$x_q = as.matrix( cbind( x_i_cont , x_i_cat))

  # -- 1.2. Individual snail data
  data_list$x_i = as.matrix( cbind( individual_data[,grep("cont_covar_", colnames(individual_data))], individual_data[,grep("cat_covar_", colnames(individual_data))]))
  if(bulinus == T){
    if(parasite == 0){ # Both parasite species
      data_list$y_i = cbind(individual_data$molec_conf_haem, individual_data$molec_conf_hybrid)
      data_list$n_i = 2
      data_list$x_i = array(unlist(list(data_list[["x_i"]], data_list[["x_i"]])), dim = c(nrow(data_list[["x_i"]]), ncol(data_list[["x_i"]]), 2))
    }
    if(parasite == 1){ # Haemotobium
      data_list$y_i = as.matrix(individual_data$molec_conf_haem)
      n_i = 1
    }
    if(parasite == 2){ # Hybrid
      data_list$y_i = as.matrix(individual_data$molec_conf_hybrid)
      data_list$n_i = 1
    }
  }

  if(bulinus == F){
    data_list$y_i = as.matrix(individual_data$molec_conf_mansoni)
    data_list$n_i <- 1
  }

  data_list$n_p = dim(data_list$y_i)[2] # Number of parasites
  data_list$n_i = dim(data_list$y_i)[1] # Number of individual snail observations

  data_list$s_i = site_df[match( as.character( individual_data$site ), as.character(site_df$site) ) , 2]
  data_list$q_i = individual_data$quad_no
  data_list$v_i = village_df[match( gsub( '[[:digit:]]+', "", as.character( individual_data$site ) ), as.character(village_df$village) ) , 2] # Get village identifiers

  # # -- 1.3. SPDE components for each site
  # k_q = data$k_q,
  # dim_M_s = as.matrix(data$dim_M[,c("nrow", "ncol")]),
  # M0_s = data$M0_s,
  # M1_s = data$M1_s,
  # M2_s = data$M2_s

  # Index missing data
  data_list$n_x_q_cont <- length(grep("cont_covar_", colnames(data_list$x_q)))


  if( min(data_list$s_q) == 1){ data_list[["s_q"]] = data_list$s_q - 1}
  if( min(data_list$s_q) == 0){ data_list[["s_q"]] = data_list$s_q}

  if( min(data_list$s_i) == 1){ data_list[["s_i"]] = data_list$s_i - 1}
  if( min(data_list$s_i) == 0){ data_list[["s_i"]] = data_list$s_i}

  data_list[["model"]] = model

  return(data_list)
}
