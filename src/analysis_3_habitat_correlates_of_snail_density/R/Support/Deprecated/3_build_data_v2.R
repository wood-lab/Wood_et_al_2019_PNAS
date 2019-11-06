# Function to build data list for TMB. Allows for specifying which parasite you want to use if bulinus (0 = both, 1 = haem, 2 = hybrid). Bulinus selects what snail species

build_data <- function(data, model, debug = 1, bulinus = T, parasite = 0){

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

  data_list$x_q = as.matrix( cbind(density_data[,grep("cont_covar_", colnames(density_data))], density_data[,grep("cat_covar_", colnames(density_data))]))

  # Index missing data
  data_list$n_x_q_cont <- length(grep("cont_covar_", colnames(data_list$x_q)))

  data_list[["model"]] = model

  return(data_list)
}
