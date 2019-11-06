
  library( INLA )
  library( RandomFields )
  library( RANN )
  library( sp )
  library( rgdal )

  # Step 0 -- Load simulator and simulate data. See file for specifications
  source("R/read_dat.R")
  data_list <- read_dat(bulinus = 1)


  density_data <- data_list$density_data
  site_df <- data_list$site_df

  # Convert lat and long to UTM meters
  coordinates(density_data) <- c("long", "lat")
  proj4string(density_data) <- sp::CRS("+proj=longlat +datum=WGS84")
  density_data <- sp::spTransform(density_data, CRS("+proj=utm +zone=28 ellps=WGS84"))
  density_data <- as.data.frame(density_data)

  # Scale to 0 and 1
  density_data$lat_scaled <- (density_data$lat - min(density_data$lat))#/ (max(density_data$lat) - min(density_data$lat))
  density_data$long_scaled <- (density_data$long - min(density_data$long))#/ (max(density_data$long) - min(density_data$long))

  # Create and identifier for each site
  site <- site_df$site
  n_s <- nrow(site_df)


  for(i in 1:n_s){
    site_ind = site[i]
    dat_sub <- density_data[which(density_data$sampling_site == site_ind),]
    loc_xy_orig = loc_xy = data.frame(dat_sub[,c("long_scaled", "lat_scaled")])
    coordinates(loc_xy) <- c("long_scaled", "lat_scaled")
    meshbuilder()
  }

