# Function to fit an spde mesh to within-site data

fit_mesh <- function(density_spatio, density_temporal){
  library( INLA )
  library( RandomFields )
  library( RANN )
  library( sp )
  library( rgdal )

  data <- cbind(density_spatio, density_temporal)

  # Convert lat and long to UTM meters
  coordinates(data) <- c("long", "lat")
  proj4string(data) <- sp::CRS("+proj=longlat +datum=WGS84")
  data <- sp::spTransform(data, CRS("+proj=utm +zone=28 ellps=WGS84"))
  data <- as.data.frame(data)


  # Create an identifier for each site and time period
  ts_ids <- c("sampling_site", "field_mission")
  ts_df <- unique(data[,ts_ids])
  ts_df$q_ts <- c(1:nrow(ts_df)) - 1
  data <- merge(data,  ts_df, by = ts_ids)


  # Create an SPDE mesh for each site and field mission
  m0 <- list()
  m1 <- list()
  m2 <- list()

  for(i in ts_df$q_ts){
    dat_sub <- data[which(data$q_ts == i),]
    loc_xy <- data.frame(dat_sub[,c("long", "lat")])

    mesh = inla.mesh.create( loc_xy )
    spde = inla.spde2.matern(mesh)

    m0[[i+1]] = spde$param.inla$M0
    m1[[i+1]] = spde$param.inla$M1
    m2[[i+1]] = spde$param.inla$M2
  }

  # Convert sparse matrices to arrays and get dimensions
  max_nrow <- max(sapply(m0, nrow))
  max_ncol <- max(sapply(m0, ncol))

  M0_s <- array(0, dim = c(max_nrow, max_ncol, length(ts_df$q_ts)))
  M1_s <- array(0, dim = c(max_nrow, max_ncol, length(ts_df$q_ts)))
  M2_s <- array(0, dim = c(max_nrow, max_ncol, length(ts_df$q_ts)))

  dim_xy <- data.frame(matrix(NA, nrow = length(ts_df$q_ts), ncol = 2)) # Number of rows and columns
  colnames(dim_xy) = c("nrow", "ncol")
  dim_xy$q_ts <- ts_df$q_ts

  for(i in 1:length(m0)){
    dim_xy$nrow[i] <- nrow(m0[[i]])
    dim_xy$ncol[i] <- ncol(m0[[i]])
    for(x in 1:dim_xy$nrow[i]){
      for(y in 1:dim_xy$ncol[i]){
        M0_s[x, y, i] <- m0[[i]][x, y]
        M1_s[x, y, i] <- m1[[i]][x, y]
        M2_s[x, y, i] <- m2[[i]][x, y]
      }
    }
  }

  meshes <- list(M0_s, M1_s, M2_s, dim_xy)
  return(meshes)
}
