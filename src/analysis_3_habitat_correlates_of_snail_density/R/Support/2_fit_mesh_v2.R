# Function to fit an spde mesh to within-site data. Does a site specific mesh.

fit_mesh <- function(data_list, cluster = F, n_knots = 15){
  library( INLA )
  library( RandomFields )
  library( RANN )
  library( sp )
  library( rgdal )

  density_data <- data_list$density_data

  # Convert lat and long to UTM meters
  coordinates(density_data) <- c("long", "lat")
  proj4string(density_data) <- sp::CRS("+proj=longlat +datum=WGS84")
  density_data <- sp::spTransform(density_data, CRS("+proj=utm +zone=28 ellps=WGS84"))
  density_data <- as.data.frame(density_data)


  # Create an identifier for each site and time period
  ts_ids <- c("site", "field_mission")
  ts_df <- unique(density_data[,ts_ids])
  ts_df$ts_q <- c(1:nrow(ts_df)) - 1
  density_data <- merge(density_data,  ts_df, by = ts_ids)


  # Create and identifier for each site
  site <- unique(density_data$site)
  n_s <- length(site)


  # Create an SPDE mesh for each site and field mission
  m0 <- list()
  m1 <- list()
  m2 <- list()

  # Minimum number of quadrats sampled per site and field mission
  min_q_s <- as.data.frame(table(density_data$site, density_data$field_mission)) # Count of quadrats in each field mission and site
  min_q_s <- min_q_s[which(min_q_s$Freq != 0),]
  min_q_s <- aggregate(Freq ~ Var1, min_q_s, function(x) min(x)) # Minimum count

  k_means_q <- c() # Index of which cluster a data point is associated with
  k_q <- c() # Index of which vertex a data point is associated with
  mesh_n <- c()


  for(i in 1:n_s){
    site_ind = site[i]
    dat_sub <- density_data[which(density_data$site == site_ind),]
    loc_xy_orig = loc_xy = data.frame(dat_sub[,c("long", "lat")])

    # Reduce number of stations -- OPTIONAL
    if(cluster == T){

      n_knots = min_q_s$Freq[which(min_q_s$Var1 == site_ind)]

      if( n_knots < nrow(loc_xy) ){
        knots_xy = kmeans( x=loc_xy, centers=n_knots )
        # Modify data
        loc_xy = knots_xy$centers
        k_means_q = as.numeric(knots_xy$cluster) # Get cluster number
      }
    }

    # boundary.loc <- loc_xy
    # boundary <- list(
    #   inla.nonconvex.hull(coordinates(boundary.loc), 4),
    #   inla.nonconvex.hull(coordinates(boundary.loc), 11.9))
    #
    # ## Build the mesh:
    # mesh <- inla.mesh.2d(boundary=boundary,
    #                      max.edge=c(15, 29.4),
    #                      min.angle=c(30, 21),
    #                      max.n=c(48000, 16000), ## Safeguard against large meshes.
    #                      max.n.strict=c(128000, 128000), ## Don't build a huge mesh!
    #                      cutoff=3, ## Filter away adjacent points.
    #                      offset=c(4, 11.9)) ## Offset for extra boundaries, if needed.

    # Get closest vertex to observation

    mesh <- INLA::inla.mesh.create( loc_xy )
    spde = INLA::inla.spde2.matern( mesh )

    A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(loc_xy))
    close_idx <- apply(A, 1, which.max)

    mesh_n <- c(mesh_n, mesh$n)

    if( cluster == F){
      k_q = c(k_q, close_idx) # Get vertex number
    }
    if( cluster == T){
      k_q = c(k_q, mesh$idx$loc[k_means_q])
    }

    # Visualize mesh and predictive process
    filename <- paste("Report/Figures/SPDE/SPDE_mesh_",site_ind, ".png",sep = "")
    png( file=filename, width=8, height=6, res=200, units='in')
    plot(mesh, xlab = paste(site_ind))
    if( cluster == T){points( loc_xy, cex=2, pch=3, col="green", lwd=5)}
    points( loc_xy_orig, cex=1.5, pch=20 )
    if( cluster == T){legend("topright", c("Quadrats", "Clusters"), pch = 20, col = c(1, "green"), bty = "n")}
    if( cluster == F){legend("topright", c("Quadrats"), pch = 20, col = c(1), bty = "n")}
    dev.off()

    m0[[i]] = spde$param.inla$M0
    m1[[i]] = spde$param.inla$M1
    m2[[i]] = spde$param.inla$M2
  }

  # Convert sparse matrices to arrays and get dimensions
  max_nrow <- max(sapply(m0, nrow))
  max_ncol <- max(sapply(m0, ncol))

  M0_s <- array(0, dim = c(max_nrow, max_ncol, n_s))
  M1_s <- array(0, dim = c(max_nrow, max_ncol, n_s))
  M2_s <- array(0, dim = c(max_nrow, max_ncol, n_s))

  dim_xy <- data.frame(matrix(NA, nrow = n_s, ncol = 2)) # Number of rows and columns
  colnames(dim_xy) = c("nrow", "ncol")
  dim_xy$site <- site

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

  # Assign to data
  data_list[["M0_s"]] <- M0_s
  data_list[["M1_s"]] <- M1_s
  data_list[["M2_s"]] <- M2_s
  data_list[["k_q"]] <- k_q
  data_list[["dim_M"]] <- dim_xy
  data_list[["mesh_n"]] <- mesh_n
  data_list[["density_data"]] <- density_data

  # Rearange dim
  data_list[["dim_M"]] <- merge(data_list[["dim_M"]], data_list$site_df, by = "site")
  data_list[["dim_M"]] <- data_list[["dim_M"]][order(data_list[["dim_M"]]$site_no), ]
  data_list[["dim_M"]] <- data_list[["dim_M"]][, c("nrow", "ncol", "site_no", "site")]


  return(data_list)
}
