# SPDE mesh specifications build using INLA:meshbuilder() for each site
for(i in 1:n_s){
  site_ind = site[i]
  dat_sub <- density_data[which(density_data$sampling_site == site_ind),]
  loc_xy_orig = loc_xy = data.frame(dat_sub[,c("long", "lat")])
  coordinates(loc_xy) <- c("long", "lat")

  boundary.loc <- loc_xy
  boundary <- list(
    inla.nonconvex.hull(coordinates(boundary.loc), 4),
    inla.nonconvex.hull(coordinates(boundary.loc), 11.9))

  ## Build the mesh:
  mesh <- inla.mesh.2d(boundary=boundary,
                       max.edge=c(15, 29.4),
                       min.angle=c(30, 21),
                       max.n=c(48000, 16000), ## Safeguard against large meshes.
                       max.n.strict=c(128000, 128000), ## Don't build a huge mesh!
                       cutoff=3, ## Filter away adjacent points.
                       offset=c(4, 11.9)) ## Offset for extra boundaries, if needed.

  ## Plot the mesh:
  plot(mesh)
  points(loc_xy)
}
