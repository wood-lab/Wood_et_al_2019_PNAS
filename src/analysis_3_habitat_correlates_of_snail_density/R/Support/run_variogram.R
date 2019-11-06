# function to create spatial variograms for large and small scale processes based on the residuals of a model. Currently only for count data. Takes in the data and the predicted response as inputs

run_variogram <- function( data, data_list, pred, bulinus, model, wd, file_name = "NULL"){

  if(bulinus == T){
    c_q <- "snail_is_trunc_globo"
  }
  if(bulinus == F){
    c_q <- "snail_is_biomph"
  }

  # Build dataset of interest
  dat <- data$density_data
  dat[,c_q] <- as.numeric(dat[,c_q])
  dat$c_q <- data$c_q
  dat$pred <- pred$pred_cq
  dat$resid <- dat[,c_q] - dat$pred
  dat$lat <- as.numeric(dat$lat)
  dat$long <- as.numeric(dat$long)
  dat <- dat[which(!is.na(dat$lat) & !is.na(dat$long)),]

  # Fit small scale spati-temporal variograms
  site_time <- unique(dat[, c("site", "field_mission" )])
  for(i in  1:nrow(site_time)){
    dat_sub <- dat[which(dat$site == site_time$site[i] & dat$field_mission == site_time$field_mission[i]), ]
    if(sum(is.na(dat_sub$lat)) == 0){
      sp::coordinates(dat_sub) <- c("long", "lat")
      v = gstat::variogram(resid~1, dat_sub)
      dat_sub <- as.data.frame(dat_sub)

      filename <- paste(wd, "/Figures/Variograms/spatio_temporal_variogram_of_ ", file_name,".png",sep = "")
      png( file=filename, width=12, height=6, res=200, units='in')
      par(mfrow = c(1,2))
      plot(x = v$dist, y = v$gamma, main = paste0("Variogram of residuals: ", dat_sub$site[1]), xlab = "Distance (m)", ylab = "Semivariance")
      plot(x = dat_sub[,c_q], y = dat_sub$pred, main = paste0("Fitted vs observed:", dat_sub$site[1]), xlab = "Observed count", ylab = "Predicted count")
      dev.off()
    }
  }

  # Fit small scale spatial variograms
  site_time <- unique(dat[, c("site" )])
  for(i in  1:length(site_time)){
    dat_sub <- dat[which(dat$site == site_time[i]), ]
    if(sum(is.na(dat_sub$lat)) == 0){
      sp::coordinates(dat_sub) <- c("long", "lat")
      v = gstat::variogram(resid~1, dat_sub)
      dat_sub <- as.data.frame(dat_sub)

      filename <- paste(wd, "/Figures/Variograms/spatial_variogram_of_ ", file_name,".png",sep = "")
      png( file=filename, width=12, height=6, res=200, units='in')
      par(mfrow = c(1,2))
      plot(x = v$dist, y = v$gamma, main = paste0("Variogram of residuals: ", dat_sub$site[1]), xlab = "Distance (m)", ylab = "Semivariance")
      plot(x = dat_sub[,c_q], y = dat_sub$pred, main = paste0("Fitted vs observed:", dat_sub$site[1]), xlab = "Observed count", ylab = "Predicted count")
      dev.off()
    }
  }


  # Do variogram across sites
  mean_coords <- aggregate(dat[,c("lat", "long")], list(dat$site), mean) # get mean lat and lon
  names(mean_coords) <- c("site", "site_lat", "site_long")
  dat <- merge(dat, mean_coords, by = "site")

  sp::coordinates(dat) <- c("site_long", "site_lat")
  v = gstat::variogram(resid~1, dat)#, cutoff = max(spDists(dat))/2)
  dat <- as.data.frame(dat)
  plot(v)
  filename <- paste(wd, "/Figures/Variograms/among_site_variogram_", file_name,  ".png",sep = "")
  png( file=filename, width=12, height=6, res=200, units='in')
  par(mfrow = c(1,2))
  plot(x = v$dist, y = v$gamma, main = paste0("Variogram of residuals: among site"), xlab = "Distance (m)", ylab = "Semivariance")
  plot(x = dat[,c_q], y = dat$pred, main = paste0("Fitted vs observed: among site"), xlab = "Observed count", ylab = "Predicted count")
  dev.off()

}

