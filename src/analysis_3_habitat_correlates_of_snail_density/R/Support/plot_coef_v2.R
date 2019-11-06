clean_hist <- function( vals , breaks = 20, xlab, xaxt = "s" , xlim ){
  if(missing(xlim)){ xlim = range(vals, na.rm = T)}
  library(FSA)
  par(xpd=F)
  hist_dets <- hist(vals, breaks = breaks, plot = FALSE)
  hist(~vals, main=NA , breaks = breaks, ylab = NA, las= 1, xlab = xlab, ylim = c(0, max(hist_dets$counts) * 1.25), xaxt = xaxt, xlim = xlim, col = "grey")


  par(xpd=T)
  segments(x0 = par("usr")[1], y0= par("usr")[3], x1 = par("usr")[1] , y1 = par("usr")[4] , col =1) # left line
  segments(x0 = par("usr")[2], y0= par("usr")[3], x1 = par("usr")[2], y1 = par("usr")[4] , col =1) # right line
  segments(x0 = par("usr")[1], y0= par("usr")[4] , x1 = par("usr")[2], y1 = par("usr")[4] , col =1) # top line
  segments(x0 = par("usr")[1], y0= par("usr")[3], x1 = par("usr")[2], y1 = par("usr")[3], col =1) # bottom line
  par(xpd=F)
}

# Plot the result of ZIPL coefficients on the response
# Function to do backward model selection using p-values, aic, or bic

plot_coef <- function( rep, Final_Obj, Final_Terms, data_list, data, bulinus, model , fill_cont_mean, fill_cat_mode, BIC_sel, file_name = "NULL", wd = getwd()){
  pred <- Final_Obj$report()

  params_p <- Final_Terms[1, ]
  params_c <- Final_Terms[2, ]

  sig_params <- unique( c(
    colnames(Final_Terms)[which(Final_Terms[1,] != 0)],
    colnames(Final_Terms)[which(Final_Terms[2,] != 0)]))

  # Figure dimensions
  n_panels <- length(sig_params) + 1

  intercepts <- rep$value[which(names(rep$value) == "beta0")][1:2]

  # Figure dimensions
  n_panels <- length(sig_params) + 1

  # 3 X 3
  ncols <- 3
  nrows <- ceiling(n_panels/ncols)
  # Figure object
  for(p in 1:2){
    filename <- paste(wd ,"/Figures/partial_response_curve_", file_name ,c(".eps", ".tiff")[p], sep = "")
    if(p == 1){
      setEPS()
      postscript( file = filename , width=169 / 25.4, height = nrows * 50 / 25.4, family = "serif")
    }

    if(p == 2){
      tiff( file = filename , width=169 / 25.4, height = nrows * 50 / 25.4, family = "serif", units = "in", res = 300)
    }

    # Figure graphics
    nf <- layout(matrix(c(1:(nrows * ncols * 2)) , nrows * 2, ncols , byrow=FALSE), heights = rep(c(1,.6), nrows))
    par( mar=c(0, 3 , 0.5 , 0.25) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))

    # Observed vs predicted
    plot(x = data_list$c_q, y = pred$pred_cq, ylim = c(0, max(data_list$c_q, pred$pred_cq, na.rm = T)), xlim = c(0, max(data_list$c_q, pred$pred_cq, na.rm = T)), xlab = NA, ylab = NA, xaxt = "na", las = 1, pch = 16)
    mtext(side = 2, "Predicted count", line = 2.2, cex = 0.68)
    abline(0,1, lty = 2, col = "grey")
    legend("topleft", "a", bty = "n", inset = -0.035)

    par( mar=c(2.8, 3 , 0 , 0.25) )
    clean_hist(data_list$c_q, breaks = 20, xlab = NA)
    mtext("Observed count", side = 1, line = 1.6, cex = 0.75)

    # Plot covariates
    cat_covars <- grep("cat_covar", sig_params)
    cont_covars <- grep("cont_covar", sig_params)

    for( i in 1:length(sig_params)){
      # Continuous predictor
      if( i %in% cont_covars ){

        # Extract covariate values
        cov_range <- range(data_list$x_q[, sig_params[i]], na.rm = T)
        val <- seq( from = min(cov_range), to = max(cov_range), by = 0.05)

        cov_origin_range <- range( data$density_data[,sig_params[i]] * data$dat_summ[which(data$dat_summ$names == sig_params[i]),"sd_vals"] + data$dat_summ[which(data$dat_summ$names == sig_params[i]),"mean_mode"], na.rm = TRUE)

        # Get partial response
        mu_p <- 1 / (1 + exp( -intercepts[1] + val * params_p[which(names(params_p) == sig_params[i])] ))
        mu_c <- exp( intercepts[2] + val * params_c[which(names(params_c) == sig_params[i])])
        c_hat <- mu_c * mu_p

        # Plot predicted count
        par( mar=c(0, 3 , 0.5 , 0.25))
        plot(y = c_hat, x = val, ylab = NA, xlab = colnames(data_list$x_q)[sig_params[i]], type = "l", xaxt = "na", las = 1)
        legend("topleft", letters[i+1], bty = "n", inset = -0.035)
        mtext(side = 2, "Partial response", line = 2.2, cex = 0.68)

        # Plot histogram
        par( mar=c(2.7, 3 , 0 , 0.25) )
        clean_hist(data_list$x_q[, sig_params[i]], breaks = 20, xlab = NA, xaxt = "na")
        axis( side = 1, at = round(seq( from = min(cov_range), to = max(cov_range), length.out = 5), 1), labels = round(seq( from = min(cov_origin_range), to = max(cov_origin_range), length.out = 5), 1))
        mtext(gsub("cont_covar_", "", sig_params[i]), side = 1, line = 1.6, cex = 0.75)
      }


      # Categorical predictor
      if( i %in% cat_covars){

        # Get values
        val <- c(0, 1)
        mu_p <- 1 / (1 + exp( -intercepts[1] + val * params_p[which(names(params_p) == sig_params[i])] ))
        mu_c <- exp( intercepts[2] + val * params_c[which(names(params_c) == sig_params[i])])
        c_hat <- mu_c * mu_p

        # Plot predicted count
        par( mar=c(0, 3 , 0.5 , 0.25))
        plot(y = c_hat, x = val, ylab = NA, xlab = colnames(data_list$x_q)[sig_params[i]], type = "n", lwd = 6, xaxt = "na", las = 1, xlim = c(-0.25, 1.25))
        mtext(side = 2, "Partial response", line = 2.2, cex = 0.68)

        for(j in 1:length(val)){
          rect(val[j] - 0.12, 0, val[j] + 0.12, c_hat[j])
        }

        legend("topleft", letters[i+1], bty = "n", inset = -0.035)

        # Plot histogram
        par( mar=c(2.7, 3 , 0 , 0.25) )
        clean_hist(data_list$x_q[, sig_params[i]], breaks = 4, xlab = NA, xaxt = "n", xlim = c(-0.15, 1.15))
        mtext(gsub("cat_covar_", "", sig_params[i]), side = 1, line = 1.6, cex = 0.75)
        axis(1, at  = c(0.1,0.9), labels = c(0, 1))
      }
    }

    dev.off()
  }
}
