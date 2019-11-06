# Plot the result of ZIPL coefficients on the response
# Function to do backward model selection using p-values, aic, or bic

plot_coef <- function( Final_Opt, Final_Obj, Final_Terms, data_list, bulinus, model , fill_cont_mean, fill_cat_mode){


  params_p <- Final_Terms[1, which(Final_Terms[1,] != 0)]
  params_c <- Final_Terms[2, which(Final_Terms[2,] != 0)]

  sig_params <- names(params_c)

  intercepts <- Final_Opt$SD$par.fixed[which(names(Final_Opt$SD$par.fixed) == "beta0")]


  filename <- paste("Report/Figures/partial_response_curve_", c("biom", "bulinus")[as.numeric(bulinus)+1],"model",  model, c("AIC", "BIC")[as.numeric(BIC_sel)+1], c("no_mean_cont_fill", "mean_cont_fill")[as.numeric(fill_cont_mean)+1], c("no_mode_cat_fill", "mode_cat_fill")[as.numeric(fill_cat_mode)+1],".png", sep = "")

  for( i in 1:length(sig_params)){


    png( file = filename , width=8, height=6, res=200, units='in')

    if(length( unique(data_list$x_cq[, sig_params[i]])) > 3 ){
    cov_range <- range(data_list$x_cq[, sig_params[i]], na.rm = T)
    val <- seq( from = min(cov_range), to = max(cov_range), by = 0.05)

    mu_p <- 1 / (1 + exp( -intercepts[1] + val * params_p[sig_params[i]] ))
    mu_c <- exp( intercepts[2] + val * params_c[sig_params[i]])

    c_hat <- mu_c * mu_p

    plot(y = c_hat, x = val, ylab = "Predicted count", xlab = colnames(data_list$x_cq)[sig_params[i]], type = "l")
    }

    if( length( unique(data_list$x_cq[, sig_params[i]])) <= 3){
      val <- c(0, 1)

      mu_p <- 1 / (1 + exp( -intercepts[1] + val * params_p[sig_params[i]] ))
      mu_c <- exp( intercepts[2] + val * params_c[sig_params[i]])

      c_hat <- mu_c * mu_p

      plot(y = c_hat, x = val, ylab = "Predicted count", xlab = colnames(data_list$x_cq)[sig_params[i]], type = "h", lwd = 6, xaxt = "n")
      axis(1, at  = c(0,1))
    }
    dev.off()
  }
}
