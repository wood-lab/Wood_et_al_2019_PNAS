# Preliminary data analysis of model

source("R/read_dat.R")




# Step 0 -- Load data, create mesh, and assign to list for TMB
data_bulinus <- read_dat(bulinus = T)
data_biomph <- read_dat(bulinus = F)
length(data_biomph$c_q[which(data_biomph$c_q > 0)])
# Histogram of catches

png( file="Report/histograms_of_data.png", width=8, height=10, res=200, units='in')
par(mfrow = c(2,2))
hist(data_bulinus$c_q, breaks = 200, xlab = "Count of Bulinus", xlim = c(0, 50), main = NA)
hist(data_bulinus$y_i, xlab = paste0("0 = not infected (", nrow(data_bulinus$y_i) - sum(data_bulinus$y_i), "), 1 = infected (", sum(data_bulinus$y_i), ")"), ylab = "Count", main = NA, breaks = 2)

hist(data_biomph$c_q, breaks = 30, xlab = "Count of Biomphalaria", xlim = c(0, 50), main = NA)
hist(data_biomph$y_i, xlab = paste0("0 = not infected (", length(data_biomph$y_i) - sum(data_biomph$y_i), "), 1 = infected (", sum(data_biomph$y_i), ")"), ylab = "Count", main = NA, breaks = 2)

dev.off()

