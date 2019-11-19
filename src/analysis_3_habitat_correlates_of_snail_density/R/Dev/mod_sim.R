# Function to simulate data for the wood-adams delta model

# s = site; q = quadrat; t = sampling period; i = individual snail

# nsites = number of sites sampled
# nquad = number of quadrats within each site sampled
# nt = number of times surveyed
# sigma_c = SD of dispersion
# sigma_sc = SD of random effect of site s from mean count
# sigma_sp = SD of random effect of site s from log-odds of ecounter probability
# mu_c = average log count per quadrat across sites
# mu_p = average log-odds of encounter probability across sites
# mu_psi = average log-odds of probability of infection
# beta_c = impact of local density on log-odds probability of infection
# sigma_psi = SD of random effect of site s from log-odds of infection probability
# sigma_cor = SD of spatial autocorrelation
# scale = ength-scale of the process (how close the points need to be to affect one another), similar to the range
# n_y = number of lat spatial dimensions
# n_x = number of long dimensions


sim_dat <- function(nsites = 32, nquad = 15, nt = 6, sigma_c = 0.01, sigma_sc = 0.02, sigma_sp = 0.1, mu_c = 1, mu_p = 0, mu_psi = -3, beta_c = 0.1, sigma_psi = 0.01, sigma_cor = 1, scale = 2, test = F, n_x = 10, n_y = 10){
  library(RandomFields)
  library(raster)

  if(test == T){
    nsites = 32
    nquad = 15
    nt = 6
    sigma_c = 0.01
    sigma_sc = 0.02
    sigma_sp = 0.1
    mu_c = 1
    mu_p = 0
    mu_psi = -3
    beta_c = 0.1
    sigma_psi = 0.01
    sigma_cor = 1
    scale = 2
    n_x = 10
    n_y = 10
  }

  # Step 0 - Set up spatial array
  dim = c("n_x"=n_x, "n_y"=n_y)
  loc_xy = expand.grid("x"=1:dim['n_x'], "y"=1:dim['n_y'])
  RMmodel = RMgauss(var=sigma_cor, scale=scale)

  # Simulate spatial process for each time step and random sample
  quadrat_t <- array(NA, dim = c(nquad, nt, nsites)) # List of quadrats sampled in each survey
  omega_cq <- array(NA, dim = c(nquad, nt, nsites)) # Spatial autocorrelation for probability of encounter
  omega_pq = array(NA, dim = c(nquad, nt, nsites)) # Spatial autocorrelation for non-zero count
  omega_y = array(NA, dim = c(nquad, nt, nsites)) # Spatial autocorrelation for infection prob

  for(t in 1:nt){
    for(s in 1:nsites){
      quadrat_t[ , t, s] <- sample( c(1:nrow(loc_xy)), nquad, replace = F) # Random sample sites for time t and site s
      omega_cq[ , t, s] = (RFsimulate(model=RMmodel, x=loc_xy[,'x'], y=loc_xy[,'y'])@data[,1])[quadrat_t[ , t, s]] # Spatial autocorrelation for probability of encounter
      omega_pq[ , t, s] = (RFsimulate(model=RMmodel, x=loc_xy[,'x'], y=loc_xy[,'y'])@data[,1])[quadrat_t[ , t, s]] # Spatial autocorrelation for non-zero count
      omega_y[ , t, s] = array(RFsimulate(model=RMmodel, x=loc_xy[,'x'], y=loc_xy[,'y'])@data[,1])[quadrat_t[ , t, s]] # Spatial autocorrelation for infection prob
    }
  }

  # Step 1 - Simulate data for local density
  # ------------------------------------------------------------------
  sim_quad <- data.frame(s_q = rep(1:nsites, each = nquad * nt), q_q = rep(rep(1:nquad, each = nt), nsites ), t_q = rep(rep(1:nt, nquad), nsites))
  sim_quad$id_q <- c(1:nrow(sim_quad)) # Unique identifier for each quadrat sampled

  mat <- matrix(1:nrow(loc_xy), ncol = n_x, nrow = n_y)
  sim_quad$lat_q <- NA
  sim_quad$lon_q <- NA
  for(i in 1:nrow(sim_quad)){
    sim_quad$lat_q[i] <- which(mat == quadrat_t[sim_quad$q_q[i], sim_quad$t_q[i], sim_quad$s_q[i]], arr.ind=TRUE)[1] # y location of sample
    sim_quad$lon_q[i] <- which(mat == quadrat_t[sim_quad$q_q[i], sim_quad$t_q[i], sim_quad$s_q[i]], arr.ind=TRUE)[2] # x location of sample
  }

  # Simulate overdispersion
  gamma_q = rnorm(nrow(sim_quad), mean = 0, sd = sigma_c) # Deviation for each observation away from the predicted count for each quadrat
  epsilon_sc = rnorm(nsites, mean = 0, sd = sigma_sc) # Random deviations for each site away from the mean predicted count
  epsilon_sp = rnorm(nsites, mean = 0, sd = sigma_sp) # Random deviations for each site away from the average probability of 0 snails

  # Get linear predictors
  sim_quad$p_qs <- NA
  sim_quad$c_hat_qs <- NA
  sim_quad$c_q <- NA

  for(i in 1:nrow(sim_quad)){
    sim_quad$p_qs[i] = 1/(1 + exp( - (mu_p + epsilon_sp[sim_quad$s_q[i]] + omega_pq[sim_quad$q_q[i], sim_quad$t_q[i], sim_quad$s_q[i]]))) # Expected probability of zero-encounter for quadrat q in site s and time t
    sim_quad$c_hat_qs[i] = exp(mu_c + gamma_q[i] + epsilon_sc[sim_quad$s_q[i]] + omega_cq[sim_quad$q_q[i], sim_quad$t_q[i], sim_quad$s_q[i]]) # Expected count for quadrat q in site s and time t

    # Predict local density
    sim_quad$c_q[i] = rbinom(1, 1, (1 - sim_quad$p_qs[i])) * rpois(1, sim_quad$c_hat_qs[i])
  }

  # Step 2 - Simulate data for individual infection
  # ------------------------------------------------------------------

  sim_y <- data.frame(s_i = rep(sim_quad$s_q, sim_quad$c_q), q_i = rep(sim_quad$q_q, sim_quad$c_q), id_i = rep(sim_quad$id_q, sim_quad$c_q), t_i = rep(sim_quad$t_q, sim_quad$c_q), c_q = rep(sim_quad$c_q, sim_quad$c_q))
  # s_i = site ID for each individual
  # q_i = quadrat ID for each individual (not unique between sites)
  # id_i = unique ID for each quadrat among sites
  # t_i = time of sample
  # c_q = count of snails in quadrat

  sim_y$lat_q <- NA
  sim_y$lon_q <- NA
  for(i in 1:nrow(sim_y)){
    sim_y$lat_q[i] <- which(mat == quadrat_t[sim_y$q_i[i], sim_y$t_i[i], sim_y$s_i[i]], arr.ind=TRUE)[1] # y location of sample
    sim_y$lon_q[i] <- which(mat == quadrat_t[sim_y$q_i[i], sim_y$t_i[i], sim_y$s_i[i]], arr.ind=TRUE)[2] # x location of sample
  }

  epsilon_psi = rnorm(nsites, mean = 0, sd = sigma_psi) # Random deviations for each site away from the average probability of infection

  sim_y$psi_i <- NA
  sim_y$y_i <- NA

  for(i in 1:nrow(sim_y)){
  sim_y$psi_i[i] = 1/(1 + exp( - (mu_psi + beta_c * sim_y$c_q[i] + epsilon_psi[sim_y$s_i[i]] + omega_y[sim_y$q_i[i], sim_y$t_i[i], sim_y$s_i[i]]))) # Expected probability of infection for quadrat q in site s and individual i

  # Simulate probability of infection
  sim_y$y_i[i] <- rbinom(1, 1, sim_y$psi_i[i])
  }

  dat <- list(sim_quad = sim_quad, sim_y = sim_y)
  return(dat)
}


# Plot dat
plot_dat <- function(site = 1, t = 1, n_x = 10, n_y = 10){

  library(raster)

  data <- sim_dat()

  count <- data$sim_quad
  infection <- data$sim_y

  # Subset
  count <- count[which(count$s_q == site & count$t_q == t),]
  infection <- infection[which(infection$s_i == site & infection$t_i == t),]

  count_mat <- matrix(NA, ncol = n_x, nrow = n_y)
  infection_mat <- matrix(NA, ncol = n_x, nrow = n_y)

  for(i in 1:nrow(count)){
    count_mat[count$lat_q[i], count$lon_q[i]] <- count$c_q[i]
  }

  for(i in 1:nrow(infection)){
  infection_mat[infection$lat_q[i], infection$lon_q[i]] <- 0
  }

  for(i in 1:nrow(infection)){
    infection_mat[infection$lat_q[i], infection$lon_q[i]] <- infection_mat[infection$lat_q[i], infection$lon_q[i]]  + infection$y_i[i]
  }

  png( file="sim_dat.png", width=8, height=4, res=200, units='in')
  par(mfrow = c(1,2))
  plot(raster(count_mat), main = "Number of snails")
  plot(raster(infection_mat), main = "Number of infected snails")
  dev.off()
}

plot_dat()
