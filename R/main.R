
#----------------------------------------------------------------#
#----------------------------------------------------------------#
#              PRESENCE ABSENCE METROPOLIS HASTINGS              #
#----------------------------------------------------------------#
#----------------------------------------------------------------#

# The following commands ensure that package dependencies are listed in the
# NAMESPACE file.

#' @import fftwtools
#' @useDynLib silverblaze
#' @importFrom Rcpp evalCpp
#' @import graphics
#' @import stats
#' @import utils
NULL

#-------------------------------------------------------------------
#' @title PAjoint
#'
#' @description Calculate the logged joint probability of observing the data.
#'   This acts as the posterior draw. Also extract the logged density of the
#'   bivariate normal at each trap location.
#'
#' @param theta_x	The mean longitude for the source we condition on
#' @param theta_y	The mean latitude for the source we condition on
#' @param sigma The standard deviation of the bivariate normal (source) we
#'   condition on
#' @param data The trap data consisting of locations in lat/long and associated
#'   densities
#' @param trap_rad The trap radius we use to approximate the area under the
#'   bivariate normal (in kilometres)
#' @param exp_pop The expected population over our entire area
#' @param priorLat The mean Latitude of the prior on source locations - normally
#'   the mean of hit locations
#' @param priorLon The mean Longitude of the prior on source locations -
#'   normally the mean of hit locations
#' @param prior Specify which parameter a prior is being set for
#' @param meanSigprior The mean of the prior on sigma
#' @param priorSD The standard devitation for any prior all are set to normals
#'
#' @export
#' @examples
#' trap_sim <- cbind(runif(10), runif(10), sample(0:5, 10, replace = TRUE))
#' PAjoint(theta_x = 0, theta_y = 0, sigma = 10, data = trap_sim, trap_rad = 2, exp_pop = 25,
#' priorLat = 0.5, priorLon = 0.5, prior = "sourcePrior"
#' , meanSigprior = NULL, priorSD = 10)

PAjoint <- function(theta_x, theta_y, sigma, data, trap_rad, exp_pop, priorLat, priorLon, prior = "None", meanSigprior, priorSD) {
  
  # extract useful parameters
  K <- length(theta_x)  # K clusters
  n <- nrow(data)       # n traps
  
	# sum poisson params across sources with EQUAL waitings dictated by log(K)
	heights <- mapply(function(x, y) {
	    bvnorm(data[,1], data[,2], x, y, sd = sigma, log = TRUE)
	  }, x = theta_x, y = theta_y)
	z <- log(pi) + 2*log(trap_rad) + apply(heights, 1, logSum) - log(K)
  
	# log-likelihood
  log_lambda <- log(exp_pop) + z
  likeli_val  <- sum( data[,3]*log_lambda - exp(log_lambda) - lfactorial(data[,3]) )
	if (prior == "sourcePrior") {
		prior_val <- sum(bvnorm(x = theta_x, y = theta_y, mu_x = priorLon, mu_y = priorLat, sd = priorSD, log = TRUE))
	} else if (prior == "sigPrior") {
    prior_val <- dnorm(sigma, mean = meanSigprior, sd = priorSD, log = TRUE)
  } else if (prior == "None") {
    prior_val <- 0
	}
  posterior_val <- prior_val + likeli_val
  
  # return list
  return(list(posterior_val = posterior_val, z = z))
}

#------------------------------------------------
#' @title Propose parameter values
#'
#' @description Propose new parameter values (specifically a source location)
#'   from a bivaraite normal.
#'
#' @param n The number of new points to propose
#' @param mu1 The centre of the bivariate norm to draw from (longitude)
#' @param mu2 The centre of the bivariate norm to draw from (latitude)
#' @param sd The standard deviation of the bivariate norm to draw from
#'   (kilometres)
#'
#' @export
#' @examples
#' proposal(n = 16, mu1 = 0, mu2 = 0, sd = 5)

proposal <- function(n, mu1, mu2, sd) {
  points <- rnorm_sphere(n = n, centre_lon = mu1, centre_lat = mu2, sigma = sd)
  cbind(longitude = points$longitude, latitude = points$latitude)
}

#------------------------------------------------
#' @title PAmcmc
#'
#' @description Run the Metropolis-Hastings algorithm with Gibbs sampler to fit 
#'   the source locations, sigma and or population density
#'
#' @param trapData The main input into the model, consisting of a matrix or data
#'   frame with three columns: trap longitude, latitude and counts
#' @param burnin The MCMC's burn in; the number of iterations before the main 
#'   sampling begins
#' @param iterations The number of iterations in MCMC, this is on top of the 
#'   burn in
#' @param clusters The number of clusters to be searched for
#' @param trap_rad The trap radius in km
#' @param s_prior_shape The shape parameter for the prior on the population 
#'   density
#' @param s_prior_rate The rate parameter for the prior on the population 
#'   density
#' @param priorLon The prior mean for the source locations (Longitude)
#' @param priorLat The prior mean for the source locations (Latitude). Uses mean
#'   of data by default.
#' @param lon_minMax The longitudinal limits of the data (already defined)
#' @param lat_minMax The latitudinal limits of the data (already defined)
#' @param produce_surface Option to produce smoothed surface from posterior
#'   draws
#'
#' @export
#' @examples
#' trap_sim <- cbind(runif(10), runif(10), sample(0:5, 10, replace = TRUE))
#' PAmcmc(trapData = trap_sim, burnin = 100, iterations = 25e2, clusters = 1, trap_rad = 1,
#' s_prior_shape = 1, s_prior_rate = 1)

PAmcmc <- function(trapData, burnin = 100, iterations = 1e3, clusters = 1, trap_rad = 1, s_prior_shape = 1, s_prior_rate = 1, priorLon = mean(trapData[,1]), priorLat = mean(trapData[,2]), lon_minMax = range(trapData[,1]), lat_minMax = range(trapData[,2]), produce_surface = TRUE) {
  
  # start timer
  start <- Sys.time()
  
  # process data
  trapData <- as.matrix(trapData)
  trap_density <- sum(trapData[,3])
  hit_data <- subset(trapData, trapData[,3]>0)
  
  # define sigma prior mean based on minimum distance between points
  d <- dist_gc(hit_data[,1:2])
  diag(d) <- NA
  d[d==0] <- NA
  distance_min <- apply(d, 1, min, na.rm = TRUE)
  sigma_prior_mean <- 0.5*mean(distance_min, na.rm = TRUE)
  
  # get diagonal distance from corner to corner of domain
  areaKm <- latlon_to_bearing(origin_lat = lat_minMax[1], origin_lon = lon_minMax[1], dest_lat = lat_minMax[2], dest_lon = lon_minMax[2])$gc_dist
  
  # define MCMC parameters
	# Robbins-Monro Step Constant
  RMcM <- 5
  RMcS <- 5.579
  
	# initial proposal standard deviations
  propSD_mu <- rep(1, clusters) #rep(areaKm, clusters)  # TODO - swap back
  propSD_sig <- 1 #10*areaKm  # TODO - swap back
  
	### propose starting values
  # source locations
	old_theta <- proposal(clusters, mu1 = 0, mu2 = 0, sd = 2*areaKm)
  colnames(old_theta) <- c("Long", "Lat")
  
  # sigma
  old_sigma <- 2 #abs(rnorm(1, mean = sigma_prior_mean, sd = propSD_sig)) # TODO - swap back
  
	# population density
	new_pop_den <- 1 #rgamma(1, shape = s_prior_shape, rate = s_prior_rate)  # TODO - swap back
  
	# calculate intial joint probability
  old_posterior <- PAjoint(theta_x = old_theta[,1], theta_y = old_theta[,2], sigma = old_sigma, data = trapData, trap_rad = trap_rad, exp_pop = new_pop_den, priorLat = priorLat, priorLon = priorLon, prior = "sourcePrior", meanSigprior = sigma_prior_mean, priorSD = areaKm)
  
  ### create objects for storing results
  joint_lat <- joint_lon <- matrix(NA, ncol = clusters, nrow = iterations)
  joint_sig <- rep(NA, iterations)
  pop_densities <- rep(NA, iterations)
  a_r_rate_theta <- matrix(0, ncol = clusters, nrow = iterations)
  a_r_rate_sig <- rep(0, iterations)
  propSD_mu_vals <- matrix(NA, ncol = clusters, nrow = iterations)
  propSD_sig_vals <- rep(NA, iterations)
  
  # run MCMC
  for (i in 1:iterations) {
    
    # report current iteration
    if ((i %% 500)==0) {
      cat(paste("  iteration:", i, "\n"))
    }
    
    # loop through clusters
    for (j in 1:clusters) {
      
      # propose new theta, accept or reject. proposal() returns two values, a Lon and Lat
      new_theta <- old_theta
      new_theta[j,] <- proposal(1, mu1 = old_theta[j,1], mu2 = old_theta[j,2], sd = propSD_mu[j])
      
      # calculate new joint probability
      new_posterior <- PAjoint(theta_x = new_theta[,1], theta_y = new_theta[,2], sigma = old_sigma, data = trapData, trap_rad = trap_rad, exp_pop = new_pop_den, priorLat = priorLat, priorLon = priorLon, prior = "sourcePrior", meanSigprior = sigma_prior_mean, priorSD = areaKm)
      
      # Metropolis-Hastings step
      ratio <- new_posterior$posterior_val - old_posterior$posterior_val
      if (log(runif(1)) < ratio) {
        old_theta[j,] <- new_theta[j,]
        old_posterior <- new_posterior
        a_r_rate_theta[i, j] <- 1
      }
      if (i <= burnin) {
        theta_shift <- ((RMcM*0.766)^a_r_rate_theta[i,j])*(-RMcM*0.234)^(1 - a_r_rate_theta[i,j])/sqrt(i)
        
        # Robbins Monro step for source location proposal
        propSD_mu[j] <- propSD_mu[j] + theta_shift
        propSD_mu[j] <- max(1e-6, propSD_mu[j])
      }
      
    }  # end loop over clusters
    
    if (FALSE) {  # TODO - remove
    # propose new sigma, accept or reject
    new_sigma <- abs(rnorm(1, mean = old_sigma, sd = propSD_sig))
    
    # calculate new joint probability
    new_posterior <- PAjoint(theta_x = old_theta[,1], theta_y = old_theta[,2], sigma = new_sigma, data = trapData, trap_rad = trap_rad, exp_pop = new_pop_den, priorLat = priorLat, priorLon = priorLon, prior = "sigPrior", meanSigprior = sigma_prior_mean, priorSD = areaKm)
    
    # Metropolis-Hastings step
    ratio <- new_posterior$posterior_val - old_posterior$posterior_val
    if (log(runif(1)) < ratio) {
      old_sigma <- new_sigma
      old_posterior <- new_posterior
      a_r_rate_sig[i] <- 1
    }
    if (i <= burnin) {
      sig_shift <- ( (RMcS*0.766^(a_r_rate_sig[i]))*(-RMcS*0.234)^(1 - a_r_rate_sig[i])  )/sqrt(i)
      propSD_sig <- propSD_sig + sig_shift    # Robbins Monro step for sigma proposal
      propSD_sig <- max(1e-6, propSD_sig)
    }
    }  # end temporary FALSE statement
    
    #update s by Gibbs sampling
    #new_z <- PAjoint(theta_x = old_theta[,1], theta_y = old_theta[,2], sigma = old_sigma, data = trapData, trap_rad = trap_rad, exp_pop = new_pop_den, priorLat = priorLat, priorLon = priorLon, prior = "None", meanSigprior = sigma_prior_mean)$z
    #new_pop_den <- rgamma(1, shape = s_prior_shape + trap_density, rate = s_prior_rate + exp(logSum(new_z)))
    
    # store values of this iteration
    joint_lon[i,] <- old_theta[,1]
    joint_lat[i,] <- old_theta[,2]
    joint_sig[i] <- old_sigma
    pop_densities[i] <- new_pop_den
    propSD_mu_vals[i,] <- propSD_mu
    propSD_sig_vals[i] <- propSD_sig
  }
  
  # trim burn-in
  lon_draws_burnin <- joint_lon[1:burnin,]
  lat_draws_burnin <- joint_lat[1:burnin,]
  lon_draws_samples <- joint_lon[-(1:burnin),]
  lat_draws_samples <- joint_lat[-(1:burnin),]
  
  # produce smooth surface using geoSmooth()
  rawSurface <- lon_seq <- lat_seq <- NULL
  if (produce_surface) {
    lon_seq <- seq(lon_minMax[1], lon_minMax[2], (lon_minMax[2] - lon_minMax[1])/500)
    lat_seq <- seq(lat_minMax[1], lat_minMax[2], (lat_minMax[2] - lat_minMax[1])/500)
    rawSurface <- geoSmooth(lon_draws_samples, lat_draws_samples, breaks_lon = lon_seq, breaks_lat = lat_seq, lambda = NULL)
    cat("\n")
  }
  
  # end timer
  end <- Sys.time()
  time_elapsed <- end - start
  print(time_elapsed)
  
  # return as list
  ret <- list(lon_draws_burnin = lon_draws_burnin,
              lat_draws_burnin = lat_draws_burnin,
              lon_draws_samples = lon_draws_samples,
              lat_draws_samples = lat_draws_samples,
              propSD_mu_vals = propSD_mu_vals,
              a_r_rate_theta = a_r_rate_theta,
              lon_seq  = lon_seq,
              lat_seq = lat_seq,
              joint_sig = joint_sig,
              a_r_rate_sig = a_r_rate_sig,
              propSDsig_vals = propSD_sig_vals,
              pop_densities = pop_densities,
              rawSurface = rawSurface,
              trap_rad = trap_rad,
              Clusters = clusters,
              popPriorSh = s_prior_shape ,
              popPriorRa = s_prior_rate,
              sig_init = sigma_prior_mean,
              areaKm = areaKm)
  return(ret)
}
