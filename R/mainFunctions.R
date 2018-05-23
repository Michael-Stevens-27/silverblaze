#----------------------------------------------------------------#
#----------------------------------------------------------------#
#              PRESENCE ABSENCE METROPOLIS HASTINGS              #
#----------------------------------------------------------------#
#----------------------------------------------------------------#

# The following commands ensure that package dependencies are listed in the NAMESPACE file.

#' @useDynLib silverblaze
#' @importFrom Rcpp evalCpp
#' @importFrom plot3D hist3D
NULL

#------------------------------------------------
#' Dummy function
#'
#' This is a dummy function
#'
#' @param x Some parameter
#'
#' @export
#' @examples
#' dummy1()

dummy1 <- function() {
	print("dummy1")
}

#------------------------------------------------
#' logSum
#'
#' To avoid underflow issues from R, scale the product of small values by xmax.
#'
#' @param x A vector whose values are in log space
#'
#' @export
#' @examples
#' x <- seq(5e-49, 5e-48, 1e-49)
#' log(prod(x))
#' logSum(x)

logSum <- function(x){
  xmax <- max(x)
  xmax + log(sum(exp(x-xmax)))
}

#------------------------------------------------
#' latlon_to_cartesian
#'
#' Convert from latitude/longitude to cartesian co-ordinates
#'
#' @param centre_lat The centre latitude
#' @param centre_lon The centre longitude
#' @param data_lat   The data latitude
#' @param data_lon   The data longitude
#'
#' @export
#' @examples
#' centre_latLon <- c(0,1)
#' data_latLon <- c(0,1)
#' latlon_to_cartesian(centre_lat = centre_latLon[1], centre_lon = centre_latLon[2], data_lat = data_latLon[1], data_lon = data_latLon[2])
#  From B. VERITY

latlon_to_cartesian <- function(centre_lat, centre_lon, data_lat, data_lon) {

  # calculate bearing and great circle distance of data relative to centre
  data_trans <- latlon_to_bearing(centre_lat, centre_lon, data_lat, data_lon)

  # use bearing and distance to calculate cartesian coordinates
  theta <- data_trans$bearing*2*pi/360
  d <- data_trans$gc_dist
  data_x <- d*sin(theta)
  data_y <- d*cos(theta)

  return(list(x=data_x, y=data_y))
}

#------------------------------------------------
#' latlon_to_bearing
#'
#' Calculate distance between latitude/longitude points
#'
#' @param origin_lat The origin latitude
#' @param origin_lon The origin longitude
#' @param dest_lat   The destination latitude
#' @param dest_lon   The destination longitude
#'
#' @export
#' @examples
#' origin_latLon <- c(0,0)
#' dest_latLon <- c(1,1)
#' latlon_to_bearing(origin_lat = origin_latLon[1], origin_lon = origin_latLon[2], dest_lat = dest_latLon[1], dest_lon = dest_latLon[2])
#  From B. VERITY

latlon_to_bearing <- function(origin_lat, origin_lon, dest_lat, dest_lon) {
  # convert input arguments to radians
  origin_lat <- origin_lat*2*pi/360
  dest_lat <- dest_lat*2*pi/360
  origin_lon <- origin_lon*2*pi/360
  dest_lon <- dest_lon*2*pi/360

  delta_lon <- dest_lon-origin_lon

  # calculate bearing and great circle distance
  bearing <- atan2(sin(delta_lon)*cos(dest_lat), cos(origin_lat)*sin(dest_lat)-sin(origin_lat)*cos(dest_lat)*cos(delta_lon))
  gc_angle <- acos(sin(origin_lat)*sin(dest_lat) + cos(origin_lat)*cos(dest_lat)*cos(delta_lon))
  gc_angle[is.nan(gc_angle)] <- 0
  # convert bearing from radians to degrees measured clockwise from due north, and convert gc_angle to great circle distance via radius of earth (km)
  bearing <- bearing*360/(2*pi)
  bearing <- (bearing+360)%%360
  earthRad <- 6371
  gc_dist <- earthRad*gc_angle

  return(list(bearing=bearing, gc_dist=gc_dist))
}

#
# # # install_github("Michael-Stevens-27/silverblaze", ref = "master")
# # # library(silverblaze)
# #
# # # rm(list = ls()) #remove all objects
# #
# # # set.seed(5) # throws out ratio error
# # Draw from normal distribution converted to spherical coordinate system. Points are first drawn from an ordinary cartesian 2D normal distribution.
# # The distances to points are then assumed to be great circle distances, and are combined with a random bearing from the point {centre_lat, centre_lon}
# # to produce a final set of lat/lon points. Note that this is not a truly spherical normal distribution, as the domain of the distribution is not the sphere -
# # rather it is a transformation from one coordinate system to another that is satisfactory when the curvature of the sphere is not severe.
# # From B. VERITY
#
# rnorm_sphere <- function(n, centre_lat, centre_lon, sigma) {
# 	x <- rnorm(n,sd=sigma)
# 	y <- rnorm(n,sd=sigma)
# 	output <- cartesian_to_latlon(centre_lat, centre_lon, x, y)
# 	return(output)
# }
#
# # Calculate destination lat/lon given an origin, a bearing and a great circle distance of travel
# # Note that bearing should be in degrees relative to due north, and gc_dist should be in units of kilometres
# # From B. VERITY
#
# bearing_to_latlon <- function(origin_lat, origin_lon, bearing, gc_dist) {
#
#   # convert origin_lat, origin_lon and bearing from degrees to radians
#   origin_lat <- origin_lat*2*pi/360
#   origin_lon <- origin_lon*2*pi/360
#   bearing <- bearing*2*pi/360
#
#   # calculate new lat/lon using great circle distance
#   earthRad <- 6371
#   new_lat <- asin(sin(origin_lat)*cos(gc_dist/earthRad) + cos(origin_lat)*sin(gc_dist/earthRad)*cos(bearing))
#   new_lon <- origin_lon + atan2(sin(bearing)*sin(gc_dist/earthRad)*cos(origin_lat), cos(gc_dist/earthRad)-sin(origin_lat)*sin(new_lat))
#
#   # convert new_lat and new_lon from radians to degrees
#   new_lat <- new_lat*360/(2*pi)
#   new_lon <- new_lon*360/(2*pi)
#
#   return(list(longitude=new_lon, latitude=new_lat))
# }
#
# # Convert Cartesian co-ordinates to lat long
# # From B. VERITY
#
# cartesian_to_latlon <- function(centre_lat, centre_lon, data_x, data_y) {
#
#   # calculate angle and euclidian distance of all points relative to origin
#   d <- sqrt(data_x^2+data_y^2)
#   theta <- atan2(data_y,data_x)
#
#   # convert theta to bearing relative to due north
#   theta <- theta*360/(2*pi)
#   theta <- (90-theta)%%360
#
#   # use bearing and great circle distance to calculate lat/lon relative to an origin point
#   data_trans <- bearing_to_latlon(centre_lat, centre_lon, theta, d)
#
#   return(list(longitude=data_trans$longitude, latitude=data_trans$latitude))
# }
#
# # calculate pairwise distance
# # From B. VERITY
#
# pairwise_distance <- function(points){
# 	distance <- matrix(NA, ncol = length(points[,1]), nrow = length(points[,1]))
#
# 	for(i in 1:length(points[,1]))
# 	{
# 	for(j in 1:i)
# 		{
# 		x1 <- points[i, 1]
# 		y1 <- points[i, 2]
# 		x2 <- points[j, 1]
# 		y2 <- points[j, 2]
# 		dist <- latlon_to_bearing(y1, x1, y2, x2)$gc_dist
# 		distance[i,j] <- dist
# 		}
# 	}
# 	diag(distance)  <- NA
# 	distance[distance==0] <- NA
# 	distance_min <- apply(distance, 1, min, na.rm = TRUE)
# 	distance_min[distance_min =="Inf"] <- NA
# 	SD_TR <- list(distance = distance, distance_min = distance_min)
# 	return(SD_TR)
# }

# #-------------------------------------------------------------------
# # Presence Absence Functions
#
# # Extract the height on the bivariate normal centred at (mu_x, mu_y)
# Prob_data <- function(x, y, mu_x, mu_y, Sd) {
#   param <- dnorm(latlon_to_bearing(y, x, y, mu_x)$gc_dist, mean=0, sd=Sd, log=TRUE) + dnorm(latlon_to_bearing(y, mu_x, mu_y, mu_x)$gc_dist, mean=0, sd=Sd, log=TRUE)
#   return(param)
# }
#
# # x is in log space, return log(sum(exp(x))) to reduce underflow issues
# logSum <- function(x){
#   xmax <- max(x)
#   xmax + log(sum(exp(x-xmax)))
# }
#
#
# # joint probability
# post <- function(theta_x, theta_y, prior_Sd = NULL, sigma, data, trap_rad, exp_pop) {
#   # extract useful parameters
#   K <- length(theta_x)  # K clusters
#   n <- nrow(data)       # n traps
# 
#   # distinguish between a varying or constant sigma
#   # if(length(sigma) == 1)
#   # {
#   #   sigma <- rep(sigma, K)
#   # }
#
#   # # log-likelihood - V.S.
#   # heights <- matrix(NA, nrow = n, ncol = K)
#   # for (k in 1:K) {
#   #   heights[,k] <- dnorm(latlon_to_bearing(data[,2], data[,1], data[,2], theta_x[k])$gc_dist, mean=0, sd=sigma[k], log=TRUE) + dnorm(latlon_to_bearing(data[,2], data[,1], theta_y[k], data[,1])$gc_dist, mean=0, sd=sigma[k], log=TRUE)
#   # }
#   # log-likelihood
#   heights <- matrix(NA, nrow = n, ncol = K)
#   for (k in 1:K) {
#     heights[,k] <- dnorm(latlon_to_bearing(data[,2], data[,1], data[,2], theta_x[k])$gc_dist, mean=0, sd=sigma, log=TRUE) + dnorm(latlon_to_bearing(data[,2], data[,1], theta_y[k], data[,1])$gc_dist, mean=0, sd=sigma, log=TRUE)
#   }
#   # Sum poisson params across sources with EQUAL waitings dictated by log(K)
#   z <- log(pi) + 2*log(trap_rad) + apply(heights, 1, logSum) - log(K)
#   log_lambda <- log(exp_pop) + z
#   likeli_val  <- sum( data[,3]*log_lambda - exp(log_lambda) - lfactorial(data[,3]) )
#   # posterior_val <- prior_val + likeli_val # prior currently undefined
#   posterior_val <- likeli_val
#   return(list(posterior_val = posterior_val, z = z))
# }
#
# # draw from proposal
# proposal <- function(n, mu1, mu2, Sd) {
#   points <- rnorm_sphere(n = n, centre_lon = mu1, centre_lat = mu2, sigma = Sd)
#   cbind(longitude = points$longitude, latitude = points$latitude)
# }
#
# # A function that simulates and collects all observed and unobserved data
# PAsim <- function(n_true = 100, sigma_true = 1, long_minMax = c(0, 0.5), lat_minMax = c(0, 0.5), n_traps = 100, trap_const = 0.3, K_true = 3, single_count = F, plotting = F, plotRail = 0.05, trap_spacing = "random", trap_clusters = NULL, bias = 0)
#   {
#   C <- 1
#   while(C < 2)
#   {
#   # define search/area to generate data within
# 	total_area <- 1
#   trap_rad_true <- trap_const*sigma_true
#   fakeSlong <- NULL
#   fakeSlat <- NULL
#
# 	# Draw source locations from a uniform prior (i.e randomly)
#   if(bias == 1)
#   {
#   source_lat <- runif(K_true - 2, lat_minMax[1], 0.5*(lat_minMax[2]- lat_minMax[1]) )
#   source_long <- runif(K_true - 2, long_minMax[1], 0.5*(long_minMax[2]- long_minMax[1] ) )
#   fakeSlong <- runif(1, 0.5*(long_minMax[2]- long_minMax[1]), long_minMax[2])
#   fakeSlat <- runif(1, lat_minMax[1], 0.5*(lat_minMax[2]- lat_minMax[1]))
#   source_lat <- c(source_lat, runif(1, 0.5*(lat_minMax[2]- lat_minMax[1]), lat_minMax[2]), fakeSlat)
#   source_long <- c(source_long, runif(1, long_minMax[1], 0.5*(long_minMax[2]- long_minMax[1])), fakeSlong)
#   # Draw number of observations from a Poisson with rate n_true * total_area and allocate to sources
#   n_obs <- rpois(1, n_true*total_area)
#   alloc <- c(sample(1:(K_true - 1), size = n_obs, replace = T), 0)
#     } else{
#     source_lat <- runif(K_true, lat_minMax[1], lat_minMax[2])
#     source_long <- runif(K_true, long_minMax[1], long_minMax[2])
#     n_obs <- rpois(1, n_true*total_area)
#     alloc <- c(sample(1:(K_true), size = n_obs, replace = T), 0)
#   }
#   source_loc <- data.frame(source_long = source_long, source_lat = source_lat)
#   perSource <-  mapply(function(x){length(which(x == alloc))}, x = 1:K_true)
#
#   # Distinguish between constant or variable sigma
#   # if(length(sigma_true) == 1)
#   # {
#   #   sigma_true <- rep(sigma_true, K_true)
#   # }
#
#   # Draw observation locations from a bivariate norm with standard deviation sigma_true NEEDS EDITING
#   crime_long <- unlist(mapply(rnorm_sphere, n = perSource, centre_lat = source_lat, centre_lon = source_long, sigma = sigma_true)[1,])
#   crime_lat <- unlist(mapply(rnorm_sphere, n = perSource, centre_lat = source_lat, centre_lon = source_long, sigma = sigma_true)[2,])
# 	crime_loc <- data.frame(longitude = crime_long, latitude = crime_lat)
#
#   # set trap spacing to random, uniform, or clustered/uniform (trap_clusters dictates the number of traps in each uniform location)
#   if(trap_spacing == "random")
#   {
#     trap_lat <- runif(n_traps, lat_minMax[1], lat_minMax[2] - (lat_minMax[2]/2 - lat_minMax[1]/2)*bias)
#     #trap_lon <- runif(n_traps, long_minMax[1], long_minMax[2] - (long_minMax[2]/2 - long_minMax[1]/2)*bias)
#     trap_lon <- runif(n_traps, long_minMax[1], long_minMax[2])
#     trap_loc <- data.frame(trap_lon = trap_lon, trap_lat = trap_lat)
#
#   } else if(trap_spacing == "uniform")
#   {
#     trap_lat <- seq(lat_minMax[1], lat_minMax[2]/(2^bias), (lat_minMax[2]/(2^bias) - lat_minMax[1])/sqrt(n_traps))
#     trap_lon <- seq(long_minMax[1], long_minMax[2], (long_minMax[2] - long_minMax[1])/sqrt(n_traps))
#     trap_loc <- expand.grid(trap_lon, trap_lat)
#
#   } else if(trap_spacing == "cluster") {
#
#     trap_lat <- seq(lat_minMax[1], lat_minMax[2]/(2^bias), (lat_minMax[2]/(2^bias) - lat_minMax[1])/sqrt(n_traps))
#     trap_lon <- seq(long_minMax[1], long_minMax[2], (long_minMax[2] - long_minMax[1])/sqrt(n_traps))
#     trap_loc <- expand.grid(trap_lon, trap_lat)
#     c_trap <- seq(0, 2*pi, 2*pi/trap_clusters)
#     # Inter cluster distance between traps needs specifying
# 		polx <- 2*sin(c_trap)
# 		poly <- 2*cos(c_trap)
# 		trap_lon <- mapply(function(x, y){cartesian_to_latlon(y, x, poly, polx)$longitude}, x = trap_loc[,1], y = trap_loc[,2])
# 		trap_lat <- mapply(function(x, y){cartesian_to_latlon(y, x, poly, polx)$latitude}, x = trap_loc[,1], y = trap_loc[,2])
#     trap_loc <- cbind(c(trap_lon), c(trap_lat))
#
#   } else { stop("trap_spacing must be set to 'uniform','random' or 'cluster'")}
# 	# allocate raw data to traps
#
#   all_dist <- mapply(function(x, y){latlon_to_bearing(trap_loc[,2], trap_loc[,1], x, y)$gc_dist}, x = crime_loc$latitude, y = crime_loc$longitude)
#   all_dist <- t(all_dist) # rows are crimes, columns are traps, entries are distances between crimes and traps
#   trapCounts <- colSums(all_dist<=trap_rad_true)  # concatonate crimes into traps given the distance is within trap radius
#   is_observed <- rowSums(all_dist<=trap_rad_true)  # count the number of observed observations
#   if(length(which(trapCounts > 0)) > 4) # check there are more than 5 observations overall
#   #if(TRUE)
#   {
#   #allow for single or multiple trappings of the same observation
# 	if(single_count == T)
# 	{
#   if(max(is_observed) == 1) # check if no trap contains more than one observation
#   {} else if(length(which(is_observed > 1)) > 1) { # check that there is more than one trap with multiple observations
#     single_densities <- table(apply(all_dist[which(is_observed > 1),], FUN = which.min, MARGIN = 1))
#     index <- strtoi(names(single_densities))
#     trapCounts[index] <- single_densities
#   } else { # if there is a single trap with multiple observations, then replace with a single obs
#   single_index <- which.min(all_dist[which(is_observed > 1),])
#   trapCounts[single_index] <- 1
#   }
# 	}
#
# 	# make final trap data
# 	trap_data <- data.frame(longitude=trap_loc[,1], latitude=trap_loc[,2], count=trapCounts)
#   hit_data <- subset(trap_data, trap_data$count >0 )
#   miss_data <- subset(trap_data, trap_data$count ==0 )
#
#   source_hits_dist <- t(mapply(function(x, y){latlon_to_bearing(hit_data[,2], hit_data[,1], x, y)$gc_dist}, x = source_loc$source_lat, y = source_loc$source_long)) # rows are sources, columns are hits, entries are distances between sources and hits
#   nega_source1 <- rowSums(source_hits_dist <= sigma_true)
#   nega_source2 <- rowSums(source_hits_dist <= 2*sigma_true)
#   nega_source3 <- rowSums(source_hits_dist <= 3*sigma_true)
#
# 	# plot if plotting returns true
# 	if(plotting == T)
# 	{
# 		c_seq <- seq(0, 2*pi, 0.25)
# 		datax <- trap_rad_true*sin(c_seq)
# 		datay <- trap_rad_true*cos(c_seq)
# 		circle_longs <- mapply(function(x, y){cartesian_to_latlon(y, x, datay, datax)$longitude}, x = trap_data$longitude, y = trap_data$latitude)
# 		circle_lats <- mapply(function(x, y){cartesian_to_latlon(y, x, datay, datax)$latitude}, x = trap_data$longitude, y = trap_data$latitude)
#
# 		plot(0,0, type = "n", xlim = c(long_minMax[1] - 10*plotRail, long_minMax[2] + 10*plotRail), ylim = c(lat_minMax[1] - 10*plotRail, lat_minMax[2] + 10*plotRail), xlab = "Longitude", ylab = "Latitude")
# 		for(i in 1:length(trap_data[,1]))
# 		{
# 		polygon(circle_longs[,i], circle_lats[,i])
# 		}
# 		miss <- subset(trap_data, trap_data$count == 0)
# 		hits <- subset(trap_data, trap_data$count > 0)
# 		points(hits$longitude, hits$latitude, pch = 20, col = "green", cex = hits[,3])
# 		points(miss$longitude, miss$latitude, pch = 4, col = "red")
# 		points(crime_loc$longitude, crime_loc$latitude, pch = 18)
# 		points(source_loc$source_long, source_loc$source_lat, col = "blue", pch = 20)
# 	}
# 	return(list(trap_data = trap_data, source_loc = source_loc, crime_loc = crime_loc, n_obs = n_obs, n_true = n_true, sigma_true = sigma_true, hit_data = hit_data, miss_data = miss_data, perSource = perSource, long_minMax = long_minMax, lat_minMax= lat_minMax,
#               n_traps = n_traps, trap_rad_true = trap_rad_true, K_true = K_true, plotRail = plotRail, nega_source1 = nega_source1, nega_source2 = nega_source2, nega_source3 = nega_source3, fakeSlong = fakeSlong, fakeSlat = fakeSlat))
#   C <- C + 1
# }   else{} # sigma_true <- sigma_true[1]}
#
#   }
# }
#
# # #----------------------------------------------------------------
# # MCMC
# PAmcmc <- function(modelSigma, sentinel_data, burnin = 100, iterations = 25e2, Clusters = 1, prop_const = 5, trap_rad = 1, s_prior_shape = 1, s_prior_rate = 1, chains = NULL)
#   {
#   start <- Sys.time()
#
#   # define model parameters
#   trap_density <- sum(sentinel_data$trap_data[,3])
#   sig_init <- 0.5*mean(pairwise_distance(sim$hit_data)$distance_min, na.rm = TRUE)
#
#   # define MCMC parameters
#   propSD <- prop_const*modelSigma
#   PropSD_sig <- prop_const*modelSigma
#
#   ### propose starting values
#   # source locations
#   old_theta <- matrix(NA, nrow = Clusters, ncol = 2 )
#   colnames(old_theta) <- c("Long", "Lat")
#   for (k in 1:Clusters) {
#     old_theta[k,] <- proposal(1, mu1 = 0, mu2 = 0, Sd = propSD)
#   }
#   # sigma
#   old_sigma <- abs(rnorm(1, mean = sig_init, sd = PropSD_sig))
#   old_posterior <- post(theta_x = old_theta[,1], theta_y = old_theta[,2], sigma = old_sigma, data = sentinel_data$trap_data, trap_rad = trap_rad, exp_pop = trap_density)
#
#   # population density
#   new_pop_den <- rgamma(1, shape = s_prior_shape, rate = s_prior_rate)
#
#   ### create objects for storing results
#   joint_lat <- joint_long <- matrix(NA, ncol = Clusters, nrow = iterations)
#   joint_sig <- rep(NA, iterations)
#   pop_densities <- rep(NA, iterations)
#   a_r_rate_theta <- rep(0, iterations)
#   a_r_rate_sig <- rep(0, iterations)
#   propSD_vals <- rep(NA, iterations)
#
#   # run MCMC
#   for (i in 1:iterations) {
#     # report current iteration
#     if ((i %% 500)==0) {
#       cat(paste("  iteration:", i, "\n"))
#     }
#     for (CL in 1:Clusters) {
#     # propose new theta, accept or reject!
#       new_theta <- old_theta
#       new_theta[CL,] <- proposal(1, mu1 = old_theta[CL,1], mu2 = old_theta[CL,2], Sd = propSD) # returns two values a Lon and Lat
#
#       # calculate new joint probability
#       new_posterior <- post(theta_x = new_theta[,1], theta_y = new_theta[,2], sigma = old_sigma, data = sentinel_data$trap_data, trap_rad = trap_rad, exp_pop = new_pop_den)
#
#       # Metropolis-Hastings step
#       ratio <- new_posterior$posterior_val - old_posterior$posterior_val
#       if (log(runif(1)) < ratio) {
#         old_theta[CL,] <- new_theta[CL,]
#         old_posterior <- new_posterior
#         a_r_rate_theta[i] <- 1
#       }
#       if(i < burnin)
#       {
#         theta_shift <- ( (0.234^(a_r_rate_theta[i]))*(-0.766)^(1 - a_r_rate_theta[i])  )/sqrt(i)
#         propSD <- propSD + theta_shift            # Robbins Monro step for source location proposal
#         propSD <- max(1e-4, propSD)
#       }
#     }
#     # propose new sigma, accept or reject!
#     new_sigma <- abs(rnorm(1, mean = old_sigma, sd = PropSD_sig))
#
#     # calculate new joint probability
#     new_posterior <- post(theta_x = old_theta[,1], theta_y = old_theta[,2], sigma = new_sigma, data = sentinel_data$trap_data, trap_rad = trap_rad, exp_pop = new_pop_den)
#
#     # Metropolis-Hastings step
#     ratio <- new_posterior$posterior_val - old_posterior$posterior_val
#     if (log(runif(1)) < ratio) {
#       old_sigma <- new_sigma
#       old_posterior <- new_posterior
#       a_r_rate_sig[i] <- 1
#     }
#     if(i < burnin)
#     {
#       sig_shift <- ( (0.234^(a_r_rate_sig[i]))*(-0.766)^(1 - a_r_rate_sig[i])  )/sqrt(i)
#       PropSD_sig <- PropSD_sig + sig_shift    # Robbins Monro step for sigma proposal
#       PropSD_sig <- max(1e-4, PropSD_sig)
#     }
#     # update s by Gibbs sampling
#     new_z <- post(theta_x = old_theta[,1], theta_y = old_theta[,2], sigma = old_sigma, data = sentinel_data$trap_data, trap_rad = trap_rad, exp_pop = new_pop_den)$z
#     new_pop_den <- rgamma(1, shape = s_prior_shape + trap_density, rate = s_prior_rate + exp(logSum(new_z)))
#
#     # store values of this iteration
#     joint_long[i,] <- old_theta[,1]
#     joint_lat[i,] <- old_theta[,2]
#     joint_sig[i] <- old_sigma
#     pop_densities[i] <- new_pop_den
#     propSD_vals[i] <- propSD
#   }
#
#   long_seq <- seq(sentinel_data$long_minMax[1], sentinel_data$long_minMax[2], (sentinel_data$long_minMax[2] - sentinel_data$long_minMax[1])/500)
#   lat_seq <- seq(sentinel_data$lat_minMax[1], sentinel_data$lat_minMax[2], (sentinel_data$lat_minMax[2] - sentinel_data$lat_minMax[1])/500)
#   #rawSurface <- geoSmooth(joint_long, joint_lat, breaks_lon = long_seq, breaks_lat = lat_seq, lambda = NULL)
#
#   end <- Sys.time()
#   time_elapsed <- end - start
#   print(time_elapsed)
#   # return ... rawSurface
#
#   return(list(joint_long = joint_long, joint_lat = joint_lat, joint_sig = joint_sig, pop_densities = pop_densities, propSD_vals = propSD_vals, long_seq  = long_seq, lat_seq = lat_seq, modelSigma = modelSigma, trap_rad = trap_rad))
#           }
#
# # --------------------------------------------------------------------------------------------------------------------------------
# # set.seed(1)
# # ---------------------------------------
# # Diagnostics code
# # ---------------------------------------
# par(mfrow = c(1,2))
# sim_number <- 20
# TR_seq <- seq(0.2, 1, 0.05)
# diagnosticsData <- matrix(NA, nrow = sim_number*length(TR_seq), ncol = 5)
# colnames(diagnosticsData) <- c("Population", "Estimated_Population", "True_Sigma", "Estimated_Sigma", "Constant")
#
# true_pop.sim  <- 100
# # true_pop.sim <- sample(50:100, 1) # randomly draw true population
# sigma_true.sim <- 10  # randomly draw true sigma
# # sigma_true.sim <- 1             # fix the true sigma value
# K_true.sim <- 2  # randomly draw a number of sources
# traps.sim <- 36
# # traps.sim <- sample(20:100, 1)    # randomly draw a number of traps
# search_cluster.sim <- K_true.sim           # Search for this number of clusters (Presence Absence)
# MCMC_burn <- 1e3
# MCMC_iterations <- 5e3
#
# lat_minMax.sim <- c(0,1) # Set spatial extent to london
# long_minMax.sim <- c(0,1) # set spatial extent to london
#
# # set up sim params
#
# for(i in 1:length(TR_seq))
# {
#   for(j in 1:sim_number)
#   {
#     trap_rad_constant <- TR_seq[i]  # Set the constant, multiplied with sigma to obtain trap radius
#     # sim data given these values
#     sim <- PAsim(n_true = true_pop.sim, trap_const = trap_rad_constant, plotting = F, sigma_true = sigma_true.sim, plotRail = 0.1, K_true = K_true.sim, lat_minMax = lat_minMax.sim, long_minMax = long_minMax.sim,
#                  n_traps = traps.sim, single_count = T, trap_spacing = "uniform", bias = 0)  #, trap_clusters = 4)
#
#     # Define a prior on sigma
#     # priorSigma_mean <- 0.5*mean(pairwise_distance(sim$hit_data)$distance_min, na.rm = TRUE)
#
#     ### Run the presence absence with the fitted sigma and number of clusters from the DPM
#     pa <- PAmcmc(modelSigma = sigma_true.sim, sentinel_data = sim, burnin = MCMC_burn, iterations = MCMC_iterations, Clusters = K_true.sim, prop_const = 3, trap_rad = trap_rad_constant*sigma_true.sim, s_prior_shape = 0.1, s_prior_rate = 0.001, chains = 2)
#
#     est <- hist(pa$pop_densities, 200)
#     est <- mean(est$breaks[which(est$counts == max(est$counts))])
#     est2 <- hist(pa$joint_sig, 200)
#     est2 <- mean(est2$breaks[which(est2$counts == max(est2$counts))])
#     diagnosticsData[sim_number*(i - 1) + j  , 1] <- sim$n_obs
#     diagnosticsData[sim_number*(i - 1) + j  , 2] <- est
#     diagnosticsData[sim_number*(i - 1) + j  , 3] <- sigma_true.sim
#     diagnosticsData[sim_number*(i - 1) + j  , 4] <- est2
#     diagnosticsData[sim_number*(i - 1) + j  , 5] <- TR_seq[i]
#     }
#     print(TR_seq[i])
#     }
#
#
# df <- as.data.frame(diagnosticsData)
# df$Constant <- as.factor(df$Constant)
# save(df, file = "diagnostics_Data.Rdata")
