
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.

#' @useDynLib silverblaze
#' @importFrom Rcpp evalCpp
#' @import roxygen2
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

#----------------------------------------------------------------#
#----------------------------------------------------------------#
#              PRESENCE ABSENCE METROPOLIS HASTINGS              #
#----------------------------------------------------------------#
#----------------------------------------------------------------#
# library(devtools)
# library(plot3D)
# library(RgeoProfile)
# library(MASS)

# install_github("Michael-Stevens-27/silverblaze", ref = "master")
# library(silverblaze)

# rm(list = ls()) #remove all objects

# set.seed(5) # set seed for random number generation

# Calculate destination lat/lon given an origin, a bearing and a great circle distance of travel
# Note that bearing should be in degrees relative to due north, and gc_dist should be in units of kilometres
# From B. VERITY

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
# # Calculate distance between lat long points
# # From B. VERITY
#
# latlon_to_bearing <- function(origin_lat, origin_lon, dest_lat, dest_lon) {
#   if(origin_lat == dest_lat)
#   {
#   dest_lat <- dest_lat + 1e-6
#   }
#   else if(origin_lon == dest_lon)
#   {
#   dest_lon <- dest_lon + 1e-6
#   }
#   # convert input arguments to radians
#   origin_lat <- origin_lat*2*pi/360
#   dest_lat <- dest_lat*2*pi/360
#   origin_lon <- origin_lon*2*pi/360
#   dest_lon <- dest_lon*2*pi/360
#
#   delta_lon <- dest_lon-origin_lon
#
#   # calculate bearing and great circle distance
#   bearing <- atan2(sin(delta_lon)*cos(dest_lat), cos(origin_lat)*sin(dest_lat)-sin(origin_lat)*cos(dest_lat)*cos(delta_lon))
#   gc_angle <- acos(sin(origin_lat)*sin(dest_lat) + cos(origin_lat)*cos(dest_lat)*cos(delta_lon))
#
#   # convert bearing from radians to degrees measured clockwise from due north, and convert gc_angle to great circle distance via radius of earth (km)
#   bearing <- bearing*360/(2*pi)
#   bearing <- (bearing+360)%%360
#   earthRad <- 6371
#   gc_dist <- earthRad*gc_angle
#
#   return(list(bearing=bearing, gc_dist=gc_dist))
# }
#
# #-------------------------------------------------------------------
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
# # joint probability - ADD IN EXPECTED POPULATION DENSITY HERE
# post <- function(theta_x, theta_y, prior_long_mean, prior_lat_mean, prior_Sd, sigma, data, trap_rad, exp_pop) {
#
#   # extract useful parameters
#   K <- length(theta_x)  # K clustrers
#   n <- nrow(data)       # n traps
#
#   # prior log-probability of source i.e. the height on the bivariate normal at the mean of the
#   prior_val <- Prob_data(x = theta_x, y= theta_y, mu_x = prior_long_mean, mu_y = prior_lat_mean, Sd = prior_Sd)
#
#   # log-likelihood
#   heights <- matrix(NA, n, K)
#   for (k in 1:K) {
#     heights[,k] <- dnorm(latlon_to_bearing(data[,2], data[,1], data[,2], theta_x[k])$gc_dist, mean=0, sd=sigma, log=TRUE) + dnorm(latlon_to_bearing(data[,2], data[,1], theta_y[k], data[,1])$gc_dist, mean=0, sd=sigma, log=TRUE)
#   }
#   # Sum poisson params across sources with EQUAL waitings dictated by log(K)
#   z <- log(pi) + 2*log(trap_rad) + apply(heights, 1, logSum) - log(K)
#   log_lambda <- log(exp_pop) + z
#   likeli_val  <- sum( data[,3]*log_lambda - exp(log_lambda) - lfactorial(data[,3]) )
#   posterior_val <- prior_val + likeli_val
#
#   return(list(posterior_val = posterior_val, z = z))
# }
#
# # draw from proposal - CURRENTLY rDPM, NEEDS CHANGING
# proposal <- function(n, mu1, mu2, Sd) {
#   x <- rDPM(1, sigma=Sd, tau=0, priorMean_longitude=mu1, priorMean_latitude=mu2, alpha = 0 )
#   cbind(x$longitude, x$latitude)
# }
#
# # A function that simulates and collects all observed and unobserved data
#
# PAsim <- function(n_true = 100, sigma_true = 0.005, tau_true = NULL, alpha_true = NULL, limits = 1e-1, n_traps = 100, trap_rad_true = 2, K_true = 3, single_count = F, plotting = F)
# {
# 	# define search/area to generate data within
# 	lat_minmax <- long_minmax <- c(-limits, limits)
# 	total_area <- 1
# 	# total_area <- latlon_to_bearing(lat_minmax[1], 0, lat_minmax[2], 0)$gc_dist*latlon_to_bearing(0, long_minmax[1], 0, long_minmax[2])$gc_dist
#
# 	# Draw source locations from a uniform prior (i.e randomly)
# 	source_lat <- runif(K_true, lat_minmax[1], lat_minmax[2])
# 	source_long <- runif(K_true, long_minmax[1], long_minmax[2])
# 	source_loc <- data.frame(source_long = source_long, source_lat = source_lat)
#
# 	# Draw number of observations from a Poisson with rate n_true * total_area and allocate to sources
# 	n_obs <- rpois(1, n_true*total_area)
# 	alloc <- sample(1:K_true, size = n_obs, replace = T)
# 	perSource <- table(alloc)
#
# 	# Draw observation locations from a bivariate norm with standard deviation sigma_true NEEDS EDITING
# 	crime_lat <- unlist(mapply(rnorm, n = perSource, source_lat, sd =  sigma_true))
# 	crime_long <- unlist(mapply(rnorm, n = perSource, source_long, sd = sigma_true))
# 	crime_loc <- data.frame(latitude = crime_lat, longitude = crime_long)
#
# 	# generate traps
#
# 	trap_lat <- runif(n_traps, -limits, limits)
# 	trap_lon <- runif(n_traps, -limits, limits)
# 	trap_loc <- data.frame(trap_lon = trap_lon, trap_lat = trap_lat)
#
# 	# allocate raw data to traps
# 	all_dist <- t(mapply(function(x, y){latlon_to_bearing(trap_lat, trap_lon, x, y)$gc_dist}, x = crime_loc$latitude, y = crime_loc$longitude)) # rows are crimes, columns are traps, entries are distances between crimes and traps
# 	trapCounts <- colSums(all_dist<=trap_rad_true)    # concatonate crimes into traps given the distance is within trap radius
# 	is_observed <- rowSums(all_dist<=trap_rad_true) # Show which crimes were observed
#
# 	# allow for single or multiple trappings of the same observation
# 	if(single_count == T)
# 	{
# 	single_densities <- table(apply(all_dist[which(is_observed > 1),], FUN = which.min, MARGIN = 1))
# 	index <- strtoi(names(single_densities))
# 	trapCounts[index] <- single_densities
# 	}
#
# 	# make final trap data
# 	trap_data <- data.frame(longitude=trap_lon, latitude=trap_lat, count=trapCounts)
#
# 	# plot if plotting returns true
# 	if(plotting == T)
# 	{
# 		c_seq <- seq(0, 4*pi, 0.5)
# 		datax <- trap_rad_true*sin(c_seq)
# 		datay <- trap_rad_true*cos(c_seq)
#
# 		circle_longs <- mapply(function(x, y){cartesian_to_latlon(y, x, datay, datax)$longitude}, x = trap_data$longitude, y = trap_data$latitude)
# 		circle_lats <- mapply(function(x, y){cartesian_to_latlon(y, x, datay, datax)$latitude}, x = trap_data$longitude, y = trap_data$latitude)
#
# 		plot(0,0, type = "n", xlim = c(min(crime_long), max(crime_long)), ylim = c(min(crime_lat), max(crime_lat)), xlab = "Longitude", ylab = "Latitude")
# 		for(i in 1:length(trap_data[,1]))
# 		{
# 		lines(circle_longs[,i], circle_lats[,i])
# 		}
# 		miss <- subset(trap_data, trap_data$count == 0)
# 		hits <- subset(trap_data, trap_data$count > 0)
# 		points(hits$longitude, hits$latitude, pch = 20, col = "green")
# 		points(miss$longitude, miss$latitude, pch = 20, col = "red")
# 		points(crime_loc$longitude, crime_loc$latitude, pch = 4)
# 		points(source_loc$source_long, source_loc$source_lat, col = "blue", pch = 20)
# 	}
# 	return(list(trap_data = trap_data, source_loc = source_loc, crime_loc = crime_loc, n_obs = n_obs, n_true = n_true, sigma_true = sigma_true, tau_true = tau_true,
# 							alpha_true = alpha_true, limits = limits, n_traps = n_traps, trap_rad_true = trap_rad_true, K_true = K_true))
# }
#
# # run an example
# PAsim(plotting = T, trap_rad_true = 1)

# #----------------------------------------------------------------
# MCMC


#----------------------------------------------------------------
