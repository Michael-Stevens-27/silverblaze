
################## TEMP FUNCTIONS TO BE SORTED
#
# dts <- function(x, df, scale=1, log=FALSE) {
#   ret <- lgamma((df+1)/2)-lgamma(df/2)-0.5*log(pi*df*scale^2) - ((df+1)/2)*log(1 + x^2/(df*scale^2))
#   if (!log) { ret <- exp(ret) }
#   return(ret)
# }
#
# latlon_to_cartesian <- function(centre_lat, centre_lon, data_lat, data_lon) {
#
#   # calculate bearing and great circle distance of data relative to centre
#   data_trans <- latlon_to_bearing(centre_lat, centre_lon, data_lat, data_lon)
#
#   # use bearing and distance to calculate cartesian coordinates
#   theta <- data_trans$bearing*2*pi/360
#   d <- data_trans$gc_dist
#   data_x <- d*sin(theta)
#   data_y <- d*cos(theta)
#
#   return(list(x=data_x, y=data_y))
# }
#
# bin2D <- function(x, y, x_breaks, y_breaks) {
#
#   # get number of breaks in each dimension
#   nx <- length(x_breaks)
#   ny <- length(y_breaks)
#
#   # create table of binned values
#   tab1 <- table(findInterval(x, x_breaks), findInterval(y, y_breaks))
#
#   # convert to dataframe and force numeric
#   df1 <- as.data.frame(tab1, stringsAsFactors=FALSE)
#   names(df1) <- c("x", "y", "count")
#   df1$x <- as.numeric(df1$x)
#   df1$y <- as.numeric(df1$y)
#
#   # subset to within breaks range
#   df2 <- subset(df1, x>0 & x<nx & y>0 & y<ny)
#
#   # fill in matrix
#   mat1 <- matrix(0,ny-1,nx-1)
#   mat1[cbind(df2$y, df2$x)] <- df2$count
#
#   # calculate cell midpoints
#   x_mids <- (x_breaks[-1]+x_breaks[-nx])/2
#   y_mids <- (y_breaks[-1]+y_breaks[-ny])/2
#
#   # return output as list
#   output <- list(x_mids=x_mids, y_mids=y_mids, z=mat1)
#   return(output)
# }
#
# geoSmooth <- function (longitude, latitude, breaks_lon, breaks_lat, lambda = NULL)
# {
#     cells_lon <- length(breaks_lon) - 1
#     cells_lat <- length(breaks_lat) - 1
#     centre_lon <- mean(breaks_lon)
#     centre_lat <- mean(breaks_lat)
#     cellSize_lon <- diff(breaks_lon[1:2])
#     cellSize_lat <- diff(breaks_lat[1:2])
#     surface_raw <- bin2D(longitude, latitude, breaks_lon, breaks_lat)$z
#     if (all(surface_raw == 0)) {
#         stop("chosen lat/long window contains no posterior draws")
#     }
#     railSize_lon <- cells_lon
#     railSize_lat <- cells_lat
#     railMat_lon <- matrix(0, cells_lat, railSize_lon)
#     railMat_lat <- matrix(0, railSize_lat, cells_lon + 2 * railSize_lon)
#     surface_normalised <- surface_raw/sum(surface_raw)
#     surface_normalised <- cbind(railMat_lon, surface_normalised,
#         railMat_lon)
#     surface_normalised <- rbind(railMat_lat, surface_normalised,
#         railMat_lat)
#     f1 = fftw2d(surface_normalised)
#     cellSize_trans <- latlon_to_cartesian(centre_lat, centre_lon,
#         centre_lat + cellSize_lat, centre_lon + cellSize_lon)
#     cellSize_trans_lon <- cellSize_trans$x
#     cellSize_trans_lat <- cellSize_trans$y
#     kernel_lon <- cellSize_trans_lon * c(0:floor(ncol(surface_normalised)/2),
#         floor((ncol(surface_normalised) - 1)/2):1)
#     kernel_lat <- cellSize_trans_lat * c(0:floor(nrow(surface_normalised)/2),
#         floor((nrow(surface_normalised) - 1)/2):1)
#     kernel_lon_mat <- outer(rep(1, length(kernel_lat)), kernel_lon)
#     kernel_lat_mat <- outer(kernel_lat, rep(1, length(kernel_lon)))
#     kernel_s_mat <- sqrt(kernel_lon_mat^2 + kernel_lat_mat^2)
#     if (is.null(lambda)) {
#         lambda_step <- min(cellSize_trans_lon, cellSize_trans_lat)/5
#         lambda_vec <- lambda_step * (1:100)
#     }
#     else {
#         lambda_vec <- lambda
#     }
#     cat("Smoothing posterior surface")
#     flush.console()
#     logLike <- -Inf
#     for (i in 1:length(lambda_vec)) {
#         if (i > 1) {
#             cat(".")
#             flush.console()
#         }
#         lambda_this <- lambda_vec[i]
#         kernel <- dts(kernel_s_mat, df = 3, scale = lambda_this)
#         f2 = fftw2d(kernel)
#         f3 = f1 * f2
#         f4 = Re(fftw2d(f3, inverse = T))/length(surface_normalised)
#         f5 <- f4 - surface_normalised * dts(0, df = 3, scale = lambda_this)
#         f5[f5 < 0] <- 0
#         f5 <- f5/sum(f4)
#         f6 <- surface_normalised * log(f5)
#         if (sum(f6, na.rm = T) < logLike) {
#             (break)()
#         }
#         logLike <- sum(f6, na.rm = T)
#     }
#     if (is.null(lambda)) {
#         cat(paste("\nmaximum likelihood lambda = ", round(lambda_this,
#             3), sep = ""))
#     }
#     f4 <- f4[, (railSize_lon + 1):(ncol(f4) - railSize_lon)]
#     f4 <- f4[(railSize_lat + 1):(nrow(f4) - railSize_lat), ]
#     return(f4)
# }
###############################################
#------------------------------------------------
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

#----------------------------------------------------------------#
#----------------------------------------------------------------#
#              PRESENCE ABSENCE METROPOLIS HASTINGS              #
#----------------------------------------------------------------#
#----------------------------------------------------------------#
#
# library(devtools)
# library(plot3D)
# library(RgeoProfile)
# library(MASS)
# library(ggplot2)
# library(fftwtools)
#
# # install_github("Michael-Stevens-27/silverblaze", ref = "master")
# # library(silverblaze)
#
# # rm(list = ls()) #remove all objects
#
# # set.seed(5) # throws out ratio error
#
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
# # Calculate distance between lat long points
# # From B. VERITY
#
# latlon_to_bearing <- function(origin_lat, origin_lon, dest_lat, dest_lon) {
#   # if(origin_lat == dest_lat)
#   # {
#   # dest_lat <- dest_lat + 1e-6
#   # }
#   # else if(origin_lon == dest_lon)
#   # {
#   # dest_lon <- dest_lon + 1e-6
#   # }
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
# # joint probability -
# post <- function(theta_x, theta_y, prior_long_mean, prior_lat_mean, prior_Sd = NULL, sigma, data, trap_rad, exp_pop) {
#
#   # extract useful parameters
#   K <- length(theta_x)  # K clustrers
#   n <- nrow(data)       # n traps
#
#   # prior log-probability of source - currently undefined
#   # prior_val <- Prob_data(x = theta_x, y= theta_y, mu_x = prior_long_mean, mu_y = prior_lat_mean, Sd = prior_Sd)
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
#   # posterior_val <- prior_val + likeli_val # prior currently undefined
#   posterior_val <- likeli_val
#   return(list(posterior_val = posterior_val, z = z))
# }
#
# # draw from proposal - CURRENTLY rDPM, NEEDS CHANGING
# proposal <- function(n, mu1, mu2, Sd) {
#   points <- rnorm_sphere(n = n, centre_lon = mu1, centre_lat = mu2, sigma = Sd)
#   cbind(longitude = points$longitude, latitude = points$latitude)
# }
#
# # A function that simulates and collects all observed and unobserved data
# PAsim <- function(n_true = 100, sigma_true = 1, tau_true = NULL, alpha_true = NULL, long_minMax = c(0, 0.5), lat_minMax = c(0, 0.5), n_traps = 100, trap_rad_true = 0.5,
#                   K_true = 3, single_count = F, plotting = F, plotRail = 0.05)
# {
#   # define search/area to generate data within
# 	total_area <- 1
# 	# total_area <- latlon_to_bearing(lat_minmax[1], 0, lat_minmax[2], 0)$gc_dist*latlon_to_bearing(0, long_minmax[1], 0, long_minmax[2])$gc_dist
#
# 	# Draw source locations from a uniform prior (i.e randomly)
# 	source_lat <- runif(K_true, lat_minMax[1], lat_minMax[2])
# 	source_long <- runif(K_true, long_minMax[1], long_minMax[2])
# 	source_loc <- data.frame(source_long = source_long, source_lat = source_lat)
#
# 	# Draw number of observations from a Poisson with rate n_true * total_area and allocate to sources
# 	n_obs <- rpois(1, n_true*total_area)
# 	alloc <- sample(1:K_true, size = n_obs, replace = T)
# 	perSource <- table(alloc)
#
# 	# Draw observation locations from a bivariate norm with standard deviation sigma_true NEEDS EDITING
#   crime_long <- unlist(mapply(rnorm_sphere, n = perSource, centre_lat = source_lat, centre_lon = source_long, sigma = sigma_true)[1,])
#   crime_lat <- unlist(mapply(rnorm_sphere, n = perSource, centre_lat = source_lat, centre_lon = source_long, sigma = sigma_true)[2,])
# 	crime_loc <- data.frame(longitude = crime_long, latitude = crime_lat)
#
# 	# generate traps
#   trap_lat <- runif(n_traps, lat_minMax[1], lat_minMax[2])
# 	trap_lon <- runif(n_traps, long_minMax[1], long_minMax[2])
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
#   hit_data <- subset(trap_data, trap_data$count >0 )
#   miss_data <- subset(trap_data, trap_data$count ==0 )
#
# 	# plot if plotting returns true
# 	if(plotting == T)
# 	{
# 		c_seq <- seq(0, 2*pi, 0.5)
# 		datax <- trap_rad_true*sin(c_seq)
# 		datay <- trap_rad_true*cos(c_seq)
#
# 		circle_longs <- mapply(function(x, y){cartesian_to_latlon(y, x, datay, datax)$longitude}, x = trap_data$longitude, y = trap_data$latitude)
# 		circle_lats <- mapply(function(x, y){cartesian_to_latlon(y, x, datay, datax)$latitude}, x = trap_data$longitude, y = trap_data$latitude)
#
# 		plot(0,0, type = "n", xlim = c(long_minMax[1] - plotRail, long_minMax[2] + plotRail), ylim = c(lat_minMax[1] - plotRail, lat_minMax[2] + plotRail), xlab = "Longitude", ylab = "Latitude")
# 		for(i in 1:length(trap_data[,1]))
# 		{
# 		polygon(circle_longs[,i], circle_lats[,i])
# 		}
# 		miss <- subset(trap_data, trap_data$count == 0)
# 		hits <- subset(trap_data, trap_data$count > 0)
# 		points(hits$longitude, hits$latitude, pch = 20, col = "green")
# 		points(miss$longitude, miss$latitude, pch = 20, col = "red")
# 		points(crime_loc$longitude, crime_loc$latitude, pch = 4)
# 		points(source_loc$source_long, source_loc$source_lat, col = "blue", pch = 20)
# 	}
# 	return(list(trap_data = trap_data, source_loc = source_loc, crime_loc = crime_loc, n_obs = n_obs, n_true = n_true, sigma_true = sigma_true, tau_true = tau_true, hit_data = hit_data, miss_data = miss_data,
# 							alpha_true = alpha_true, long_minMax = long_minMax, lat_minMax= lat_minMax, n_traps = n_traps, trap_rad_true = trap_rad_true, K_true = K_true, plotRail = plotRail))
# }
#
# set.seed(6)
# prac_sim <- PAsim(trap_rad_true = 0.23, plotting= T, sigma_true = 1, plotRail = 0.05, K_true = 5, lat_minMax = c(51.4574, 51.5574), long_minMax = c(-0.1777,-0.0777))
# # #----------------------------------------------------------------
# # MCMC
#
# # define model parameters
# sigma <- 0.5*mean(pairwise_distance(prac_sim$trap_data[,1:2][which(prac_sim$trap_data[,3]>0),])$distance_min, na.rm = T)
# tau <- 10*sigma
# prior_long_mean <- mean(prac_sim$trap_data[,1][which(prac_sim$trap_data[,3]>0)])
# prior_lat_mean <- mean(prac_sim$trap_data[,2][which(prac_sim$trap_data[,3]>0)])
# K <- 5
# small_const <- 0.2
# # trap_rad <- trap_rad_true
# trap_rad <- small_const*sigma
# trap_density <- sum(prac_sim$trap_data[,3])
#
# # define MCMC parameters
# burnin <- 500
# iterations <- 1e4
# iteration_vec <- 1:iterations
# propSD <- tau
# start_propSD <- tau
#
# # propose starting values
# old_theta <- matrix(NA, nrow=K, ncol=2)
# for (k in 1:K) {
#   old_theta[k,] <- proposal(1, mu1 = prior_long_mean, mu2 = prior_lat_mean, Sd = propSD)
# }
# old_posterior <- post(theta_x = old_theta[,1], theta_y = old_theta[,2], prior_long_mean = prior_long_mean, prior_lat_mean = prior_lat_mean,
#                       sigma = sigma, data = prac_sim$trap_data, trap_rad = trap_rad, exp_pop = trap_density)
#
# new_pop_den <- rgamma(1, trap_density, 1)
#
# # create objects for storing results
# all_theta_lat <- all_theta_long <- matrix(NA, ncol = K, nrow = iterations)
# pop_densities <- rep(NA, iterations)
# a_r_rate <- rep(0, iterations)
# sigma_vals <- rep(NA, iterations)
# #
# # run MCMC
# for (i in 1:iterations) {
#   #  i <- 1
#   # report current iteration
#   if ((i %% 5e2)==0) {
#     cat(paste("  iteration:", i, "\n"))
#   }
#   # propose new theta and expected population
#   for (k in 1:K) {
#     new_theta <- old_theta
#     new_theta[k,] <- proposal(1, mu1 = old_theta[k,1], mu2 = old_theta[k,2], Sd = propSD)
#
#     new_z <- post(theta_x = new_theta[,1], theta_y = new_theta[,2], prior_long_mean = prior_long_mean, prior_lat_mean = prior_lat_mean,
#                   sigma = sigma, data = prac_sim$trap_data, trap_rad = trap_rad, exp_pop = new_pop_den)$z
#
#     new_pop_den <- rgamma(1, shape = 2*trap_density, rate = 1 + (prac_sim$n_traps)*exp(logSum(new_z)))
#
#     # calculate new joint probability
#     new_posterior <- post(theta_x = new_theta[,1], theta_y = new_theta[,2], prior_long_mean = prior_long_mean, prior_lat_mean = prior_lat_mean,
#                           sigma = sigma, data = prac_sim$trap_data, trap_rad = trap_rad, exp_pop = new_pop_den)
#
#     # Metropolis-Hastings step
#     ratio <- new_posterior$posterior_val - old_posterior$posterior_val
#     if (log(runif(1)) < ratio) {
#       old_theta[k,] <- new_theta[k,]
#       old_posterior <- new_posterior
#       a_r_rate[i] <- 1
#     }
#     if(i < burnin)
#     {
#       if(a_r_rate[i] == 1)
#       {propSD <- propSD + (1 - 0.23)/sqrt(iteration_vec[i])}
#       else{propSD <- propSD - 0.23/sqrt(iteration_vec[i])
#         propSD <- max(1e-4, propSD)}
#         print(propSD)
#     }
#     propSD <- propSD
#   }
#   # store values of this iteration
#   all_theta_long[i,] <- old_theta[,1]
#   all_theta_lat[i,] <- old_theta[,2]
#   pop_densities[i] <- new_pop_den
#   sigma_vals[i] <- propSD
# }
#
#   lines(all_theta_long, all_theta_lat, col = "cyan")
#
#   prac_sim$long_minMax
#
#   long_seq <- seq(prac_sim$long_minMax[1], prac_sim$long_minMax[2], 0.0005)
#   lat_seq <- seq(prac_sim$lat_minMax[1], prac_sim$lat_minMax[2], 0.0005)
#   temp_smooth <- geoSmooth(all_theta_long, all_theta_lat, breaks_lon = long_seq, breaks_lat = lat_seq, lambda = NULL)
#
# # contour(x = seq(prac_sim$long_minMax[1], prac_sim$long_minMax[2] - 0.001, 0.001), y= seq(prac_sim$lat_minMax[1], prac_sim$lat_minMax[2]- 0.001, 0.001), t(temp_smooth), add = T)
#  trap_data <- geoData(prac_sim$hit_data$longitude, prac_sim$hit_data$latitude)
#  source_dataM <- geoDataSource(prac_sim$source_loc$source_long, prac_sim$source_loc$source_lat)
#  temp_param <- geoParams(data = trap_data, sigma_mean = 1, sigma_var = 0, longitude_cells = length(long_seq) - 1, latitude_cells = length(lat_seq) - 1 )
#  temp_param$output$longitude_minMax <- c(prac_sim$long_minMax[1], prac_sim$long_minMax[2])
#  temp_param$output$latitude_minMax <- c(prac_sim$lat_minMax[1], prac_sim$lat_minMax[2])
#  temp_profile <- geoProfile(temp_smooth)
#
#  x11()
#  geoPlotMap(params = temp_param, data = trap_data, surface = temp_profile, crimeCol = "green", source = source_dataM, sourceCol = "red", breakPercent = seq(0,50,5))
