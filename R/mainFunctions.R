#----------------------------------------------------------------#
#----------------------------------------------------------------#
#              PRESENCE ABSENCE METROPOLIS HASTINGS              #
#----------------------------------------------------------------#
#----------------------------------------------------------------#

# The following commands ensure that package dependencies are listed in the NAMESPACE file.

#' @useDynLib silverblaze
#' @importFrom Rcpp evalCpp
NULL

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
#' latlon_to_cartesian(centre_lat = centre_latLon[1], centre_lon = centre_latLon[2],
#' data_lat = data_latLon[1], data_lon = data_latLon[2])
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
#' latlon_to_bearing(origin_lat = origin_latLon[1], origin_lon = origin_latLon[2],
#' dest_lat = dest_latLon[1], dest_lon = dest_latLon[2])
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

#------------------------------------------------
#' bearing_to_latlon
#'
#' Calculate destination lat/lon given an origin, a bearing and a greater circle distance of travel
#'
#' @param origin_lat The origin latitude
#' @param origin_lon The origin longitude
#' @param bearing   The angle in degrees relative to due north
#' @param gc_dist   The greater circle distance in kilometres
#'
#' @export
#' @examples
#' origin_latLon <- c(0,0)
#' some_bearing <- 50
#' some_gc_dist <- 27
#' bearing_to_latlon(origin_lat = origin_latLon[1], origin_lon = origin_latLon[2],
#' bearing = some_bearing, gc_dist = some_gc_dist)
#  From B. VERITY

bearing_to_latlon <- function(origin_lat, origin_lon, bearing, gc_dist) {

  # convert origin_lat, origin_lon and bearing from degrees to radians
  origin_lat <- origin_lat*2*pi/360
  origin_lon <- origin_lon*2*pi/360
  bearing <- bearing*2*pi/360

  # calculate new lat/lon using great circle distance
  earthRad <- 6371
  new_lat <- asin(sin(origin_lat)*cos(gc_dist/earthRad) + cos(origin_lat)*sin(gc_dist/earthRad)*cos(bearing))
  new_lon <- origin_lon + atan2(sin(bearing)*sin(gc_dist/earthRad)*cos(origin_lat), cos(gc_dist/earthRad)-sin(origin_lat)*sin(new_lat))

  # convert new_lat and new_lon from radians to degrees
  new_lat <- new_lat*360/(2*pi)
  new_lon <- new_lon*360/(2*pi)

  return(list(longitude=new_lon, latitude=new_lat))
}

#------------------------------------------------
#' cartesian_to_latlon
#'
#' Transform Cartesian Co-ordinates into Latitude and Longitude points assuming a Lat_Long centre
#'
#' @param centre_lat The centre latitude
#' @param centre_lon The centre longitude
#' @param data_x Cartesian x co-ordinate
#' @param data_y Cartesian y co-ordinate
#'
#' @export
#' @examples
#' centre <- c(1,5)
#' data <- c(2,3)
#' cartesian_to_latlon(centre_lat = centre[1], centre_lon = centre[2], data_x = data[1], data_y = data[2])
# From B. VERITY

cartesian_to_latlon <- function(centre_lat, centre_lon, data_x, data_y) {

  # calculate angle and euclidian distance of all points relative to origin
  d <- sqrt(data_x^2+data_y^2)
  theta <- atan2(data_y,data_x)

  # convert theta to bearing relative to due north
  theta <- theta*360/(2*pi)
  theta <- (90-theta)%%360

  # use bearing and great circle distance to calculate lat/lon relative to an origin point
  data_trans <- bearing_to_latlon(centre_lat, centre_lon, theta, d)

  return(list(longitude=data_trans$longitude, latitude=data_trans$latitude))
}

#------------------------------------------------
#' pairwise_distance
#'
#' Given two sets of points calculate the distance between each pair of points. Returns a matrix whose entries are distances between pairs of points. Also returns the nearest neighbour point in the form of the minimum rowise distance of this matrix.
#'
#' @param points two sets of x and y points all in latitude and longitude
#'
#' @export
#' @examples
#' some_points <- cbind(1:10, 1:10, 11:20, 11:20)
#' pairwise_distance(some_points)
# From B. VERITY

pairwise_distance <- function(points){
	distance <- matrix(NA, ncol = length(points[,1]), nrow = length(points[,1]))

	for(i in 1:length(points[,1]))
	{
	for(j in 1:i)
		{
		x1 <- points[i, 1]
		y1 <- points[i, 2]
		x2 <- points[j, 1]
		y2 <- points[j, 2]
		dist <- latlon_to_bearing(y1, x1, y2, x2)$gc_dist
		distance[i,j] <- dist
		}
	}
	diag(distance)  <- NA
	distance[distance==0] <- NA
	distance_min <- apply(distance, 1, min, na.rm = TRUE)
	distance_min[distance_min =="Inf"] <- NA
	SD_TR <- list(distance = distance, distance_min = distance_min)
	return(SD_TR)
}

#------------------------------------------------
#' rnorm_sphere
#'
#' Draw from normal distribution converted to spherical coordinate system. Points are first drawn from an ordinary cartesian 2D normal distribution.
#' The distances to points are then assumed to be great circle distances, and are combined with a random bearing from the point {centre_lat, centre_lon}
#' to produce a final set of lat/lon points. Note that this is not a truly spherical normal distribution, as the domain of the distribution is not the sphere -
#' rather it is a transformation from one coordinate system to another that is satisfactory when the curvature of the sphere is not severe.
#'
#' @param n The number of points to draw
#' @param centre_lat The mean latitude of the normal being drawn from
#' @param centre_lon The mean longitude of the normal being drawn from
#' @param sigma The standard deviation of the normal being drawn from
#'
#' @export
#' @examples
#' rnorm_sphere(n = 100, centre_lat = 0, centre_lon = 0, sigma = 1)
# From B. VERITY

rnorm_sphere <- function(n, centre_lat, centre_lon, sigma) {
	x <- rnorm(n,sd=sigma)
	y <- rnorm(n,sd=sigma)
	output <- cartesian_to_latlon(centre_lat, centre_lon, x, y)
	return(output)
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
#' x <- log(seq(5e-49, 5e-48, 1e-49))
#' exp(sum(x))
#' exp(logSum(x))

logSum <- function(x){
	xmax <- max(x)
	xmax + log(sum(exp(x-xmax)))
}

# Presence Absence Functions
#-------------------------------------------------------------------
#' bvnorm
#'
#' Extract the log density on the bivariate normal at (x, y) centred at (mu_x, mu_y)
#'
#' @param x 		The location of the point whose density we require (longitude)
#' @param y 		The location of the point whose density we require (latitude)
#' @param mu_x  The mean of the bivariate normal (longitude)
#' @param mu_y  The mean of the bivariate normal (latitude)
#' @param sd 		The standard deviation of the bivariate normal in kilometres
#'
#' @export
#' @examples
#' bvnorm(x = 1, y = 1, mu_x = 0, mu_y = 0, Sd = 10)

bvnorm <- function(x, y, mu_x, mu_y, Sd, log = T) {
  param <- dnorm(latlon_to_bearing(y, x, y, mu_x)$gc_dist, mean = 0, sd = Sd, log = log) + dnorm(latlon_to_bearing(y, mu_x, mu_y, mu_x)$gc_dist, mean = 0, sd = Sd, log = log)
  return(param)
}

#-------------------------------------------------------------------
#' PAjoint
#'
#' Calculate the logged joint probability of observing the data. This acts as the posterior draw. Also extract the logged density of the bivariate normal at each trap location.
#'
#' @param theta_x 					The mean longitude for the source we condition on
#' @param theta_y 					The mean latitude for the source we condition on
#' @param sigma  						The standard deviation of the bivariate normal (source) we condition on
#' @param data  						The trap data consisting of locations in lat/long and associated densities
#' @param trap_rad 					The trap radius we use to approximate the area under the bivariate normal (in kilometres)
#' @param exp_pop						The expected population over our entire area
#' @param priorLat					The mean Latitude of the prior on source locations - normally the mean of hit locations
#' @param priorLon					The mean Longitude of the prior on source locations - normally the mean of hit locations
#' @param prior 						Specify which parameter a prior is being set for
#' @param meanSigprior			The mean of the prior on sigma
#' @param propLoc 					The proposed location for a source
#' @param priorSD 					The standard devitation for any prior all are set to normals
#'
#' @export
#' @examples
#' trap_sim <- cbind(runif(10), runif(10), sample(0:5, 10, replace = T))
#' PAjoint(theta_x = 0, theta_y = 0, sigma = 10, data = trap_sim, trap_rad = 2, exp_pop = 25,
#' priorLat = 0.5, priorLon = 0.5, prior = "sourcePrior"
#' , meanSigprior = NULL, propLoc = c(0.1, 0.3), priorSD = 10)

PAjoint <- function(theta_x, theta_y, sigma, data, trap_rad, exp_pop, priorLat, priorLon, prior = NULL, meanSigprior, propLoc = NULL, priorSD = NULL) {
  # extract useful parameters
  K <- length(theta_x)  # K clusters
  n <- nrow(data)       # n traps

	# Sum poisson params across sources with EQUAL waitings dictated by log(K)
	heights <- mapply(function(x, y){bvnorm(data[,1], data[,2], x, y, Sd = sigma, log = T)}, x = theta_x, y = theta_y)
	z <- log(pi) + 2*log(trap_rad) + apply(heights, 1, logSum) - log(K)

	# log-likelihood
  log_lambda <- log(exp_pop) + z
  likeli_val  <- sum( data[,3]*log_lambda - exp(log_lambda) - lfactorial(data[,3]) )
	if(prior == "sourcePrior")
		{
		prior_val <- bvnorm(x = propLoc[1], y = propLoc[2], mu_x = priorLon, mu_y = priorLat, Sd = priorSD, log = T)
		} else if(prior == "sigPrior") {
							prior_val <- dnorm(sigma, mean = meanSigprior, sd = priorSD, log = TRUE)
							} else if(prior == "None"){
								prior_val <- 0
	}
  posterior_val <- prior_val + likeli_val
  return(list(posterior_val = posterior_val, z = z))
}

#------------------------------------------------
#' proposal
#'
#' Propose new parameter values (specifically a source location) from a bivaraite normal.
#'
#' @param n The number of new points to propose
#' @param mu1 The centre of the bivariate norm to draw from (longitude)
#' @param mu2 The centre of the bivariate norm to draw from (latitude)
#' @param Sd  The standard deviation of the bivariate norm to draw from (kilometres)
#'
#' @export
#' @examples
#' proposal(n = 16, mu1 = 0, mu2 = 0, Sd = 5)

proposal <- function(n, mu1, mu2, Sd) {
  points <- rnorm_sphere(n = n, centre_lon = mu1, centre_lat = mu2, sigma = Sd)
  cbind(longitude = points$longitude, latitude = points$latitude)
}

#------------------------------------------------
#' PAsim
#'
#' A function that simulates and collects all observed and unobserved data
#'
#' @param n_true 					The underlying population size
#' @param sigma_true			The underlying sigma value of the data
#' @param K_true					The true number of sources
#' @param n_traps					The number of traps to be used
#' @param trap_spacing		Set to either "uniform", "random", or "cluster" describes the trap configuration
#' @param trapExtent			When trap_spacing is clustered, this describes the distance from each trap to the centre of that cluster
#' @param trap_clusters		When trap_spacing is clustered, this is sets the number of traps within a cluster
#' @param long_minMax			The longitudinal extent of the data (to be set using this function)
#' @param lat_minMax			The latitudinal extent of the data (to be set using this function)
#' @param plotRail				A buffer to be placed around the data when plotting
#' @param trap_const 			This governs the trap radiu, which is the trap_const multiplied by the true sigma
#' @param single_count		Set to TRUE or FALSE this governs if a n individual can be caught by multiple traps
#' @param plotting				Set to TRUE or FALSE should the data want to be plotted
#' @param bias						Set to 0 or 1, if set to 1 two additional sources will be created, one surrounded by empty traps, and one far away from the traps
#'
#' @export
#' @examples
#' PAsim(n_true = 100, sigma_true = 1, K_true = 3, n_traps = 10, trap_spacing = "cluster",
#' trapExtent = 0.5, long_minMax = c(0, 0.1), lat_minMax = c(0, 0.1), plotRail = 0.05,
#' trap_const = 0.5, trap_clusters = 4, single_count = T, plotting = T, bias = 1)

PAsim <- function(n_true = 1000, sigma_true = 2, K_true = 3, n_traps = 64, trap_spacing = "cluster", trapExtent = 3, long_minMax = c(0, 0.5), lat_minMax = c(0, 0.5), plotRail = 0.05, trap_const = 0.5, trap_clusters = 4, single_count = T, plotting = T, bias = 0)
  {
  # define search/area to generate data within
	total_area <- 1
  trap_rad_true <- trap_const*sigma_true

  # cataegorise sources in real, fake or unseen
  sTyp <- rep("Real", K_true + 2*bias)

	# Draw source locations from a uniform prior (i.e randomly) with or without bias
  if(bias == 1)
  {
  source_lat <- runif(K_true, lat_minMax[1], lat_minMax[1] + 0.5*(lat_minMax[2] - lat_minMax[1]) )
  source_long <- runif(K_true, long_minMax[1], long_minMax[1] + 0.5*(long_minMax[2] - long_minMax[1] ) )

  # create a source within traps that has no individuals associated with it
  fakeSlong <- runif(1, long_minMax[1] + 0.99*(long_minMax[2] - long_minMax[1]), long_minMax[2])
  fakeSlat <- runif(1, lat_minMax[1], lat_minMax[1] + 0.25*(lat_minMax[2] - lat_minMax[1]))

  # create a source outside of traps that has individuals associated with it
  source_lat <- c(source_lat, runif(1, lat_minMax[1] + 0.9*(lat_minMax[2] - lat_minMax[1]), lat_minMax[2]), fakeSlat)
  source_long <- c(source_long, runif(1, long_minMax[1], long_minMax[2] - 0.5*(long_minMax[2] - long_minMax[1])), fakeSlong)

  # Draw number of observations from a Poisson with rate n_true * total_area and allocate to sources
  n_obs <- rpois(1, n_true*total_area)
  alloc <- sample(1:(K_true + 1), size = n_obs, replace = T)

	# Update the sources type object to include those extra two
	sTyp[K_true + 2] <- "Fake"
  sTyp[K_true + 1] <- "Unseen"
	 } else{
  # if there is no bias create the objects as normal
	  source_lat <- runif(K_true, lat_minMax[1], lat_minMax[2])
    source_long <- runif(K_true, long_minMax[1], long_minMax[2])
    n_obs <- rpois(1, n_true*total_area)
    alloc <- sample(1:(K_true), size = n_obs, replace = T)
  }
	source_loc <- data.frame(source_long = source_long, source_lat = source_lat, sTyp = sTyp)
	# create an object to state individuals per sources
	perSource <-  mapply(function(x){length(which(x == alloc))}, x = 1:(K_true + bias))
	if(bias ==  1)
	{
	# include the number of individuals associated with the fake source, zero.
	perSource <-  c(perSource, 0)
	}

  # Draw individuals' locations from a bivariate norm with standard deviation sigma_true NEEDS EDITING
  indiv_long <- unlist(mapply(rnorm_sphere, n = perSource, centre_lat = source_lat, centre_lon = source_long, sigma = sigma_true)[1,])
  indiv_lat <- unlist(mapply(rnorm_sphere, n = perSource, centre_lat = source_lat, centre_lon = source_long, sigma = sigma_true)[2,])
	indiv_loc <- data.frame(longitude = indiv_long, latitude = indiv_lat)

  # set trap spacing to random, uniform, or clustered/uniform (trap_clusters dictates the number of traps in each uniform location)
  if(trap_spacing == "random")
  {
    trap_lat <- runif(n_traps, lat_minMax[1], lat_minMax[2] - (lat_minMax[2]/2 - lat_minMax[1]/2)*bias)
    trap_lon <- runif(n_traps, long_minMax[1], long_minMax[2])
    trap_loc <- data.frame(trap_lon = trap_lon, trap_lat = trap_lat)

  } else {
    trap_lat <- seq(lat_minMax[1], (lat_minMax[1]*(bias) + lat_minMax[2])/(2^bias), (lat_minMax[2]/(2^bias) - lat_minMax[1])/ceiling(sqrt(n_traps)))
		trap_lon <- seq(long_minMax[1], long_minMax[2], (long_minMax[2] - long_minMax[1])/ceiling(sqrt(n_traps)))
    trap_loc <- expand.grid(trap_lon, trap_lat)
	}
  if(trap_spacing == "cluster") {
    #same as method for unifrom, but then extend to produce clustered trap structure
    c_trap <- seq(2*pi/trap_clusters, 2*pi, 2*pi/trap_clusters )
		polx <- trapExtent*sin(c_trap)
		poly <- trapExtent*cos(c_trap)
		trap_lon <- mapply(function(x, y){cartesian_to_latlon(y, x, poly, polx)$longitude}, x = trap_loc[,1], y = trap_loc[,2])
		trap_lat <- mapply(function(x, y){cartesian_to_latlon(y, x, poly, polx)$latitude}, x = trap_loc[,1], y = trap_loc[,2])
    trap_loc <- cbind(c(trap_lon), c(trap_lat))
	}
	if(length(trap_lat) == 0 )
	{
	stop("trap_spacing must be set to 'uniform','random' or 'cluster'")
	}
	# allocate raw data to traps
  all_dist <- mapply(function(x, y){latlon_to_bearing(trap_loc[,2], trap_loc[,1], x, y)$gc_dist}, x = indiv_loc$latitude, y = indiv_loc$longitude)
  all_dist <- t(all_dist) # rows are individuals, columns are traps, entries are distances between individuals and traps
  trapCounts <- colSums(all_dist<=trap_rad_true)  # concatonate individuals into traps given the distance is within trap radius
  is_observed <- rowSums(all_dist<=trap_rad_true)  # count the number of observed observations

	if(length(which(trapCounts > 0)) < 5) # check there are more than 5 observations overall
  {
	stop("Not enough observations, increase sigma_true, n_traps or reduce long_minMax/lat_minMax ")
	}
  # allow for single or multiple trappings of the same observation
	if(single_count == T)
	{
  if(max(is_observed) == 1) # check traps do not already only contain a single observation
  {} else if(length(which(is_observed > 1)) > 1) { # check that there is more than one trap with multiple observations
  TF_count <- all_dist<=trap_rad_true
  for(G in 1:length(which(is_observed > 1)))
  {
  # for each individual observed more than once, pull out the index of traps that observe it
  anti_index <- which(TF_count[which(is_observed >1)[G],] == T)
  # remove the closest trap from this index
  to_remove <- which.min(all_dist[which(is_observed > 1),][G,])
  anti_index <- anti_index[! anti_index %in% to_remove]
  # reduce the trap counts of those traps that should not be observing a particular individual
  trapCounts[anti_index] <- trapCounts[anti_index] - 1
  }

  } else { # if there is a single trap with multiple observations, then replace with a single obs
  single_index <- which.min(all_dist[which(is_observed > 1),])
  trapCounts[single_index] <- 1
  }
	}

	# make final trap data
	trap_data <- data.frame(longitude=trap_loc[,1], latitude=trap_loc[,2], count=trapCounts)
  hit_data <- subset(trap_data, trap_data$count > 0)
  miss_data <- subset(trap_data, trap_data$count ==0)

  # record number of hits and miss' within 1, 2 and 3 sd's of sources
  source_hits_dist <- t(mapply(function(x, y){latlon_to_bearing(hit_data[,2], hit_data[,1], x, y)$gc_dist}, x = source_loc$source_lat, y = source_loc$source_long)) # rows are sources, columns are hits, entries are distances between sources and hits
  hitsWSD <- cbind(rowSums(source_hits_dist <= sigma_true), rowSums(source_hits_dist <= 2*sigma_true), rowSums(source_hits_dist <= 3*sigma_true))

  source_miss_dist <- t(mapply(function(x, y){latlon_to_bearing(miss_data[,2], miss_data[,1], x, y)$gc_dist}, x = source_loc$source_lat, y = source_loc$source_long)) # rows are sources, columns are miss, entries are distances between sources and miss
  missWSD <- cbind(rowSums(source_miss_dist <= sigma_true) , rowSums(source_miss_dist <= 2*sigma_true), rowSums(source_miss_dist <= 3*sigma_true))

	# plot if plotting returns true
	if(plotting == T)
	{
		c_seq <- seq(0, 2*pi, 0.25)
		datax <- trap_rad_true*sin(c_seq)
		datay <- trap_rad_true*cos(c_seq)
		circle_longs <- mapply(function(x, y){cartesian_to_latlon(y, x, datay, datax)$longitude}, x = trap_data$longitude, y = trap_data$latitude)
		circle_lats <- mapply(function(x, y){cartesian_to_latlon(y, x, datay, datax)$latitude}, x = trap_data$longitude, y = trap_data$latitude)

		plot(0,0, type = "n", xlim = c(long_minMax[1] - plotRail, long_minMax[2] + plotRail), ylim = c(lat_minMax[1] - plotRail, lat_minMax[2] + plotRail), xlab = "Longitude", ylab = "Latitude")
		for(i in 1:length(trap_data[,1]))
		{
		polygon(circle_longs[,i], circle_lats[,i])
		}
		miss <- subset(trap_data, trap_data$count == 0)
		hits <- subset(trap_data, trap_data$count > 0)
		points(hits$longitude, hits$latitude, pch = 20, col = "green", cex = 1) # hits[,3])
		points(miss$longitude, miss$latitude, pch = 4, col = "red")
		points(indiv_loc$longitude, indiv_loc$latitude, pch = 18, cex = 0.75)
		points(source_loc$source_long, source_loc$source_lat, col = "blue", pch = 15, cex = 1.5)
	}
	return(list(trap_data = trap_data, source_loc = source_loc, indiv_loc = indiv_loc, n_obs = n_obs, n_true = n_true, sigma_true = sigma_true, hit_data = hit_data, miss_data = miss_data, perSource = perSource, long_minMax = long_minMax, lat_minMax= lat_minMax, trap_spacing = trap_spacing,
              n_traps = n_traps, trap_rad_true = trap_rad_true, K_true = K_true, plotRail = plotRail, hitsWSD = hitsWSD,  missWSD = missWSD))
  }

#------------------------------------------------
#' PAsim
#'
#' A function that simulates and collects all observed and unobserved data
#'
#' @param trapData				 The main input into the model consisting of trap longitude, latitude and densities
#' @param burnin					 The MCMC's burn in; the number of iterations before the main sampling begins
#' @param iterations			 The number of iterations in MCMC, this is on top of the burn in
#' @param Clusters				 The number of clusters to be searched for
#' @param trap_rad				 The trap radius defined as some constant times the true sigma value
#' @param s_prior_shape		 The shape parameter for the prior on the population density
#' @param s_prior_rate		 The rate parameter for the prior on the population density
#' @param priorLat				 The prior mean for the source locations (Latitude)
#' @param priorLon				 The prior mean for the source locations (Longitude)
#' @param long_minMax			 The longitudinal limits of the data (already defined)
#' @param lat_minMax			 The latitudinal limits of the data (already defined)
#'
#' @export
#' @examples
#' trap_sim <- cbind(runif(10), runif(10), sample(0:5, 10, replace = T))
#' PAmcmc(trapData = trap_sim, burnin = 100, iterations = 25e2, Clusters = 1, trap_rad = 1,
#' s_prior_shape = 1, s_prior_rate = 1, priorLat = mean(trap_sim[,2]), priorLon = mean(trap_sim[,1]),
#' long_minMax = c(min(trap_sim[,1]),max(trap_sim[,1])), lat_minMax = c(min(trap_sim[,2]),max(trap_sim[,2])))

PAmcmc <- function(trapData, burnin = 100, iterations = 25e2, Clusters = 1, trap_rad = 1, s_prior_shape = 1, s_prior_rate = 1, priorLat = NULL, priorLon = NULL, long_minMax = NULL, lat_minMax = NULL)
  {
  start <- Sys.time()

  # define model parameters
  trap_density <- sum(trapData[,3])
  hit_data <- subset(trapData, trapData[,3]>0)
  sig_init <- 0.5*mean(pairwise_distance(hit_data)$distance_min, na.rm = TRUE)
  areaKm <- latlon_to_bearing(origin_lat = lat_minMax[1], origin_lon = long_minMax[1], dest_lat = lat_minMax[2], dest_lon = long_minMax[2])$gc_dist

  # define MCMC parameters
	# Robbins-Monro Step Constant
  RMcM <- rep(5, Clusters)
  RMcS <- 5.579
	# initial proposal standard deviations
  propSD_mu <- rep(areaKm, Clusters)
  PropSD_sig <- 10*areaKm

	### propose starting values
  # source locations
	old_theta <- proposal(Clusters, mu1 = 0, mu2 = 0, Sd = 2*areaKm)
  colnames(old_theta) <- c("Long", "Lat")

  # sigma
  old_sigma <- abs(rnorm(1, mean = 0, sd = PropSD_sig))

	# population density
	new_pop_den <- rgamma(1, shape = s_prior_shape, rate = s_prior_rate)

	# calculate intial joint probability
  old_posterior <- PAjoint(theta_x = old_theta[,1], theta_y = old_theta[,2], sigma = old_sigma, data = trapData, trap_rad = trap_rad, exp_pop = new_pop_den, priorLat = priorLat, priorLon = priorLon, prior = "None", meanSigprior = sig_init)

  ### create objects for storing results
  joint_lat <- joint_long <- matrix(NA, ncol = Clusters, nrow = iterations)
  joint_sig <- rep(NA, iterations)
  pop_densities <- rep(NA, iterations)
  a_r_rate_theta <- matrix(0, ncol = Clusters, nrow = iterations)
  a_r_rate_sig <- rep(0, iterations)
  propSD_mu_vals <- matrix(NA, ncol = Clusters, nrow = iterations)
  propSDsig_vals <- rep(NA, iterations)

  # run MCMC
  for (i in 1:iterations) {
    # report current iteration
    if ((i %% 500)==0) {
      cat(paste("  iteration:", i, "\n"))
    }
    for (CL in 1:Clusters) {
    # propose new theta, accept or reject!
      new_theta <- old_theta
      new_theta[CL,] <- proposal(1, mu1 = old_theta[CL,1], mu2 = old_theta[CL,2], Sd = propSD_mu[CL]) # returns two values a Lon and Lat
      # calculate new joint probability
      new_posterior <- PAjoint(theta_x = new_theta[,1], theta_y = new_theta[,2], sigma = old_sigma, data = trapData, trap_rad = trap_rad, exp_pop = new_pop_den, priorLat = priorLat, priorLon = priorLon, prior = "sourcePrior", meanSigprior = sig_init, propLoc = new_theta[CL,], priorSD = areaKm)
      # Metropolis-Hastings step
      ratio <- new_posterior$posterior_val - old_posterior$posterior_val
      temp <- sample(1:5,1)
      ratio  <- ratio/temp
      if (log(runif(1)) < ratio) {
        old_theta[CL,] <- new_theta[CL,]
        old_posterior <- new_posterior
        a_r_rate_theta[i, CL] <- 1
      }
      if(i < burnin)
      {
        theta_shift <- ( ((RMcM*0.766)^(a_r_rate_theta[i, CL]))*(-RMcM*0.234)^(1 - a_r_rate_theta[i, CL])  )/sqrt(i)
        propSD_mu[CL] <- propSD_mu[CL] + theta_shift            # Robbins Monro step for source location proposal
        propSD_mu[CL] <- max(1e-6, propSD_mu[CL])
      }
    }
    # propose new sigma, accept or reject!
    new_sigma <- abs(rnorm(1, mean = old_sigma, sd = PropSD_sig))

    # calculate new joint probability
    new_posterior <- PAjoint(theta_x = old_theta[,1], theta_y = old_theta[,2], sigma = new_sigma, data = trapData, trap_rad = trap_rad, exp_pop = new_pop_den, priorLat = priorLat, priorLon = priorLon, prior = "sigPrior", meanSigprior = sig_init, priorSD = areaKm)
    # Metropolis-Hastings step
    ratio <- new_posterior$posterior_val - old_posterior$posterior_val
    if (log(runif(1)) < ratio) {
      old_sigma <- new_sigma
      old_posterior <- new_posterior
      a_r_rate_sig[i] <- 1
    }
    if(i < burnin)
    {
      sig_shift <- ( (RMcS*0.766^(a_r_rate_sig[i]))*(-RMcS*0.234)^(1 - a_r_rate_sig[i])  )/sqrt(i)
      PropSD_sig <- PropSD_sig + sig_shift    # Robbins Monro step for sigma proposal
      PropSD_sig <- max(1e-6, PropSD_sig)
    }
    #update s by Gibbs sampling
    new_z <- PAjoint(theta_x = old_theta[,1], theta_y = old_theta[,2], sigma = old_sigma, data = trapData, trap_rad = trap_rad, exp_pop = new_pop_den, priorLat = priorLat, priorLon = priorLon, prior = "None", meanSigprior = sig_init)$z
    new_pop_den <- rgamma(1, shape = s_prior_shape + trap_density, rate = s_prior_rate + exp(logSum(new_z)))

    # # # store values of this iteration
    joint_long[i,] <- old_theta[,1]
    joint_lat[i,] <- old_theta[,2]
    joint_sig[i] <- old_sigma
    pop_densities[i] <- new_pop_den
    propSD_mu_vals[i,] <- propSD_mu
    propSDsig_vals[i] <- PropSD_sig
  }

  long_seq <- seq(long_minMax[1], long_minMax[2], (long_minMax[2] - long_minMax[1])/500)
  lat_seq <- seq(lat_minMax[1], lat_minMax[2], (lat_minMax[2] - lat_minMax[1])/500)
  rawSurface <- geoSmooth(joint_long[burnin:iterations,], joint_lat[burnin:iterations,], breaks_lon = long_seq, breaks_lat = lat_seq, lambda = NULL)

  end <- Sys.time()
  time_elapsed <- end - start
  print(time_elapsed)
  return(list(joint_long = joint_long, joint_lat = joint_lat, propSD_mu_vals = propSD_mu_vals, a_r_rate_theta = a_r_rate_theta, long_seq  = long_seq, lat_seq = lat_seq, joint_sig = joint_sig, a_r_rate_sig= a_r_rate_sig, propSDsig_vals = propSDsig_vals,
              pop_densities = pop_densities, rawSurface = rawSurface, trap_rad = trap_rad, Clusters = Clusters, popPriorSh = s_prior_shape , popPriorRa = s_prior_rate, sig_init = sig_init, areaKm = areaKm))
}

# # # install_github("Michael-Stevens-27/silverblaze", ref = "master")
# # # library(silverblaze)
# # # rm(list = ls()) #remove all objects

# pack <- "silverblaze"
# path <- find.package(pack)
# system(paste(shQuote(file.path(R.home("bin"), "R")), "CMD", "Rd2pdf", shQuote(path)))
