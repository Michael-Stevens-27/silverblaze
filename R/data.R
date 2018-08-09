
#------------------------------------------------
#' @title PAsim
#'
#' @description A function that simulates all observed and unobserved data from 
#'   the presence-absence model.
#'
#' @param n_true The underlying population size
#' @param sigma_true The underlying sigma value of the data
#' @param K_true The true number of sources
#' @param n_traps The number of traps to be used
#' @param trap_spacing Set to either "uniform", "random", or "cluster" describes
#'   the trap configuration
#' @param trap_extent When trap_spacing is clustered, this describes the distance
#'   from each trap to the centre of that cluster
#' @param trap_clusters	When trap_spacing is clustered, this is sets the number
#'   of traps within a cluster
#' @param long_minMax The longitudinal extent of the data (to be set using this
#'   function)
#' @param lat_minMax The latitudinal extent of the data (to be set using this
#'   function)
#' @param plotRail A buffer to be placed around the data when plotting
#' @param trap_const This governs the trap radius, which is the trap_const
#'   multiplied by the true sigma
#' @param single_count Set to TRUE or FALSE this governs if a n individual can
#'   be caught by multiple traps
#' @param plotting Set to TRUE or FALSE should the data want to be plotted
#' @param bias Set to 0 or 1, if set to 1 two additional sources will be
#'   created, one surrounded by empty traps, and one far away from the traps
#'
#' @export
#' @examples
#' PAsim(n_true = 100, sigma_true = 1, K_true = 3, n_traps = 10, trap_spacing = "cluster",
#' trap_extent = 0.5, long_minMax = c(0, 0.1), lat_minMax = c(0, 0.1), plotRail = 0.05,
#' trap_const = 0.5, trap_clusters = 4, single_count = T, plotting = T, bias = 1)

PAsim <- function(n_true = 1000, sigma_true = 2, K_true = 3, n_traps = 64, trap_spacing = "cluster", trap_extent = 3, long_minMax = c(0, 0.5), lat_minMax = c(0, 0.5), plotRail = 0.05, trap_const = 0.5, trap_clusters = 4, single_count = T, plotting = T, bias = 0) {
  
  # define search/area to generate data within
	total_area <- 1
  trap_rad_true <- trap_const*sigma_true
  
  # cataegorise sources in real, fake or unseen
  sTyp <- rep("Real", K_true + 2*bias)
  
	# Draw source locations from a uniform prior (i.e randomly) with or without bias
  if (bias == 1) {
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
    alloc <- sample(1:(K_true + 1), size = n_obs, replace = TRUE)
    
  	# Update the sources type object to include those extra two
  	sTyp[K_true + 2] <- "Fake"
    sTyp[K_true + 1] <- "Unseen"
    
  } else {
    
    # if there is no bias create the objects as normal
	  source_lat <- runif(K_true, lat_minMax[1], lat_minMax[2])
    source_long <- runif(K_true, long_minMax[1], long_minMax[2])
    n_obs <- rpois(1, n_true*total_area)
    alloc <- sample(1:(K_true), size = n_obs, replace = TRUE)
  }
  
	source_loc <- data.frame(source_long = source_long, source_lat = source_lat, sTyp = sTyp)
	
	# create an object to state individuals per sources
	perSource <-  mapply(function(x){length(which(x == alloc))}, x = 1:(K_true + bias))
	
	if (bias ==  1) {
  	# include the number of individuals associated with the fake source, zero.
  	perSource <-  c(perSource, 0)
	}
  
  # Draw individuals' locations from a bivariate norm with standard deviation sigma_true NEEDS EDITING
  indiv_long <- unlist(mapply(rnorm_sphere, n = perSource, centre_lat = source_lat, centre_lon = source_long, sigma = sigma_true)[1,])
  indiv_lat <- unlist(mapply(rnorm_sphere, n = perSource, centre_lat = source_lat, centre_lon = source_long, sigma = sigma_true)[2,])
	indiv_loc <- data.frame(longitude = indiv_long, latitude = indiv_lat)
  
  # set trap spacing to random, uniform, or clustered/uniform (trap_clusters dictates the number of traps in each uniform location)
  if (trap_spacing == "random") {
    trap_lat <- runif(n_traps, lat_minMax[1], lat_minMax[2] - (lat_minMax[2]/2 - lat_minMax[1]/2)*bias)
    trap_lon <- runif(n_traps, long_minMax[1], long_minMax[2])
    trap_loc <- data.frame(trap_lon = trap_lon, trap_lat = trap_lat)

  } else {
    trap_lat <- seq(lat_minMax[1], (lat_minMax[1]*(bias) + lat_minMax[2])/(2^bias), (lat_minMax[2]/(2^bias) - lat_minMax[1])/ceiling(sqrt(n_traps)))
		trap_lon <- seq(long_minMax[1], long_minMax[2], (long_minMax[2] - long_minMax[1])/ceiling(sqrt(n_traps)))
    trap_loc <- expand.grid(trap_lon, trap_lat)
  }
	
  if (trap_spacing == "cluster") {
    
    #same as method for unifrom, but then extend to produce clustered trap structure
    c_trap <- seq(2*pi/trap_clusters, 2*pi, 2*pi/trap_clusters )
		polx <- trap_extent*sin(c_trap)
		poly <- trap_extent*cos(c_trap)
		trap_lon <- mapply(function(x, y){cartesian_to_latlon(y, x, poly, polx)$longitude}, x = trap_loc[,1], y = trap_loc[,2])
		trap_lat <- mapply(function(x, y){cartesian_to_latlon(y, x, poly, polx)$latitude}, x = trap_loc[,1], y = trap_loc[,2])
    trap_loc <- cbind(c(trap_lon), c(trap_lat))
  }
	
	if (length(trap_lat) == 0 ) {
	  stop("trap_spacing must be set to 'uniform','random' or 'cluster'")
	}
	
	# allocate raw data to traps
	# rows are individuals, columns are traps, entries are distances between individuals and traps
  all_dist <- mapply(function(x, y) {
      latlon_to_bearing(trap_loc[,2], trap_loc[,1], x, y)$gc_dist
    }, x = indiv_loc$latitude, y = indiv_loc$longitude)
  all_dist <- t(all_dist) 
  
  # concatonate individuals into traps given the distance is within trap radius
  trapCounts <- colSums(all_dist <= trap_rad_true)
  
  # count the number of observed observations
  is_observed <- rowSums(all_dist <= trap_rad_true) 

  # check there are more than 5 observations overall
	if (length(which(trapCounts > 0)) < 5) {
	  stop("Not enough observations, increase sigma_true, n_traps or reduce long_minMax/lat_minMax ")
	}
  
  # allow for single or multiple trappings of the same observation
	if (single_count == T) {
	  
	  # check traps do not already only contain a single observation
    if (max(is_observed) == 1) {
      
    } else if (length(which(is_observed > 1)) > 1) { # check that there is more than one trap with multiple observations
      TF_count <- all_dist<=trap_rad_true
      for (G in 1:length(which(is_observed > 1))) {
          
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
  # rows are sources, columns are hits, entries are distances between sources and hits
  source_hits_dist <- t(mapply(function(x, y) {
      latlon_to_bearing(hit_data[,2], hit_data[,1], x, y)$gc_dist
    }, x = source_loc$source_lat, y = source_loc$source_long))
  hitsWSD <- cbind(rowSums(source_hits_dist <= sigma_true), rowSums(source_hits_dist <= 2*sigma_true), rowSums(source_hits_dist <= 3*sigma_true))
  
  # rows are sources, columns are miss, entries are distances between sources and miss
  source_miss_dist <- t(mapply(function(x, y) {
      latlon_to_bearing(miss_data[,2], miss_data[,1], x, y)$gc_dist
    }, x = source_loc$source_lat, y = source_loc$source_long))
  missWSD <- cbind(rowSums(source_miss_dist <= sigma_true) , rowSums(source_miss_dist <= 2*sigma_true), rowSums(source_miss_dist <= 3*sigma_true))
  
	# plot if plotting returns true
	if (plotting == TRUE) {
		c_seq <- seq(0, 2*pi, 0.25)
		datax <- trap_rad_true*sin(c_seq)
		datay <- trap_rad_true*cos(c_seq)
		circle_longs <- mapply(function(x, y) {
		    cartesian_to_latlon(y, x, datay, datax)$longitude
		  }, x = trap_data$longitude, y = trap_data$latitude)
		circle_lats <- mapply(function(x, y) {
		    cartesian_to_latlon(y, x, datay, datax)$latitude
		  }, x = trap_data$longitude, y = trap_data$latitude)
    
		plot(0,0, type = "n", xlim = c(long_minMax[1] - plotRail, long_minMax[2] + plotRail), ylim = c(lat_minMax[1] - plotRail, lat_minMax[2] + plotRail), xlab = "Longitude", ylab = "Latitude")
		for (i in 1:length(trap_data[,1])) {
		  polygon(circle_longs[,i], circle_lats[,i])
		}
		miss <- subset(trap_data, trap_data$count == 0)
		hits <- subset(trap_data, trap_data$count > 0)
		points(hits$longitude, hits$latitude, pch = 20, col = "green", cex = 1) # hits[,3])
		points(miss$longitude, miss$latitude, pch = 4, col = "red")
		points(indiv_loc$longitude, indiv_loc$latitude, pch = 18, cex = 0.75)
		points(source_loc$source_long, source_loc$source_lat, col = "blue", pch = 15, cex = 1.5)
	}
  
  # return list
  ret <- list(trap_data = trap_data,
              source_loc = source_loc,
              indiv_loc = indiv_loc,
              n_obs = n_obs,
              n_true = n_true,
              sigma_true = sigma_true,
              hit_data = hit_data,
              miss_data = miss_data,
              perSource = perSource,
              long_minMax = long_minMax,
              lat_minMax= lat_minMax,
              trap_spacing = trap_spacing,
              n_traps = n_traps,
              trap_rad_true = trap_rad_true,
              K_true = K_true,
              plotRail = plotRail,
              hitsWSD = hitsWSD,
              missWSD = missWSD)
	return(ret)
}
