
#------------------------------------------------
#' @title latlon_to_cartesian
#'
#' @description Convert from latitude/longitude to cartesian co-ordinates
#'
#' @param centre_lat The centre latitude
#' @param centre_lon The centre longitude
#' @param data_lat The data latitude
#' @param data_lon The data longitude
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

  return(list(x = data_x, y = data_y))
}

#------------------------------------------------
#' @title latlon_to_bearing
#'
#' @description Calculate distance between latitude/longitude points
#'
#' @param origin_lat The origin latitude
#' @param origin_lon The origin longitude
#' @param dest_lat The destination latitude
#' @param dest_lon The destination longitude
#'
#' @export
#' @examples
#' origin_latLon <- c(0,0)
#' dest_latLon <- c(1,1)
#' latlon_to_bearing(origin_lat = origin_latLon[1], origin_lon = origin_latLon[2],
#' dest_lat = dest_latLon[1], dest_lon = dest_latLon[2])
#  From B. VERITY

latlon_to_bearing <- function(origin_lat, origin_lon, dest_lat, dest_lon) {
  
  # check for exact equality of points
  if (origin_lat == dest_lat && origin_lon == dest_lon) {
    #return(list(bearing = 0, gc_dist = 0))
  }
  
  # convert input arguments to radians
  origin_lat <- origin_lat*2*pi/360
  dest_lat <- dest_lat*2*pi/360
  origin_lon <- origin_lon*2*pi/360
  dest_lon <- dest_lon*2*pi/360
  
  delta_lon <- dest_lon-origin_lon
  
  # calculate bearing
  bearing <- atan2(sin(delta_lon)*cos(dest_lat), cos(origin_lat)*sin(dest_lat)-sin(origin_lat)*cos(dest_lat)*cos(delta_lon))
  
  # calculate great circle angle. Use temporary variable to avoid acos(>1) or
  # acos(<0), which can happen due to underflow issues
  tmp <- sin(origin_lat)*sin(dest_lat) + cos(origin_lat)*cos(dest_lat)*cos(delta_lon)
  tmp <- ifelse(tmp > 1, 1, tmp)
  tmp <- ifelse(tmp < 0, 0, tmp)
  gc_angle <- acos(tmp)
  gc_angle[is.nan(gc_angle)] <- 0
  
  # convert bearing from radians to degrees measured clockwise from due north, and convert gc_angle to great circle distance via radius of earth (km)
  bearing <- bearing*360/(2*pi)
  bearing <- (bearing+360)%%360
  earthRad <- 6371
  gc_dist <- earthRad*gc_angle
  
  return(list(bearing = bearing, gc_dist = gc_dist))
}

#------------------------------------------------
#' @title bearing_to_latlon
#'
#' @description Calculate destination lat/lon given an origin, a bearing and a
#'   greater circle distance of travel
#'
#' @param origin_lat The origin latitude
#' @param origin_lon The origin longitude
#' @param bearing The angle in degrees relative to due north
#' @param gc_dist The greater circle distance in kilometres
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
#' @title cartesian_to_latlon
#'
#' @description Transform Cartesian Co-ordinates into Latitude and Longitude
#'   points assuming a Lat_Long centre
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
#' cartesian_to_latlon(centre_lat = centre[1], centre_lon = centre[2],
#'                     data_x = data[1], data_y = data[2])
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
  
  return(list(longitude = data_trans$longitude, latitude = data_trans$latitude))
}

#------------------------------------------------
#' @title Pairwise Great-Circle Distance
#'
#' @description Analogue of \code{dist()} function, but returning great-circle
#'   distances.
#'
#' @param points a matrix with two columns specifying latitude and longitude
#'
#' @export
#' @examples
#' some_points <- cbind(1:10, 1:10)
#' dist_gc(some_points)

dist_gc <- function(points) {
  
  # check inputs
  assert_matrix(points)
  assert_eq(ncol(points), 2)
  
  # calculate distance matrix
  d <- apply(points, 1, function(x) {latlon_to_bearing(x[1], x[2], points[,1], points[,2])$gc_dist})
  diag(d) <- 0
  
	return(d)
}

#------------------------------------------------
#' @title rnorm_sphere
#'
#' @description Draw from normal distribution converted to spherical coordinate
#'   system. Points are first drawn from an ordinary cartesian 2D normal
#'   distribution. The distances to points are then assumed to be great circle
#'   distances, and are combined with a random bearing from the point
#'   {centre_lat, centre_lon} to produce a final set of lat/lon points. Note
#'   that this is not a truly spherical normal distribution, as the domain of
#'   the distribution is not the sphere - rather it is a transformation from one
#'   coordinate system to another that is satisfactory when the curvature of the
#'   sphere is not severe.
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
	x <- rnorm(n, sd = sigma)
	y <- rnorm(n, sd = sigma)
	ret <- cartesian_to_latlon(centre_lat, centre_lon, x, y)
	return(ret)
}

#------------------------------------------------
#' @title logSum
#'
#' @description Sum a vector of values together in log space, while avoiding
#'   underflow issues. Takes log(x) values as input and returns log(sum(x)) as
#'   output.
#'
#' @param logx A vector of values in log space
#'
#' @export
#' @examples
#' x <- log(seq(5e-49, 5e-48, 1e-49))
#' exp(sum(x))
#' exp(logSum(x))

logSum <- function(logx) {
	max(logx) + log(sum(exp(logx - max(logx))))
}

#-------------------------------------------------------------------
#' @title bvnorm
#'
#' @description Extract the log density on the bivariate normal at (x, y)
#'   centred at (mu_x, mu_y)
#'
#' @param x The location of the point whose density we require (longitude)
#' @param y The location of the point whose density we require (latitude)
#' @param mu_x The mean of the bivariate normal (longitude)
#' @param mu_y The mean of the bivariate normal (latitude)
#' @param sd The standard deviation of the bivariate normal (kilometres)
#' @param log Whether to return output in log space
#'
#' @export
#' @examples
#' bvnorm(x = 1, y = 1, mu_x = 0, mu_y = 0, sd = 10)

bvnorm <- function(x, y, mu_x, mu_y, sd, log = TRUE) {
  dnorm(latlon_to_bearing(y, x, y, mu_x)$gc_dist, mean = 0, sd = sd, log = log) + dnorm(latlon_to_bearing(y, mu_x, mu_y, mu_x)$gc_dist, mean = 0, sd = sd, log = log)
}

#------------------------------------------------
#' @title Produce a smooth surface using 2D kernel density smoothing
#'
#' @description Takes lon/lat coordinates, bins in two dimensions, and smooths
#'   using kernel density smoothing. Kernel densities are computed using the
#'   fast Fourier transform method, which is many times faster than simple
#'   summation when using a large number of points. Each Kernel is student's-t
#'   distributed with 3 degrees of freedom, and scaled by the bandwidth lambda.
#'   If lambda is set to \code{NULL} then the optimal value of lambda is chosen
#'   automatically using the leave-one-out maximum likelihood method.
#'
#' @param longitude longitude of input points
#' @param latitude latitude of input points
#' @param breaks_lon positions of longitude breaks
#' @param breaks_lat positions of latitude breaks
#' @param lambda bandwidth to use in posterior smoothing. If NULL then optimal 
#'   bandwidth is chosen automatically by maximum-likelihood.
#'
#' @references Barnard, Etienne. "Maximum leave-one-out likelihood for kernel density estimation." Proceedings of the Twenty-First Annual Symposium of the Pattern Recognition Association of South Africa. 2010.
#' @export

geoSmooth <- function(longitude, latitude, breaks_lon, breaks_lat, lambda = NULL) {
  
  # get properties of cells in each dimension
  cells_lon <- length(breaks_lon) - 1
  cells_lat <- length(breaks_lat) - 1
  centre_lon <- mean(breaks_lon)
  centre_lat <- mean(breaks_lat)
  cellSize_lon <- diff(breaks_lon[1:2])
  cellSize_lat <- diff(breaks_lat[1:2])
  
  # bin lon/lat values in two dimensions and check that at least one value in
  # chosen region
  surface_raw <- bin2D(longitude, latitude, breaks_lon, breaks_lat)$z
  if (all(surface_raw==0)) {
    stop('chosen lat/long window contains no posterior draws')
  }
  
  # temporarily add guard rail to surface to avoid Fourier series bleeding round
  # edges
  railSize_lon <- cells_lon
  railSize_lat <- cells_lat
  railMat_lon <- matrix(0, cells_lat, railSize_lon)
  railMat_lat <- matrix(0, railSize_lat, cells_lon + 2*railSize_lon)
  
  surface_normalised <- surface_raw/sum(surface_raw)
  surface_normalised <- cbind(railMat_lon, surface_normalised, railMat_lon)
  surface_normalised <- rbind(railMat_lat, surface_normalised, railMat_lat)
  
  # calculate Fourier transform of posterior surface
  f1 = fftw2d(surface_normalised)
  
  # calculate x and y size of one cell in cartesian space. Because of
  # transformation, this size will technically be different for each cell, but
  # use centre of space to get a middling value
  cellSize_trans <- latlon_to_cartesian(centre_lat, centre_lon, centre_lat + cellSize_lat, centre_lon + cellSize_lon)
  cellSize_trans_lon <- cellSize_trans$x
  cellSize_trans_lat <- cellSize_trans$y
  
  # produce surface over which kernel will be calculated. This surface wraps
  # around in both x and y (i.e. the kernel is actually defined over a torus).
  kernel_lon <- cellSize_trans_lon * c(0:floor(ncol(surface_normalised)/2), floor((ncol(surface_normalised) - 1)/2):1)
  kernel_lat <- cellSize_trans_lat * c(0:floor(nrow(surface_normalised)/2), floor((nrow(surface_normalised) - 1)/2):1)
  kernel_lon_mat <- outer(rep(1,length(kernel_lat)), kernel_lon)
  kernel_lat_mat <- outer(kernel_lat, rep(1,length(kernel_lon)))
  kernel_s_mat <- sqrt(kernel_lon_mat^2 + kernel_lat_mat^2)
  
  # set lambda (bandwidth) range to be explored
  if (is.null(lambda)) {
    lambda_step <- min(cellSize_trans_lon, cellSize_trans_lat)/5
    lambda_vec <- lambda_step*(1:100)
  } else {
    lambda_vec <- lambda
  }
  
  # loop through range of values of lambda
  cat('Smoothing posterior surface')
  flush.console()
  logLike <- -Inf
  for (i in 1:length(lambda_vec)) {
    
    # print dots to screen
    if (i>1) {
      cat(".")
      flush.console()
    }
    
    # calculate Fourier transform of kernel
    lambda_this <- lambda_vec[i]
    kernel <- dts(kernel_s_mat, df=3, scale=lambda_this)
    f2 = fftw2d(kernel)
    
    # combine Fourier transformed surfaces and take inverse. f4 will ultimately
    # become the main surface of interest.
    f3 = f1*f2
    f4 = Re(fftw2d(f3,inverse=T))/length(surface_normalised)
    
    # subtract from f4 the probability density of each point measured from
    # itself. In other words, move towards a leave-one-out kernel density method
    f5 <- f4 - surface_normalised*dts(0, df=3, scale=lambda_this)
    f5[f5<0] <- 0
    f5 <- f5/sum(f4)
    
    # calculate leave-one-out log-likelihood at each point on surface
    f6 <- surface_normalised*log(f5)
    
    # break if total log-likelihood is at a local maximum
    if (sum(f6,na.rm=T)<logLike) {
      break()
    }
    
    # otherwise update logLike
    logLike <- sum(f6,na.rm=T)
  }
  
  # report chosen value of lambda
  if (is.null(lambda)) {
    cat(paste('\nmaximum likelihood lambda = ', round(lambda_this,3), sep=''))
  }
  
  # remove guard rail
  f4 <- f4[,(railSize_lon+1):(ncol(f4)-railSize_lon)]
  f4 <- f4[(railSize_lat+1):(nrow(f4)-railSize_lat),]
  
  # return surface
  return(f4)
}

#------------------------------------------------
# Bin values in two dimensions
# (not exported)
#' @noRd
bin2D <- function(x, y, x_breaks, y_breaks) {
  
  # get number of breaks in each dimension
  nx <- length(x_breaks)
  ny <- length(y_breaks)
  
  # create table of binned values
  tab1 <- table(findInterval(x, x_breaks), findInterval(y, y_breaks))
  
  # convert to dataframe and force numeric
  df1 <- as.data.frame(tab1, stringsAsFactors=FALSE)
  names(df1) <- c("x", "y", "count")
  df1$x <- as.numeric(df1$x)
  df1$y <- as.numeric(df1$y)
  
  # subset to within breaks range
  df2 <- subset(df1, x>0 & x<nx & y>0 & y<ny)
  
  # fill in matrix
  mat1 <- matrix(0,ny-1,nx-1)
  mat1[cbind(df2$y, df2$x)] <- df2$count
  
  # calculate cell midpoints
  x_mids <- (x_breaks[-1]+x_breaks[-nx])/2
  y_mids <- (y_breaks[-1]+y_breaks[-ny])/2
  
  # return output as list
  return(list(x_mids = x_mids,
              y_mids = y_mids,
              z = mat1))
}

#------------------------------------------------
# Scaled Student's t distribution. Used in kernel density smoothing.
# (not exported)
#' @noRd
dts <- function(x, df, scale=1, log = FALSE) {
  ret <- lgamma((df+1)/2)-lgamma(df/2)-0.5*log(pi*df*scale^2) - ((df+1)/2)*log(1 + x^2/(df*scale^2))
  if (!log) {
    ret <- exp(ret)
  }
  return(ret)
}