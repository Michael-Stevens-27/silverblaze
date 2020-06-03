
#------------------------------------------------
#' @title Simulate data
#'
#' @description Simulate data from the same presence-absence model used in the
#'   inference step.
#'
#' @param sentinel_lon vector giving longitudes of sentinel sites.
#' @param sentinel_lat vector giving latitudes of sentinel sites.
#' @param sentinel_radius observation radius of the sentinel site (km).
#' @param K the number of sources.
#' @param source_lon_min minimum limit on source longitudes.
#' @param source_lon_max maximum limit on source longitudes.
#' @param source_lat_min minimum limit on source latitudes.
#' @param source_lat_max maximum limit on source latitudes.
#' @param source_lon manually define source longitude positions. If \code{NULL}
#'   then drawn uniformly from limits specified in \code{source_lon_min} and
#'   \code{source_lon_max}.
#' @param source_lat manually define source latitude positions. If \code{NULL}
#'   then drawn uniformly from limits specified in \code{source_lat_min} and
#'   \code{source_lat_max}.
#' @param sigma_model set as "single" to use the same dispersal distance for all
#'   sources, or "separate" to use an independently drawn dispersal distance for
#'   each source.
#' @param sigma_mean the prior mean of the parameter sigma (km).
#' @param sigma_var the prior variance of the parameter sigma (km). Set to zero
#'   to use a fixed distance.
#' @param expected_popsize the expected total number of observations (observed
#'   and unobserved) in the study area.
#' @param data_type what model we wish to simulate under - a poisson or binomial
#'   corresponding to "counts" or "prevalence"
#' @param test_rate The rate of the Poisson distribution with which we draw the 
#'   number of individuals tested at each sentinel site
#'
#' @import stats
#' @export
#' 
#' @examples
#' # State the number of sources to be generated
#' K_sim <- 3
#' # Create some sentinel site locations
#' sentinal_lon <- seq(-0.2, 0.0, l=11)
#' sentinal_lat <- seq(51.45, 51.55, l=11)
#' sentinal_grid <- expand.grid(sentinal_lon, sentinal_lat)
#' names(sentinal_grid) <- c("longitude", "latitude")
#' # Set their sentinel radius (this constant times true sigma)
#' sentinel_radius <- 0.25
#' sim1 <- sim_data(sentinal_grid$longitude,
#'                 sentinal_grid$latitude,
#'                 sigma_model = "single",
#'                 sigma_mean = 1,
#'                 sigma_var = 0.5,
#'                 sentinel_radius = sentinel_radius,
#'                 K = K_sim,
#'                 expected_popsize = 300)

sim_data <- function(sentinel_lon,
                     sentinel_lat,
                     sentinel_radius = 0.1,
                     K = 3,
                     source_lon_min = -0.2,
                     source_lon_max = 0.0,
                     source_lat_min = 51.45,
                     source_lat_max = 51.55,
                     source_lon = NULL,
                     source_lat = NULL,
                     sigma_model = "single",
                     sigma_mean = 1.0,
                     sigma_var = 0.1,
                     expected_popsize = 100,
                     data_type = "counts",
                     test_rate = 5)
                     {

  # check inputs
  assert_numeric(sentinel_lon)
  assert_numeric(sentinel_lat)
  assert_single_pos(sentinel_radius, zero_allowed = FALSE)
  assert_single_pos_int(K, zero_allowed = FALSE)
  assert_single_numeric(source_lon_min)
  assert_single_numeric(source_lon_max)
  assert_single_numeric(source_lat_min)
  assert_single_numeric(source_lat_max)
  if (is.null(source_lon)) {
    source_lon <- runif(K, source_lon_min, source_lon_max)
  }
  if (is.null(source_lat)) {
    source_lat <- runif(K, source_lat_min, source_lat_max)
  }
  assert_vector(source_lon)
  assert_numeric(source_lon)
  assert_vector(source_lat)
  assert_numeric(source_lat)
  assert_same_length(source_lon, source_lat)
  assert_length(source_lon, K)
  assert_same_length(sentinel_lon, sentinel_lat)
  assert_single_string(sigma_model)
  assert_in(sigma_model, c("single", "independent"))
  switch(sigma_model,
         "single" = {
          assert_single_pos(sigma_var, zero_allowed = TRUE)
          assert_single_pos(sigma_mean, zero_allowed = FALSE)
         },
         "independent" = {
           assert_length(sigma_mean, K)
           assert_length(sigma_var, K)
           assert_pos(sigma_mean, zero_allowed = FALSE)
           assert_pos(sigma_var, zero_allowed = TRUE)
         })
  assert_single_pos(expected_popsize, zero_allowed = FALSE)
  assert_in(data_type, c("counts", "prevalence"))
  assert_single_pos_int(test_rate, zero_allowed = FALSE)
  
  # draw total number of points
  N <- rpois(1, expected_popsize)
  if (N == 0) {
    stop("N=0 events generated")
  }
  
  # draw true allocation of all points to sources
  group <- sort(sample(K, N, replace = TRUE))
  source_N <- tabulate(group)
  
  # draw sigma
  varlog <- log(sigma_var/sigma_mean^2 + 1)
  meanlog <- log(sigma_mean) - varlog/2
  switch(sigma_model,
         "single" = {
           sigma <- rep(rlnorm(1, meanlog, sqrt(varlog)), K)
         },
         "independent" = {
           sigma <- rlnorm(K, meanlog, sqrt(varlog))
         })
  
  #-----------------------------------------------------------------------------
  if(data_type == "counts"){
    df_all <- NULL
    for (k in 1:K) {
      if (source_N[k]>0) {
        rand_k <- rnorm_sphere(source_N[k], source_lon[k], source_lat[k], sigma[k])
        df_all <- rbind(df_all, as.data.frame(rand_k))
      }
    }
    
    # draw points around sources
    # get distance between all points and sentinel sites
    gc_dist <- mapply(function(x, y) {
      lonlat_to_bearing(x, y, df_all$longitude, df_all$latitude)$gc_dist
    }, x = sentinel_lon, y = sentinel_lat)
    
    counts <- colSums(gc_dist < sentinel_radius)
    
    df_observed <- data.frame(longitude = sentinel_lon,
                              latitude = sentinel_lat,
                              counts = counts)
                          
    # add record of whether data point is observed or unobserved to df_all
    df_all$observed <- rowSums(gc_dist < sentinel_radius)
    observed_by <- as.list(apply(gc_dist, 1, function(x) which(x < sentinel_radius)))
    if (length(observed_by) == 0) {
      observed_by <- replicate(nrow(df_all), integer())
    }
    df_all$observed_by <- observed_by
    
    # create true q-matrix as proportion of points belonging to each group per sentinel site
    true_qmatrix <- t(apply(gc_dist, 2, function(x) {
      ret <- tabulate(group[x < sentinel_radius], nbins = K)
      ret <- ret/sum(ret)
      ret[is.na(ret)] <- NA
      ret
    }))
    class(true_qmatrix) <- "rgeoprofile_qmatrix"
    
    # return simulated data and true parameter values
    ret_data <- df_observed
    ret_record <- list()
    ret_record$sentinel_radius <- sentinel_radius
    ret_record$true_group <- group
    ret_record$true_qmatrix <- true_qmatrix
    ret_record$data_all <- df_all
    
    } else if (data_type == "prevalence"){
    
    # get distances from source locations to sentinel sites
    gc_dist <- mapply(function(x, y) {lonlat_to_bearing(x, y, source_lon, source_lat)$gc_dist},
                                      x = sentinel_lon, y = sentinel_lat)
    
    # calculate mean height of each sentinel site on the mixture of bivariate normals
    heights <- dnorm(gc_dist, 0, sigma)*dnorm(0, 0, sigma)
    if(is.vector(heights)){
      av_heights <- heights       
    } else{
      av_heights <- apply(heights, 2, mean)
    }
    rate <- expected_popsize*av_heights
    
    # transform the rate to a trial success probability
    binom_prob <- rate/(1 + rate)
    
    # pick how many individuals are tested at each site and use the binom_prob to draw 
    # the number of positive individuals    
    tested <- rpois(length(sentinel_lon), lambda = test_rate)
    tested[tested == 0] <- 1
    positive <- rbinom(length(sentinel_lon), tested, binom_prob)
    
    df_observed <- data.frame(longitude = sentinel_lon,
                              latitude = sentinel_lat,
                              positive = positive,
                              tested = tested)
                                  
    # return simulated data and true parameter values
    ret_data <- df_observed
    ret_record <- list()
    ret_record$binomial_probability <- binom_prob 
  }
    
  ret_record$true_source <- data.frame(longitude = source_lon, latitude = source_lat)
  ret_record$true_source_N <- source_N
  ret_record$true_sigma <- sigma
  
  ret <- list(data = ret_data,
              record = ret_record)
  
  # make custom class
  class(ret) <- "rgeoprofile_simdata"
  
  return(ret)
}
