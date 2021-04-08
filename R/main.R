
#------------------------------------------------
#' @title Check that silverblaze package has loaded successfully
#'
#' @description Simple function to check that silverblaze package has loaded
#'   successfully. Prints "silverblaze loaded successfully!" if so.
#'
#' @export

check_silverblaze_loaded <- function() {
  message("silverblaze loaded successfully!")
}

#------------------------------------------------
#' @title Bind data to project
#'
#' @description Load data into a \code{rgeoprofile_project} prior to analysis.
#'   Data must be formatted as a dataframe with the formatting described bellow.
#'
#' @param project an \code{rgeoprofile_project}, as produced by the function
#'   \code{rgeoprofile_project()}
#' @param df a dataframe with columns that must conform to the following rules:
#'   \itemize{
#'     \item for \code{data_type = "counts"}, data must have columns
#'     "longitude", "latitude" and "counts".
#'     \item for \code{data_type = "prevalence"}, data must have columns
#'     "longitude", "latitude", "tested" and "positive"
#'     \item for \code{data_type = "point-pattern"}, data must have columns
#'     "longitude", "latitude"
#'     }
#' @param data_type the type of data, either "counts", "prevalence" or "point-pattern"
#' @param name optional name of the data set to aid in record keeping
#' @param check_delete_output whether to prompt the user before overwriting
#'   existing data
#'
#' @export

bind_data <- function(project,
                      df,
                      data_type,
                      name = NULL,
                      check_delete_output = TRUE) {

  # check inputs
  assert_custom_class(project, "rgeoprofile_project")
  assert_dataframe(df)
  assert_single_string(data_type)
  assert_in(data_type, c("counts", "prevalence", "point-pattern"))
  
  if (data_type == "counts") {
    assert_in(c("longitude", "latitude", "counts"), names(df))
    assert_pos_int(df$counts, zero_allowed = TRUE)
  } else if (data_type == "prevalence") {
    assert_in(c("longitude", "latitude", "tested", "positive"), names(df))
    assert_pos_int(df$tested, zero_allowed = TRUE)
    assert_pos_int(df$positive, zero_allowed = TRUE)
    assert_leq(df$positive, df$tested)
  } else if(data_type == "point-pattern"){
    assert_in(c("longitude", "latitude"), names(df))
  }
  
  assert_numeric(df$longitude)
  assert_numeric(df$latitude)
  if (!is.null(name)) {
    assert_single_string(name)
  }
  assert_single_logical(check_delete_output)
  
  # check before overwriting existing output
  if (project$active_set > 0 && check_delete_output) {

    # ask before overwriting. On abort, return original project
    if (!user_yes_no("All existing output and parameter sets for this project will be lost. Continue? (Y/N): ")) {
      return(project)
    }

    # replace old project with fresh empty version
    project <- rgeoprofile_project()
  }
  
  # add data to project
  project$data <- list(frame = df,
                       data_type = data_type)
  
  return(project)
}

#------------------------------------------------
#' @title Make raster grid
#'
#' @description Make raster grid
#'
#' @param range_lon  min and max longitude
#' @param range_lat  min and max latitude
#' @param cells_lon  number of cells in longitude direction
#' @param cells_lat  number of cells in latitude direction
#' @param guard_rail Extend each lon lat range by a proportion of the length of said range. 
#'    E.g. a guard_rail of 0.05 increases the lon and lat range by 5 percent.
#'
#' @importFrom raster raster setValues
#' @export

raster_grid <- function (range_lon = c(-0.2, 0),
                         range_lat = c(51.45, 51.55),
                         cells_lon = 1e2,
                         cells_lat = 1e2,
                         guard_rail = 0.05) {
  
  # check inputs
  assert_numeric(range_lon)
  assert_vector(range_lon)
  assert_length(range_lon, 2)
  assert_numeric(range_lat)
  assert_vector(range_lat)
  assert_length(range_lat, 2)
  assert_single_pos_int(cells_lon)
  assert_single_pos_int(cells_lat)
  assert_numeric(guard_rail)
  assert_single_pos(guard_rail)
                      
  # define raster extent
  dlon <- guard_rail*diff(range(range_lon))
  dlat <- guard_rail*diff(range(range_lat))
  lomn <- range_lon[1] - dlon/2
  lomx <- range_lon[2] + dlon/2
  lamn <- range_lat[1] - dlat/2
  lamx <- range_lat[2] + dlat/2
  
  # create blank raster with this extent
  r <- raster::raster(xmn = lomn,
                      xmx = lomx,
                      ymn = lamn,
                      ymx = lamx,
                      ncol = cells_lon,
                      nrow = cells_lat)
  r <- raster::setValues(r, 1/(cells_lon*cells_lat))
  return(r)
}

#------------------------------------------------
#' @title Make raster from shapefile
#'
#' @description Make raster from shapefile
#'
#' @param shp shapefile to convert to raster
#' @param cells_lon number of cells in longitude direction
#' @param cells_lat number of cells in latitude direction
#'
#' @importFrom raster raster extent rasterize projectRaster values
#' @export

raster_from_shapefile <- function (shp,
                                   cells_lon = 1e2,
                                   cells_lat = 1e2) {

  # check inputs
  assert_in(class(shp), c("SpatialPolygonsDataFrame","SpatialLinesDataFrame"))

  # make raster from shapefile
  r <- raster::raster(ncol = cells_lon, nrow = cells_lat)
  raster::extent(r) <- raster::extent(shp)
  r <- raster::rasterize(shp, r)
  r <- raster::projectRaster(r, crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

  # set all non-NA values to 1 over number of non-NA cells
  raster::values(r)[!is.na(values(r))] <- 1/sum(!is.na(values(r)))

  return(r)
}

#------------------------------------------------
#' @title Create new parameter set
#'
#' @description Create a new parameter set within an \code{rgeoprofile_project}. The new
#'   parameter set becomes the active set once created.
#'
#' @param project an rgeoprofile_project, as produced by the function
#'   \code{rgeoprofile_project()}.
#' @param name an optional name for the parameter set.
#' @param spatial_prior a raster file defining the spatial prior. Precision
#'   values are taken from this raster if it is defined.
#' @param source_model choose prior type for source locations. Pick from "uniform"
#'   (default), "normal" (bivariate normal), "kernal" (KDE based on positive data) or
#'   "manual" (the current value of the raster) 
#' @param dispersal_model distribute points via a "normal", "cauchy" or 
#'   "laplace" model
#' @param sigma_model set as \code{"single"} to assume the same dispersal
#'   distance for all sources, or \code{"independent"} to assume an
#'   independently drawn dispersal distance for each source.
#' @param sigma_prior_mean the prior mean of the parameter sigma (km).
#' @param sigma_prior_sd the prior standard deviation of the parameter sigma
#'   (km). Set to 0 to use a fixed value for sigma (fixed at
#'   \code{sigma_prior_mean}).
#' @param expected_popsize_model set as \code{"single"} to assume the same
#'   number of events for all sources, or \code{"independent"} to assume an
#'   independently drawn number of events for each source.
#' @param expected_popsize_prior_mean the prior mean of the expected total
#'   population size.
#' @param expected_popsize_prior_sd the prior standard deviation of the expected
#'   total population size. Set to 0 to use a fixed value (fixed at
#'   \code{expected_popsize_prior_mean}).
#' @param sentinel_radius the observation radius of sentinel sites.
#' @param n_binom set to true or false, decide if a negative binomial model should
#'   be run for a set of over-dispersed count data.
#' @param alpha_prior_mean the prior mean alpha.
#' @param alpha_prior_sd the prior standard deviation of alpha.
#' @param weight_prior control the prior on weights for a point-pattern model
#'
#' @export

new_set <- function(project,
                    spatial_prior = NULL,
                    source_model = "uniform",
                    name = "(no name)",
                    sigma_model = "single",
                    dispersal_model = "normal",
                    sigma_prior_mean = 1,
                    sigma_prior_sd = 1,
                    expected_popsize_model = "single",
                    expected_popsize_prior_mean = 1000,
                    expected_popsize_prior_sd = 20,
                    sentinel_radius = 0.2,
                    n_binom = FALSE,
                    alpha_prior_mean = 1, 
                    alpha_prior_sd = 100,
                    weight_prior = 1) {
  
  # check inputs
  assert_custom_class(project, "rgeoprofile_project")
  
  if (!is.null(spatial_prior)) {
    assert_custom_class(spatial_prior, "RasterLayer")
  }
  assert_in(source_model, c("uniform", "normal", "kernel", "manual"))
  
  assert_in(dispersal_model, c("normal", "cauchy", "laplace"))
  assert_in(sigma_model, c("single", "independent"))
  assert_single_pos(sigma_prior_mean, zero_allowed = FALSE)
  assert_single_pos(sigma_prior_sd, zero_allowed = TRUE)
  assert_single_pos(sentinel_radius, zero_allowed = FALSE)
  
  if (project$data$data_type == "counts" | project$data$data_type == "prevalence") {
    assert_in(expected_popsize_model, c("single", "independent"))
    assert_single_pos(expected_popsize_prior_mean, zero_allowed = FALSE)
    assert_single_pos(expected_popsize_prior_sd, zero_allowed = TRUE)
    assert_single_logical(n_binom)
    
    if (n_binom == TRUE) { # negative binomial model
      assert_single_pos(alpha_prior_mean, zero_allowed = FALSE)
      assert_single_pos(alpha_prior_sd, zero_allowed = TRUE)
    }
  } else if(project$data$data_type == "point-pattern"){
    assert_single_pos(weight_prior, zero_allowed = FALSE)
  }
  
  assert_single_string(name)
  
  # make uniform spatial_prior from data limits if unspecified
  if (is.null(spatial_prior)) {
    range_lon <- range(project$data$frame$longitude)
    range_lon <- mean(range_lon) + 1.05*c(-1,1)*diff(range_lon)/2
    range_lat <- range(project$data$frame$latitude)
    range_lat <- mean(range_lat) + 1.05*c(-1,1)*diff(range_lat)/2
    spatial_prior <- raster_grid(range_lon, range_lat, cells_lon = 1e2, cells_lat = 1e2)
    values(spatial_prior) <- 1/(1e2^2)
    
  } else if (source_model == "uniform"){
    # create uniform prior based on specified raster
    values(spatial_prior)[!is.na(values(spatial_prior))] <- 1/sum(!is.na(values(spatial_prior)))
    values(spatial_prior)[is.na(values(spatial_prior))] <- 0
  
  } else if (source_model == "normal"){
    # create spatial prior based on mean of positive data
    if (project$data$data_type == "counts") {
      pos_df <- subset(project$data$frame, project$data$frame$counts > 0)
    } else if (project$data$data_type == "prevalence") {
      pos_df <- subset(project$data$frame, project$data$frame$positive > 0)
    } else if (project$data$data_type == "point-pattern") {
      pos_df <- project$data$frame
    }
      
    # mean of positive points
    source_prior_mean_lon <- mean(pos_df$longitude)
    source_prior_mean_lat <- mean(pos_df$latitude)
      
    # max distance between POSITIVE data points and source prior mean 
    # for source prior sd 
    source_prior_sd <- apply(pos_df[,1:2], 1, function(x) 
                             lonlat_to_bearing(source_prior_mean_lon,
                                               source_prior_mean_lat, 
                                               x[1], x[2])$gc_dist)
    source_prior_sd <- max(source_prior_sd)
    
    # create domain based on empty sptail_prior raster 
    lomn <- extent(spatial_prior)[1]
    lomx <- extent(spatial_prior)[2]
    lamn <- extent(spatial_prior)[3]
    lamx <- extent(spatial_prior)[4]
    
    lon_dom <- seq(lomn, lomx, l = ncol(spatial_prior))
    lat_dom <- seq(lamn, lamx, l = nrow(spatial_prior))
    grid_domain <- expand.grid(lon_dom, lat_dom)
            
    source_prior_dists <- apply(grid_domain, 1, function(x) 
                                lonlat_to_bearing(source_prior_mean_lon,
                                                  source_prior_mean_lat, 
                                                  x[1], x[2])$gc_dist)
    
    # calculate density of the bivariate normal for each cell  
    source_prior_vals <- dnorm(source_prior_dists, 0, sd = source_prior_sd)*dnorm(0, 0, sd = source_prior_sd)
    source_prior_vals <- matrix(source_prior_vals, 
                                ncol = ncol(spatial_prior), 
                                nrow = nrow(spatial_prior), 
                                byrow = T)
    source_prior_vals <- apply(source_prior_vals, 2, rev)
    
    # allocate these values to the spatial prior (masking out areas with NAs) 
    values(spatial_prior)[!is.na(values(spatial_prior))] <- 1
    values(spatial_prior) <- values(spatial_prior)*source_prior_vals
    values(spatial_prior) <- values(spatial_prior)/sum(values(spatial_prior), na.rm = TRUE)
    
  } else if (source_model == "kernel"){
    # create spatial prior based on a KDE of positive data
    if (project$data$data_type == "counts") {
      pos_df <- subset(project$data$frame, project$data$frame$counts > 0)
    } else if (project$data$data_type == "prevalence") {
      pos_df <- subset(project$data$frame, project$data$frame$positive > 0)
    } else if (project$data$data_type == "point-pattern") {
      pos_df <- project$data$frame
    }
  
    # create domain based on empty spatial_prior raster 
    lomn <- extent(spatial_prior)[1]
    lomx <- extent(spatial_prior)[2]
    lamn <- extent(spatial_prior)[3]
    lamx <- extent(spatial_prior)[4]
    lon_dom <- seq(lomn, lomx, l = ncol(spatial_prior) + 1)
    lat_dom <- seq(lamn, lamx, l = nrow(spatial_prior) + 1)
    
    # run kernel density estimator on positive data
    source_prior_vals <- kernel_smooth(pos_df$longitude,
                                       pos_df$latitude,
                                       lon_dom,
                                       lat_dom)
    source_prior_vals <- t(apply(source_prior_vals, 2, rev))
      
    # allocate these values to the spatial prior (masking out areas with NAs) 
    values(spatial_prior)[!is.na(values(spatial_prior))] <- 1
    values(spatial_prior) <- values(spatial_prior)*c(source_prior_vals)
    values(spatial_prior) <- values(spatial_prior)/sum(values(spatial_prior), na.rm = TRUE)
  
  } else if (source_model == "manual"){
    print("Using manually specified values in spatial prior")
  }

  # get average single cell area and total study area in km^2
  study_area <- sum(raster::area(spatial_prior)[])
  cell_area <- mean(raster::area(spatial_prior)[])
  
  # count current parameter sets and add one
  s <- length(project$parameter_sets) + 1
  
  # make new set active
  project$active_set <- s
  
  # create new parameter set
  project$parameter_sets[[s]] <- list(name = name,
                                      spatial_prior = spatial_prior,
                                      study_area = study_area,
                                      cell_area = cell_area, 
                                      dispersal_model = dispersal_model, 
                                      sentinel_radius = sentinel_radius,
                                      sigma_model = sigma_model,
                                      sigma_prior_mean = sigma_prior_mean,
                                      sigma_prior_sd = sigma_prior_sd,
                                      expected_popsize_model = expected_popsize_model,
                                      expected_popsize_prior_mean = expected_popsize_prior_mean,
                                      expected_popsize_prior_sd = expected_popsize_prior_sd,
                                      n_binom = n_binom,
                                      alpha_prior_mean = alpha_prior_mean,
                                      alpha_prior_sd = alpha_prior_sd,
                                      weight_prior = weight_prior)

  # name parameter set
  names(project$parameter_sets)[s] <- paste0("set", s)
  
  # create new output at all_K level
  project$output$single_set[[s]] <- list(single_K = list(), all_K = list())
  
  # name new output
  names(project$output$single_set) <- paste0("set", 1:length(project$output$single_set))
  
  # return
  return(project)
}

#------------------------------------------------
#' @title Delete parameter set
#'
#' @description Delete a given parameter set from an \code{rgeoprofile_project}.
#'
#' @param project an rgeoprofile_project, as produced by the function
#'   \code{rgeoprofile_project()}
#' @param set which set to delete. Defaults to the current active set
#' @param check_delete_output whether to prompt the user before deleting any
#'   existing output
#'
#' @export

delete_set <- function(project,
                       set = NULL,
                       check_delete_output = TRUE) {

  # check inputs
  assert_custom_class(project, "rgeoprofile_project")
  assert_single_logical(check_delete_output)

  # set index to active_set by default
  set <- define_default(set, project$active_set)

  # further checks
  assert_single_pos_int(set, zero_allowed = FALSE)
  assert_leq(set, length(project$parameter_sets))

  # check before overwriting existing output
  if (project$active_set>0 & check_delete_output) {

    # ask before overwriting. On abort, return original project
    if (!user_yes_no(sprintf("Any existing output for set %s will be deleted. Continue? (Y/N): ", set))) {
      return(project)
    }
  }

  # drop chosen parameter set
  project$parameter_sets[[set]] <- NULL

  # drop chosen output
  project$output$single_set[[set]] <- NULL

  # make new final set active
  project$active_set <- length(project$parameter_sets)

  # return
  return(project)
}

#------------------------------------------------
#' @title Run main MCMC
#'
#' @description Run the main geographc profiling MCMC. Model parameters are taken 
#'   from the current active parameter set, and MCMC parameters are passed in as
#'   arguments. All output is stored within the project.
#'
#' @details Both longitude and latitude values can be represented to a given
#'   precision level using the arguments \code{precision_lon} and
#'   \code{precision_lat} - for example, a precision of 0.01 means that values
#'   are rounded to the nearest hundredth of a degree. This allows the use of
#'   look-up tables for the likelihood calculation, which significantly speeds
#'   up the MCMC. Set to 0 to use exact values (up to C++ "double" precision)
#'   rather than using look-up tables.
#'
#' @param project an rgeoprofile_project, as produced by the function
#'   \code{rgeoprofile_project()}.
#' @param K the number of sources.
#' @param burnin the number of burn-in iterations.
#' @param samples the number of sampling iterations.
#' @param auto_converge whether convergence should be assessed automatically
#'   every \code{converge_test} iterations, leading to termination of the
#'   burn-in phase. If \code{FALSE} then the full \code{burnin} iterations are
#'   used.
#' @param converge_test test for convergence every \code{convergence_test}
#'   iterations if \code{auto_converge} is being used.
#' @param coupling_on whether to implement Metropolis coupling.
#' @param beta_manual allows manual specification of thermodynamic
#'   powers used. 
#' @param cluster option to pass in a cluster environment (see package
#'   "parallel").
#' @param pb_markdown whether to run progress bars in markdown mode, in which
#'   case they are updated once at the end to avoid large amounts of output.
#' @param store_raw whether to store raw MCMC output in addition to summary
#'   output. Setting to FALSE can considerably reduce output size in memory.
#' @param create_maps whether to create maps of posterior probability and
#'   geoprofile. Usually will want to create these maps, but the code runs much
#'   faster without this step, hence the option.
#' @param silent whether to suppress all console output.
#' @param rung_store Pick a rung whose output will be stored 
#'
#' @import parallel
#' @import coda
#' @import stats
#' @importFrom fftwtools fftw2d
#' @importFrom utils txtProgressBar
#' @importFrom raster values<- values xyFromCell xmin xmax xres ymin ymax yres extent<- extent res<- res addLayer
#'
#' @export

run_mcmc <- function(project,
                     K = 3,
                     burnin = 1e2,
                     samples = 1e3,
                     auto_converge = TRUE,
                     converge_test = 1e2,
                     coupling_on = FALSE,
                     cluster = NULL,
                     pb_markdown = FALSE,
                     store_raw = TRUE,
                     create_maps = TRUE,
                     silent = !is.null(cluster),
                     beta_manual = NULL,
                     rung_store = NULL) {
  
  # start timer
  t0 <- Sys.time()
  
  if(is.null(rung_store)){
    if(is.null(beta_manual)){
      rung_store <- rungs <- 1
      beta_manual <- 1
    } else {
      rung_store <- rungs <- length(beta_manual)
    }
  } else {
    assert_pos_int(rung_store, zero_allowed = FALSE)
    assert_leq(rung_store, length(beta_manual))
    rungs <- length(beta_manual)
  }

  # check inputs
  assert_custom_class(project, "rgeoprofile_project")
  assert_pos_int(K, zero_allowed = FALSE)
  assert_single_pos_int(burnin, zero_allowed = FALSE)
  assert_single_pos_int(samples, zero_allowed = FALSE)
  assert_greq(samples, 10, message = "at least 10 sampling iterations must be used")
  assert_single_pos_int(rungs, zero_allowed = FALSE)
  assert_single_logical(auto_converge)
  assert_single_pos_int(converge_test, zero_allowed = FALSE)
  assert_single_logical(coupling_on)
  
  if (!is.null(beta_manual)) {
    assert_vector(beta_manual)
    assert_bounded(beta_manual)
    assert_eq(beta_manual, sort(beta_manual),
              message = "beta_manual must be increasing from left to right")
    assert_eq(beta_manual[length(beta_manual)], 1.0,
              message = "final value of beta_manual (i.e. cold chain) must equal 1.0")
  }
  
  if (!is.null(cluster)) {
    assert_custom_class(cluster, "cluster")
  }
  assert_single_logical(pb_markdown)
  assert_single_logical(silent)
  
  beta_vec <- beta_manual
  rungs <- length(beta_vec)
  
  # get active set
  s <- project$active_set
  if (s == 0) {
    stop("no active parameter set")
  }
  
  # ---------- create argument lists ----------
  
  # data list depending on data_type
  if (project$data$data_type == "counts") {
    args_data <- list(longitude = project$data$frame$longitude,
                      latitude = project$data$frame$latitude,
                      counts = project$data$frame$counts,
                      tested = -1,
                      positive = -1,
                      data_type = 1)
  } else if (project$data$data_type == "prevalence") {
    args_data <- list(longitude = project$data$frame$longitude,
                      latitude = project$data$frame$latitude,
                      counts = -1,
                      tested = project$data$frame$tested,
                      positive = project$data$frame$positive,
                      data_type = 2)
  } else if (project$data$data_type == "point-pattern") {
    args_data <- list(longitude = project$data$frame$longitude,
                      latitude = project$data$frame$latitude,
                      counts = -1,
                      tested = -1,
                      positive = -1,
                      data_type = 3)
  }
  
  # input arguments list
  args_inputs <- list(burnin = burnin,
                      samples = samples,
                      beta_vec = beta_vec,
                      auto_converge = auto_converge,
                      converge_test = converge_test,
                      coupling_on = coupling_on,
                      pb_markdown = pb_markdown,
                      silent = silent,
                      rung_store = rung_store)
  
  # extract spatial prior object
  spatial_prior <- project$parameter_sets[[s]]$spatial_prior
  spatial_prior_values <- log(raster::values(spatial_prior))
  spatial_prior_values[is.na(spatial_prior_values)] <- 0
  
  # initialise sources based on prior
  source_init <- raster::xyFromCell(spatial_prior, sample(x = 1:ncell(spatial_prior), 
                                                          size = max(K), 
                                                          prob = values(spatial_prior), 
                                                          replace = TRUE))
  source_init_lon <- source_init[,1]                                                          
  source_init_lat <- source_init[,2]  
                                                        
  # convert sigma_model to numeric
  sigma_model_numeric <- match(project$parameter_sets[[s]]$sigma_model, c("single", "independent"))
  fixed_sigma_model <- project$parameter_sets[[s]]$sigma_prior_sd == 0
  
  # convert expected_popsize_model to numeric
  if (project$data$data_type == "counts" | project$data$data_type == "prevalence") {
    expected_popsize_model_numeric <- match(project$parameter_sets[[s]]$expected_popsize_model, c("single", "independent"))
  } else if (project$data$data_type == "point-pattern") {
    expected_popsize_model_numeric <- 2
  }
  
  # convert dispersal model and count type to numeric
  dispersal_model_numeric <- match(project$parameter_sets[[s]]$dispersal_model, c("normal", "cauchy", "laplace"))
  count_type_numeric <- match(project$parameter_sets[[s]]$n_binom, c(TRUE, FALSE))
  
  # misc properties list
  args_properties <- list(min_lon = raster::xmin(spatial_prior),
                          max_lon = raster::xmax(spatial_prior),
                          res_lon = raster::xres(spatial_prior),
                          n_lon = ncol(spatial_prior),
                          min_lat = raster::ymin(spatial_prior),
                          max_lat = raster::ymax(spatial_prior),
                          res_lat = raster::yres(spatial_prior),
                          n_lat = nrow(spatial_prior),
                          spatial_prior_values = spatial_prior_values,
                          source_init_lon = source_init_lon,
                          source_init_lat = source_init_lat,
                          sigma_model_numeric = sigma_model_numeric,
                          fixed_sigma_model = fixed_sigma_model,
                          expected_popsize_model_numeric = expected_popsize_model_numeric,
                          dispersal_model_numeric = dispersal_model_numeric,
                          count_type_numeric = count_type_numeric)
  
  # combine parameters, inputs and properties into single list
  args_model <- c(project$parameter_sets[[s]], args_inputs, args_properties)
  
  # R functions to pass to Rcpp
  args_functions <- list(test_convergence = test_convergence,
                         update_progress = update_progress)
  
  # define final argument list over all K
  parallel_args <- list()
  for (i in 1:length(K)) {
    
    # create progress bars
    pb_burnin <- txtProgressBar(min = 0, max = burnin, initial = NA, style = 3)
    pb_samples <- txtProgressBar(min = 0, max = samples, initial = NA, style = 3)
    args_progress <- list(pb_burnin = pb_burnin,
                          pb_samples = pb_samples)
    
    # incporporate arguments unique to this K
    args_model$K <- K[i]
    
    # create argument list
    parallel_args[[i]] <- list(args_data = args_data,
                               args_model = args_model,
                               args_functions = args_functions,
                               args_progress = args_progress)
  }
  
  # ---------- run MCMC ----------
  
  # split into parallel and serial implementations
  if (!is.null(cluster)) { # run in parallel
    parallel::clusterEvalQ(cluster, library(silverblaze))
    output_raw <- parallel::clusterApplyLB(cl = cluster, parallel_args, run_mcmc_cpp)
  } else { # run in serial
    output_raw <- lapply(parallel_args, run_mcmc_cpp)
  }
  
  #return(output_raw)
  
  #------------------------
  
  # begin processing results
  if (!silent) {
    message("Processing results\n")
  }
  
  # loop through K
  ret <- list()
  all_converged <- TRUE
  for (i in seq_along(K)) {
    
    # create name lists
    group_names <- sprintf("group%s", seq_len(K[i]))
    rung_names <- sprintf("rung%s", seq_len(rungs))
    
    # ---------- raw mcmc results ----------
    
    # get iteration at which each rung converged
    convergence_iteration <- output_raw[[i]]$convergence_iteration
    
    # get loglikelihood in coda::mcmc format
    loglike_burnin <- coda::mcmc(t(rcpp_to_mat(output_raw[[i]]$loglike_burnin))[1:convergence_iteration,,drop = FALSE])
    loglike_sampling <- coda::mcmc(t(rcpp_to_mat(output_raw[[i]]$loglike_sampling)))
    colnames(loglike_sampling) <- colnames(loglike_burnin) <- rung_names
    
    # get source lon and lat in coda::mcmc format
    source_lon_burnin <- coda::mcmc(rcpp_to_mat(output_raw[[i]]$source_lon_burnin)[1:convergence_iteration,,drop = FALSE])
    source_lat_burnin <- coda::mcmc(rcpp_to_mat(output_raw[[i]]$source_lat_burnin)[1:convergence_iteration,,drop = FALSE])
    colnames(source_lon_burnin) <- colnames(source_lat_burnin) <- group_names
    
    source_lon_sampling <- coda::mcmc(rcpp_to_mat(output_raw[[i]]$source_lon_sampling))
    source_lat_sampling <- coda::mcmc(rcpp_to_mat(output_raw[[i]]$source_lat_sampling))
    colnames(source_lon_sampling) <- colnames(source_lat_sampling) <- group_names
    
    # get matrix of realised sources
    source_realised_burnin <- rcpp_to_mat(output_raw[[i]]$source_realised_burnin)[1:convergence_iteration,,drop = FALSE]
    source_realised_sampling <- rcpp_to_mat(output_raw[[i]]$source_realised_sampling)
    
    # get sigma in coda::mcmc format
    sigma_burnin <- coda::mcmc(rcpp_to_mat(output_raw[[i]]$sigma_burnin)[1:convergence_iteration,,drop = FALSE])
    sigma_sampling <- coda::mcmc(rcpp_to_mat(output_raw[[i]]$sigma_sampling))
    if (args_model$sigma_model == "single") {
      sigma_burnin <- sigma_burnin[, 1, drop = FALSE]
      sigma_sampling <- sigma_sampling[, 1, drop = FALSE]
      colnames(sigma_burnin) <- colnames(sigma_sampling) <- "all_groups"
    } else {
      colnames(sigma_burnin) <- colnames(sigma_sampling) <- group_names
    }
    
    # split method based on data type
    # get expected_popsize in coda::mcmc format
    expected_popsize_burnin <- coda::mcmc(rcpp_to_mat(output_raw[[i]]$ep_burnin)[1:convergence_iteration, ,drop = FALSE])
    expected_popsize_sampling <- coda::mcmc(rcpp_to_mat(output_raw[[i]]$ep_sampling))
    
    if (project$data$data_type == "counts" | project$data$data_type == "prevalence") {
      
      if (args_model$expected_popsize_model == "single") {
        expected_popsize_burnin <- expected_popsize_burnin[, 1, drop = FALSE]
        expected_popsize_sampling <- expected_popsize_sampling[, 1, drop = FALSE]
        colnames(expected_popsize_burnin) <- colnames(expected_popsize_sampling) <- "all_groups"
      } else {
        colnames(expected_popsize_burnin) <- colnames(expected_popsize_sampling) <- group_names
      }
      
      if (args_model$n_binom == TRUE) {  # negative binomial model
        
        # get alpha in coda::mcmc format
        alpha_burnin <- coda::mcmc(rcpp_to_mat(output_raw[[i]]$alpha_burnin)[1:convergence_iteration, ,drop = FALSE])
        alpha_sampling <- coda::mcmc(rcpp_to_mat(output_raw[[i]]$alpha_sampling))

      } else {
        alpha_burnin <- alpha_sampling <- NULL
      }
    
    } else if (project$data$data_type == "point-pattern") {
      alpha_burnin <- alpha_sampling <- NULL
    } else {
      stop("invalid data type in output")
    }
    
    # ---------- summary results ----------
    
    # get 95% credible intervals over sampling and burnin loglikelihoods
    loglike_intervals_burnin <- as.data.frame(t(apply(loglike_burnin, 2, quantile_95)))
    loglike_intervals_sampling <- as.data.frame(t(apply(loglike_sampling, 2, quantile_95)))
    
    # get 95% credible intervals over sigma
    sigma_intervals <- as.data.frame(t(apply(sigma_sampling, 2, quantile_95)))
    
    # get 95% credible intervals over expected_popsize
    expected_popsize_intervals <- as.data.frame(t(apply(expected_popsize_sampling, 2, quantile_95)))
    
    # split method based on data type
    if (project$data$data_type == "counts" | project$data$data_type == "prevalence") {
      
      # get 95% credible intervals over negative binomial parameters
      if (args_model$n_binom == TRUE) {
        alpha_intervals <- as.data.frame(t(quantile_95(alpha_sampling)))
      } else{
        alpha_intervals <- NULL
      }
    
    } else if (project$data$data_type == "point-pattern") {
      alpha_intervals <- NULL
    }
    
    # process Q-matrix
    # TODO - does anything need doing in the condition else if(project$data$data_type == "point-pattern") ?
    # TODO - should this really be class rgeoprofile_qmatrix? Are we anticipating making use of any defined special methods?
    qmatrix <- rcpp_to_mat(output_raw[[i]]$qmatrix)/samples
    if (project$data$data_type == "counts") {
      qmatrix[project$data$frame$counts == 0,] <- NA
    } else if (project$data$data_type == "prevalence") {
      qmatrix[project$data$frame$positive == 0,] <- NA
    }
    colnames(qmatrix) <- group_names
    class(qmatrix) <- "rgeoprofile_qmatrix"
    
    # get distribution of realised K
    realised_K <- tabulate(rowSums(source_realised_sampling), nbins = K[i])
    
    # create empty raster with correct properties
    raster_empty <- raster()
    raster::extent(raster_empty) <- raster::extent(spatial_prior)
    raster::res(raster_empty) <- raster::res(spatial_prior)
    
    # get breaks for kernel smoothing
    breaks_lon <- seq(raster::xmin(raster_empty), raster::xmax(raster_empty), raster::xres(raster_empty))
    breaks_lat <- seq(raster::ymin(raster_empty), raster::ymax(raster_empty), raster::yres(raster_empty))
    
    # produce posterior probability surface rasters
    prob_surface_split <- raster()
    prob_surface_mat <- 0
    for (k in seq_len(K[i])) {
      if (create_maps) {
        
        # get prob surface for this K by smoothing
        prob_surface_split_mat <- kernel_smooth(source_lon_sampling[,k],
                                                source_lat_sampling[,k],
                                                breaks_lon,
                                                breaks_lat)
        
        # flip and normalise
        prob_surface_split_mat <- prob_surface_split_mat[nrow(prob_surface_split_mat):1,]
        prob_surface_split_mat <- prob_surface_split_mat / sum(prob_surface_split_mat)
        
      } else {
        
        # store dummy surface
        prob_surface_split_mat <- matrix(NA, length(breaks_lat) - 1, length(breaks_lon) - 1)
      }
      
      # add as raster layer
      prob_surface_split_k <- setValues(raster_empty, prob_surface_split_mat)
      raster::values(prob_surface_split_k)[raster::values(spatial_prior) == 0] <- NA
      prob_surface_split <- raster::addLayer(prob_surface_split, prob_surface_split_k)
      
      # add to combined surface matrix
      prob_surface_mat <- prob_surface_mat + prob_surface_split_mat / K[i]
    }
    
    # make combined raster
    prob_surface <- setValues(raster_empty, prob_surface_mat)
    values(prob_surface)[raster::values(spatial_prior) == 0] <- NA
    
    # get probability surface over realised sources only
    prob_surface_realised_mat <- 0
    if (create_maps & K[i] > 1) { 
      
      # get prob surface by smoothing
      prob_surface_realised_mat <- kernel_smooth(source_lon_sampling[source_realised_sampling == TRUE],
                                                 source_lat_sampling[source_realised_sampling == TRUE],
                                                 breaks_lon,
                                                 breaks_lat)
      
      # flip and normalise
      prob_surface_realised_mat <- prob_surface_realised_mat[nrow(prob_surface_realised_mat):1,]
      prob_surface_realised_mat <- prob_surface_realised_mat / sum(prob_surface_realised_mat)
    } else {
      
      # store dummy surface
      prob_surface_realised_mat <- matrix(NA, length(breaks_lat) - 1, length(breaks_lon) - 1)
    }
    
    # make raster
    prob_surface_realised <- setValues(raster_empty, prob_surface_realised_mat)
    values(prob_surface_realised)[raster::values(spatial_prior) == 0] <- NA

    # produce geoprofile rasters
    geoprofile_split <- raster()
    geoprofile_mat <- 0
    for (k in seq_len(K[i])) {
      
      if (create_maps) {
        
        # make geoprofile matrix from probability surface
        geoprofile_split_mat <- rank(values(prob_surface_split[[k]]), ties.method = "first")
        geoprofile_split_mat <- 100 * (1 - geoprofile_split_mat/max(geoprofile_split_mat, na.rm = TRUE))
        
      } else {
        
        # store dummy surface
        geoprofile_split_mat <- matrix(NA, length(breaks_lat) - 1, length(breaks_lon) - 1)
      }
      
      # add as raster layer
      geoprofile_split_k <- setValues(raster_empty, geoprofile_split_mat)
      geoprofile_split <- addLayer(geoprofile_split, geoprofile_split_k)
    }
    
    # make combined raster
    geoprofile_mat <- rank(values(prob_surface), ties.method = "first", na.last = FALSE)
    geoprofile_mat <- 100 * (1 - geoprofile_mat / max(geoprofile_mat, na.rm = TRUE))
    geoprofile <- setValues(raster_empty, geoprofile_mat)
    values(geoprofile)[raster::values(spatial_prior) == 0] <- NA
    
    # get groprofile over realised sources only
    geoprofile_realised_mat <- rank(values(prob_surface_realised), ties.method = "first", na.last = FALSE)
    geoprofile_realised_mat <- 100 * (1 - geoprofile_realised_mat / max(geoprofile_realised_mat, na.rm = TRUE))
    geoprofile_realised <- setValues(raster_empty, geoprofile_realised_mat)
    values(geoprofile_realised)[raster::values(spatial_prior) == 0] <- NA
    
    # get whether rungs have converged
    converged <- output_raw[[i]]$rung_converged
    if (all_converged && any(!converged)) {
      all_converged <- FALSE
    }
    
    # ---------- ESS ----------
    
    # get ESS, unless using a fixed sigma model
    ESS <- apply(loglike_sampling, 2, function(x) {
      tc <- tryCatch(effectiveSize(x), error = function(e) e, warning = function(w) w)
      if (is(tc, "error")) {
        return(NA)
      } else {
        return(coda::effectiveSize(x))
      }
    })
    names(ESS) <- rung_names
    
    # ---------- model comparison statistics ----------
    
    # ---------- DIC ----------
    mu <- mean(loglike_sampling[,ncol(loglike_sampling)])
    sigma_sq <- var(loglike_sampling[,ncol(loglike_sampling)])
    DIC_gelman <- -2*mu + 4*sigma_sq
    
    # ---------- acceptance rates ----------
    # process acceptance rates
    source_accept_burnin <- matrix(unlist(output_raw[[i]]$source_accept_burnin), ncol = K[i], nrow = rungs, byrow = T)/convergence_iteration
    source_accept_sampling <- matrix(unlist(output_raw[[i]]$source_accept_sampling), ncol = K[i], nrow = rungs, byrow = T)/samples
    colnames(source_accept_burnin) <- colnames(source_accept_sampling) <- group_names
    rownames(source_accept_burnin) <- rownames(source_accept_sampling) <- rung_names
        
    sigma_accept_burnin <- matrix(unlist(output_raw[[i]]$sigma_accept_burnin), ncol = K[i], nrow = rungs, byrow = T)/convergence_iteration
    sigma_accept_sampling <- matrix(unlist(output_raw[[i]]$sigma_accept_sampling), ncol = K[i], nrow = rungs, byrow = T)/samples
    colnames(sigma_accept_burnin) <- colnames(sigma_accept_sampling) <- group_names
    rownames(sigma_accept_burnin) <- rownames(sigma_accept_sampling) <- rung_names
    
    # if prevelance or independent expected popsize model return acceptance rates
    if (project$data$data_type == "prevalence" | expected_popsize_model_numeric == 2) {
      
      expected_popsize_accept_burnin <- matrix(unlist(output_raw[[i]]$ep_accept_burnin), ncol = K[i], nrow = rungs, byrow = T)/convergence_iteration
      expected_popsize_accept_sampling <- matrix(unlist(output_raw[[i]]$ep_accept_sampling), ncol = K[i], nrow = rungs, byrow = T)/samples
      colnames(expected_popsize_accept_burnin) <- colnames(expected_popsize_accept_sampling) <- group_names
      rownames(expected_popsize_accept_burnin) <- rownames(expected_popsize_accept_sampling) <- rung_names
    
      if (args_model$n_binom == TRUE) {
         alpha_accept_burnin <- unlist(output_raw[[i]]$alpha_accept_burnin)/convergence_iteration
         alpha_accept_sampling <- unlist(output_raw[[i]]$alpha_accept_sampling)/samples
       } else {
         alpha_accept_burnin <- alpha_accept_sampling <- NULL
       }
      
    } else {
      alpha_accept_burnin <- alpha_accept_sampling <- NULL
      expected_popsize_accept_burnin <- expected_popsize_accept_sampling <- NULL
    } 
    # get Metropolis coupling acceptance rates
    coupling_accept_burnin <- output_raw[[i]]$coupling_accept_burnin/(convergence_iteration)
    coupling_accept_sampling <- output_raw[[i]]$coupling_accept_sampling/(samples)
      
    # ---------- save arguments ----------
    
    output_args <- list(burnin = burnin,
                        samples = samples,
                        auto_converge = auto_converge,
                        converge_test = converge_test,
                        pb_markdown = pb_markdown,
                        rungs = rungs,
                        silent = silent)
    
    # ---------- save results ----------
    
    # add to project
    project$output$single_set[[s]]$single_K[[K[i]]] <- list()
    
    project$output$single_set[[s]]$single_K[[K[i]]]$summary <- list(loglike_intervals_burnin = loglike_intervals_burnin,
                                                                    loglike_intervals_sampling = loglike_intervals_sampling,
                                                                    prob_surface_split = prob_surface_split,
                                                                    prob_surface = prob_surface,
                                                                    prob_surface_realised = prob_surface_realised,
                                                                    geoprofile_split = geoprofile_split,
                                                                    geoprofile = geoprofile,
                                                                    geoprofile_realised = geoprofile_realised,
                                                                    qmatrix = qmatrix,
                                                                    realised_K = realised_K,
                                                                    sigma_intervals = sigma_intervals, 
                                                                    expected_popsize_intervals = expected_popsize_intervals,
                                                                    alpha_intervals = alpha_intervals, 
                                                                    ESS = ESS,
                                                                    DIC_gelman = DIC_gelman,
                                                                    converged = converged,
                                                                    source_accept_burnin = source_accept_burnin,
                                                                    source_accept_sampling = source_accept_sampling,
                                                                    sigma_accept_burnin = sigma_accept_burnin,
                                                                    sigma_accept_sampling = sigma_accept_sampling,
                                                                    expected_popsize_accept_burnin = expected_popsize_accept_burnin,
                                                                    expected_popsize_accept_sampling = expected_popsize_accept_sampling,
                                                                    alpha_accept_burnin = alpha_accept_burnin,
                                                                    alpha_accept_sampling = alpha_accept_sampling,
                                                                    coupling_accept_burnin = coupling_accept_burnin,
                                                                    coupling_accept_sampling = coupling_accept_sampling,
                                                                    beta_vec = beta_vec)
                                                                    
    if (store_raw) {
      project$output$single_set[[s]]$single_K[[K[i]]]$raw <- list(loglike_burnin = loglike_burnin,
                                                                  source_lon_burnin = source_lon_burnin,
                                                                  source_lat_burnin = source_lat_burnin,
                                                                  source_realised_burnin = source_realised_burnin,
                                                                  sigma_burnin = sigma_burnin,
                                                                  expected_popsize_burnin = expected_popsize_burnin,
                                                                  alpha_burnin = alpha_burnin,
                                                                  loglike_sampling = loglike_sampling,
                                                                  source_lon_sampling = source_lon_sampling,
                                                                  source_lat_sampling = source_lat_sampling,
                                                                  source_realised_sampling = source_realised_sampling,
                                                                  sigma_sampling = sigma_sampling,
                                                                  expected_popsize_sampling = expected_popsize_sampling,
                                                                  alpha_sampling = alpha_sampling)
    
    }
    
    project$output$single_set[[s]]$single_K[[K[i]]]$function_call <- list(args = output_args,
                                                                          call = match.call())

  } # end loop over K

  # name output over K
  K_all <- length(project$output$single_set[[s]]$single_K)
  names(project$output$single_set[[s]]$single_K) <- paste0("K", 1:K_all)
  
  # ---------- tidy up and end ----------
  
  # reorder qmatrices
  project <- align_qmatrix(project)
  
  # get matrix of realised K over all model K (assuming K > 1)
  if (K_all > 1 & project$data$data_type == "prevalence"){
    realised_K_all <- mapply(function(x) {
      ret <- rep(NA, K_all)
      tmp <- x$summary$realised_K
      if (!is.null(tmp)) {
        ret[seq_along(tmp)] <- tmp
      }
      return(ret)
    }, project$output$single_set[[s]]$single_K)
    colnames(realised_K_all) <- sprintf("model_K%s", seq_len(K_all))
    rownames(realised_K_all) <- sprintf("realised_K%s", seq_len(K_all))
    project$output$single_set[[s]]$all_K$realised_K <- realised_K_all
  } else if (K_all == 1){
    project$output$single_set[[s]]$all_K$realised_K <- NULL
  }
  
  # run ring-search prior to MCMC
  if (sum(project$data$frame$counts) > 0 | sum(project$data$frame$positive) > 0) {
    ringsearch <- ring_search(project, spatial_prior)
    project$output$single_set[[s]]$all_K$ringsearch <- ringsearch
  }  
  
  # get DIC over all K
  DIC_gelman <- mapply(function(x) {
    ret <- x$summary$DIC_gelman
    if (is.null(ret)) {
      return(NA)
    } else {
      return(ret)
    }
  }, project$output$single_set[[s]]$single_K)
  DIC_gelman <- as.vector(unlist(DIC_gelman))
  project$output$single_set[[s]]$all_K$DIC_gelman <- data.frame(K = 1:length(DIC_gelman), DIC_gelman = DIC_gelman)
  
  # end timer
  if (!silent) {
    tdiff <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (tdiff < 60) {
      message(sprintf("Total run-time: %s seconds", round(tdiff, 2)))
    } else {
      message(sprintf("Total run-time: %s minutes", round(tdiff/60, 2)))
    }
  }
  
  # warning if any rungs in any MCMCs did not converge
  if (!all_converged && !silent) {
    message("\n**WARNING** at least one MCMC run did not converge\n")
  }
  
  # return invisibly
  invisible(project)
}

#------------------------------------------------
# align qmatrices over all K
#' @noRd
align_qmatrix <- function(project) {
  
  # get active set
  s <- project$active_set
  
  # extract objects of interest
  x <- project$output$single_set[[s]]$single_K
  
  # find values with output
  null_output <- mapply(function(y) {is.null(y$summary$qmatrix)}, x)
  w <- which(!null_output)
  
  # set template to first qmatrix
  template_qmatrix <- x[[w[1]]]$summary$qmatrix
  n <- nrow(template_qmatrix)
  c <- ncol(template_qmatrix)
  positive_sentinels <- which(!is.na(template_qmatrix[,1]))
  
  # loop through output
  best_perm <- NULL
  for (i in w) {
    
    # expand template
    qmatrix <- unclass(x[[i]]$summary$qmatrix)
    template_qmatrix <- cbind(template_qmatrix, matrix(0, n, i-c))
    
    # calculate cost matrix
    cost_mat <- matrix(0,i,i)
    for (k1 in 1:i) {
      for (k2 in 1:i) {
        cost_mat[k1,k2] <- sum(qmatrix[positive_sentinels,k1] * (log(qmatrix[positive_sentinels,k1]+1e-100) - log(template_qmatrix[positive_sentinels,k2]+1e-100)))
      }
    }
    
    # get lowest cost permutation
    best_perm <- call_hungarian(cost_mat)$best_matching + 1
    best_perm_order <- order(best_perm)
    
    # reorder qmatrix
    group_names <- paste0("group", 1:ncol(qmatrix))
    qmatrix <- qmatrix[, best_perm_order, drop = FALSE]
    colnames(qmatrix) <- group_names
    
    # reorder raw output
    if (!is.null(x[[i]]$raw)) {
      
      # reorder source_lon_burnin and source_lat_burnin
      source_lon_burnin <- x[[i]]$raw$source_lon_burnin[, best_perm_order, drop = FALSE]
      source_lat_burnin <- x[[i]]$raw$source_lat_burnin[, best_perm_order, drop = FALSE]
      colnames(source_lon_burnin) <- colnames(source_lat_burnin) <- group_names
      project$output$single_set[[s]]$single_K[[i]]$raw$source_lon_burnin <- source_lon_burnin
      project$output$single_set[[s]]$single_K[[i]]$raw$source_lat_burnin <- source_lat_burnin
      
      # reorder source_lon_sampling and source_lat_sampling
      source_lon_sampling <- x[[i]]$raw$source_lon_sampling[, best_perm_order, drop = FALSE]
      source_lat_sampling <- x[[i]]$raw$source_lat_sampling[, best_perm_order, drop = FALSE]
      colnames(source_lon_sampling) <- colnames(source_lat_sampling) <- group_names
      project$output$single_set[[s]]$single_K[[i]]$raw$source_lon_sampling <- source_lon_sampling
      project$output$single_set[[s]]$single_K[[i]]$raw$source_lat_sampling <- source_lat_sampling
      
      # reorder sigma
      sigma_burnin <- x[[i]]$raw$sigma_burnin
      sigma_sampling <- x[[i]]$raw$sigma_sampling
      if (ncol(sigma_sampling) > 1) {
        sigma_burnin <- sigma_burnin[, best_perm_order, drop = FALSE]
        sigma_sampling <- sigma_sampling[, best_perm_order, drop = FALSE]
        colnames(sigma_burnin) <- colnames(sigma_sampling) <- group_names
        project$output$single_set[[s]]$single_K[[i]]$raw$sigma_burnin <- sigma_burnin
        project$output$single_set[[s]]$single_K[[i]]$raw$sigma_sampling <- sigma_sampling
      }
      
      # reorder expected popsize
      expected_popsize_burnin <- x[[i]]$raw$expected_popsize_burnin
      expected_popsize_sampling <- x[[i]]$raw$expected_popsize_sampling
      if (ncol(expected_popsize_sampling) > 1) {
        expected_popsize_burnin <- expected_popsize_burnin[, best_perm_order, drop = FALSE]
        expected_popsize_sampling <- expected_popsize_sampling[, best_perm_order, drop = FALSE]
        colnames(expected_popsize_burnin) <- colnames(expected_popsize_sampling) <- group_names
        project$output$single_set[[s]]$single_K[[i]]$raw$expected_popsize_burnin <- expected_popsize_burnin
        project$output$single_set[[s]]$single_K[[i]]$raw$expected_popsize_sampling <- expected_popsize_sampling
      }
    }
    
    # reorder prob_surface_split
    prob_surface_split <- x[[i]]$summary$prob_surface_split
    layer_names <- names(prob_surface_split)
    prob_surface_split <- prob_surface_split[[best_perm_order]]
    names(prob_surface_split) <- layer_names
    project$output$single_set[[s]]$single_K[[i]]$summary$prob_surface_split <- prob_surface_split
    
    # reorder geoprofile_split
    geoprofile_split <- x[[i]]$summary$geoprofile_split
    layer_names <- names(geoprofile_split)
    geoprofile_split <- geoprofile_split[[best_perm_order]]
    names(geoprofile_split) <- layer_names
    project$output$single_set[[s]]$single_K[[i]]$summary$geoprofile_split <- geoprofile_split
    
    # reorder sigma_intervals
    sigma_intervals <- x[[i]]$summary$sigma_intervals[best_perm_order,,drop = FALSE]
    rownames(sigma_intervals) <- group_names
    project$output$single_set[[s]]$single_K[[i]]$summary$sigma_intervals <- sigma_intervals

    # reorder source_accept
    source_accept_sampling <- x[[i]]$summary$source_accept_sampling[, best_perm_order, drop = FALSE]
    colnames(source_accept_sampling) <- group_names
    project$output$single_set[[s]]$single_K[[i]]$summary$source_accept_sampling <- source_accept_sampling
    
    # reorder sigma_accept
    sigma_accept_sampling <- x[[i]]$summary$sigma_accept_sampling[, best_perm_order, drop = FALSE]
    colnames(sigma_accept_sampling) <- group_names
    project$output$single_set[[s]]$single_K[[i]]$summary$sigma_accept_sampling <- sigma_accept_sampling
    
    # reorder expected_popsize_accept
    if( project$data$data_type == "prevalence" | 
        project$parameter_sets[[s]]$n_binom == TRUE | 
        project$parameter_sets[[s]]$expected_popsize_model == "independent") {
      expected_popsize_accept_sampling <- x[[i]]$summary$expected_popsize_accept_sampling[, best_perm_order, drop = FALSE]
      colnames(expected_popsize_accept_sampling) <- group_names
      project$output$single_set[[s]]$single_K[[i]]$summary$expected_popsize_accept_sampling <- expected_popsize_accept_sampling
    }
    
    # qmatrix becomes template for next level up
    template_qmatrix <- qmatrix
    
    # store result
    class(qmatrix) <- "rgeoprofile_qmatrix"
    project$output$single_set[[s]]$single_K[[i]]$summary$qmatrix <- qmatrix
  }
  
  # return modified project
  return(project)
}

#------------------------------------------------
#' @title Optimise beta values for MCMC rungs
#'
#' @description repeatedly run burn in phase of MCMC to find vector of beta values 
#'              that achieve a required coupling acceptance rate between all rungs
#'
#' @param proj An rgeoprofile project
#' @param K The value or values of K to optimise beta for
#' @param target_acceptance The minimum acceptance probability to be reached between
#'                          MCMC chains before the process terminates 
#' @param max_iterations Run analyses a maximum number of times before terminating
#' @param beta_init An initial set of beta values to run - a sequence from zero
#'                  to one
#' @param silent print output from console
#' @param ... extra arguments to be sent to the run_mcmc function 
#'
#' @export

optimise_beta <- function(proj,
                          K = 3,
                          target_acceptance = 0.4,
                          max_iterations = 1e2,
                          beta_init = seq(0, 1, l = 10),
                          silent = FALSE,
                          ...) {
  
  # check inputs (only those not checked by later functions)
  assert_single_pos(target_acceptance)
  assert_bounded(target_acceptance)
  assert_single_pos_int(max_iterations, zero_allowed = FALSE)
  
  # get arguments from ellipsis
  args_list <- list(...)
  
  # check arguments not doubly defined
  if ("beta_manual" %in% names(args_list)) {
    stop("cannot define beta_manual in optimise_beta() function")
  }
  
  # initialise beta_vec
  beta_vec <- beta_init
  
  # iteratively improve beta_vec
  ret <- list()
  for (j in 1:max_iterations) {
    
    # run MCMC using beta_vec
    proj <- run_mcmc(project = proj,
                     K = K,
                     beta_manual = beta_vec,
                     silent = silent,
                     ...)
    
    # get current rungs
    rungs <- length(beta_vec)
    
    # get coupling rate of burnin phase
    coupling_burnin <- get_output(proj, "coupling_accept_burnin", K)
    
    # store results of this iteration
    ret$beta_vec <- c(ret$beta_vec, list(beta_vec))
    ret$coupling <- c(ret$coupling, list(coupling_burnin))
    
    # count how many probabilities are below target
    under_target <- sum(coupling_burnin < target_acceptance)
    
    # report progress to console
    if (!silent) {
      message("--------------------")
      message(sprintf("K = %s", K))
      message(sprintf("iteration = %s", j))
      message(sprintf("rungs = %s", rungs))
      message("beta values:")
      message(paste(signif(beta_vec, digits = 3), collapse = ", "))
      message("coupling acceptance:")
      message(paste(coupling_burnin, collapse = ", "))
      message("--------------------")
      message(paste0(under_target, " out of ", rungs - 1, 
                     " probabilities are under target (", target_acceptance, ")" )
)

    }
    
    # return if all coupling over target threshold
    if (all(coupling_burnin >= target_acceptance)) {
      return(ret)
    }
    
  
    # get index of those acceptance values less than target
    update <- which(coupling_burnin < target_acceptance)
    
    # loop over index of acceptance prob less than target
    for(i in 1:length(update)){
      
        # access index that needs changing
        update_single <- update[i]
        
        # calculate beta value in the middle of the two that produced the low transition probability 
        beta_increase <- 0.5*(beta_vec[update_single + 1] - beta_vec[update_single])
        new_beta <- beta_vec[update_single] + beta_increase
        
        # append this new value in beta_vec
        beta_vec <- append(x = beta_vec, values = new_beta, after = update_single)
        
        # update vector of indexes to account for new length of beta_vec
        update <- update + 1
    }
    
    
  }  # end loop over iterations
  
  # if reached this point then not converged within max_iterations
  warning(sprintf("optimise_beta() did not find solution within %s iterations", max_iterations))
  return(ret)

}

#------------------------------------------------
# ring-search
#' @importFrom raster extent<- extent res<- res setValues flip
#' @noRd
ring_search <- function(project, r) {
  
  # avoid "no visible binding" error
  counts <- positive <- NULL
  
  # check that there is at least one positive observation
  if (sum(project$data$frame$counts) == 0 & sum(project$data$frame$positive) == 0) {
    stop("ring search not possible: no positive counts")
  }
  
  # extract sentinel locations with at least one observation
  if(project$data$data_type == "counts"){
    data <- subset(project$data$frame, counts > 0)
  } else if(project$data$data_type == "prevalence"){
    data <- subset(project$data$frame, positive > 0)
  }
    
  sentinel_lon <- data$longitude
  sentinel_lat <- data$latitude

  # get breaks, midpoints, and coordinates of all points in grid
  breaks_lon <- seq(xmin(r), xmax(r), xres(r))
  breaks_lat <- seq(ymin(r), ymax(r), yres(r))
  midpoints_lon <- breaks_lon[-1] - xres(r)/2
  midpoints_lat <- breaks_lat[-1] - yres(r)/2
  coords <- expand.grid(midpoints_lon, midpoints_lat)

  # get distance between all cells and sentinels
  d <- apply(cbind(sentinel_lon, sentinel_lat), 1,
             function(x) lonlat_to_bearing(coords[,1], coords[,2], x[1], x[2])$gc_dist)

  # get minimum distance to any sentinel
  d_min <- apply(d, 1, min)

  # convert min distances to hitscore percentages
  hs <- rank(d_min)/length(d_min) * 100

  # get into raster format
  ret <- raster()
  raster::extent(ret) <- raster::extent(r)
  raster::res(ret) <- raster::res(r)
  ret <- raster::setValues(ret, hs)
  ret <- raster::flip(ret, 2)

  return(ret)
}

#------------------------------------------------
#' @title Calculate Gini coefficient
#'
#' @description Calculate Gini coefficient
#'
#' @param hs dataframe of hitscores
#'
#' @export

gini <- function(hs) {

  # check inputs
  assert_dataframe(hs)

  # drop lon/lat columns
  hs <- hs[ , !names(hs) %in% c("longitude", "latitude"), drop = FALSE]

  # get properties
  ns <- nrow(hs)
  hs_names <- colnames(hs)

  # get sorted values
  hs_sort <- apply(hs, 2, function(x) c(0, sort(x, na.last = TRUE)/100, 1))

  # get areas using trapezoidal rule
  y <- c(0:ns/ns, 1)
  hs_area <- apply(hs_sort, 2, function(x) sum(0.5*(y[-1]+y[-length(y)])*(x[-1]-x[-length(x)])) )

  # get Gini coefficient
  ret <- (hs_area - 0.5)/0.5

  # message if any NAs
  if (any(is.na(ret))) {
    message("Hitscores contain NA values (most likely due to sources outside the search area) leading to NA Gini coefficients")
  }

  return(ret)
}

#------------------------------------------------
#' @title Calculate realised sources
#'
#' @description Produce a plot that indicates the suitable number of realised 
#'              sources based upon sampling from the qmatrix.
#'
#' @param proj An rgeoprofile project
#' @param n_samples how many times we sample from the qmatrix
#' @param K the value of K to check. Note, there is no qmatrix for K = 1, all 
#'          points are allocated to a single source 
#'
#' @export

realised_sources <- function(proj, 
                             n_samples = 20, 
                             K = 2) {

  # check inputs
  assert_custom_class(proj, "rgeoprofile_project")
  assert_single_pos_int(n_samples, zero_allowed = FALSE)
  assert_greq(n_samples, 10, message = "at least 10 sampling iterations must be used")
  assert_single_pos_int(K, zero_allowed = FALSE)

  # get qmatrix from project output and clean
  qmat <- get_output(proj, "qmatrix", K = K, type = "summary")
  qmat <- qmat[complete.cases(qmat),,drop = FALSE]
  
  # loop over samples, for each sample, use the qmatrix to decide which source
  # this data point belongs to, then, make a note of the number of UNIQUE sources
  # for this sample
  ret <- mapply(function(i) {
    group_allocation <- apply(qmat, 1, function(x) sample(length(x), 1, prob = x))
    ret <- length(unique(group_allocation))
    return(ret)
  }, seq_len(n_samples))
  
  return(ret)
}
