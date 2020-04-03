
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
#'   Data must be formatted as a dataframe with samples in rows and loci in
#'   columns. If individuals are polyploid then multiple rows can be used per
#'   sample. Ploidy is allowed to vary between samples, and can be specified in
#'   multiple ways.
#'
#' @param project an \code{rgeoprofile_project}, as produced by the function
#'   \code{rgeoprofile_project()}
#' @param df a dataframe with columns that must conform to the following rules:
#'   \itemize{
#'     \item for \code{data_type = "counts"}, data must have columns
#'     "longitude", "latitude" and "counts".
#'     \item for \code{data_type = "prevalence"}, data must have columns
#'     "longitude", "latitude", "tested" and "positive"
#'     }
#' @param data_type the type of data, either "counts" or "prevalence".
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
  assert_in(data_type, c("counts", "prevalence"))
  if (data_type == "counts") {
    assert_in(c("longitude", "latitude", "counts"), names(df))
    assert_pos_int(df$counts, zero_allowed = TRUE)
  } else if (data_type == "prevalence") {
    assert_in(c("longitude", "latitude", "tested", "positive"), names(df))
    assert_pos_int(df$tested, zero_allowed = TRUE)
    assert_pos_int(df$positive, zero_allowed = TRUE)
    assert_leq(df$positive, df$tested)
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
#' @param guard_rail Extend each lon lat range by a proportion of the length of said range. E.g. a guard_rail of 0.05 increases the lon and lat range by 5 percent.
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

  # make raster grid
  dlon <- guard_rail*diff(range(range_lon))
  dlat <- guard_rail*diff(range(range_lat))
  lomn <- range_lon[1] - dlon/2
  lomx <- range_lon[2] + dlon/2
  lamn <- range_lat[1] - dlat/2
  lamx <- range_lat[2] + dlat/2
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
#' @param sentinel_radius the observation radius of sentinel sites.
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
#'
#' @export

new_set <- function(project,
                    spatial_prior = NULL,
                    sentinel_radius = 0.2,
                    sigma_model = "single",
                    sigma_prior_mean = 1,
                    sigma_prior_sd = 1,
                    expected_popsize_model = "single",
                    expected_popsize_prior_mean = 100,
                    expected_popsize_prior_sd = 10,
                    name = "(no name)") {
  
  # check inputs
  assert_custom_class(project, "rgeoprofile_project")
  if (!is.null(spatial_prior)) {
    assert_custom_class(spatial_prior, "RasterLayer")
  }
  assert_single_pos(sentinel_radius, zero_allowed = FALSE)
  assert_in(sigma_model, c("single", "independent"))
  assert_single_pos(sigma_prior_mean, zero_allowed = FALSE)
  assert_single_pos(sigma_prior_sd, zero_allowed = TRUE)
  assert_in(expected_popsize_model, c("single", "independent"))
  assert_single_pos(expected_popsize_prior_mean, zero_allowed = FALSE)
  assert_single_pos(expected_popsize_prior_sd, zero_allowed = TRUE)
  assert_single_string(name)
  
  # make spatial_prior from data limits if unspecified
  if (is.null(spatial_prior)) {
    range_lon <- range(project$data$frame$longitude)
    range_lon <- mean(range_lon) + 1.05*c(-1,1)*diff(range_lon)/2
    range_lat <- range(project$data$frame$latitude)
    range_lat <- mean(range_lat) + 1.05*c(-1,1)*diff(range_lat)/2
    spatial_prior <- raster_grid(range_lon, range_lat)
  }
  
  # get total study area in km^2
  study_area <- sum(raster::area(spatial_prior)[])
  
  # count current parameter sets and add one
  s <- length(project$parameter_sets) + 1
  
  # make new set active
  project$active_set <- s
  
  # create new parameter set
  project$parameter_sets[[s]] <- list(name = name,
                                      spatial_prior = spatial_prior,
                                      study_area = study_area,
                                      sentinel_radius = sentinel_radius,
                                      sigma_model = sigma_model,
                                      sigma_prior_mean = sigma_prior_mean,
                                      sigma_prior_sd = sigma_prior_sd,
                                      expected_popsize_model = expected_popsize_model,
                                      expected_popsize_prior_mean = expected_popsize_prior_mean,
                                      expected_popsize_prior_sd = expected_popsize_prior_sd)
  
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
#' @description Run the main RgeoProfile MCMC. Model parameters are taken from
#'   the current active parameter set, and MCMC parameters are passed in as
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
#' @param rungs the number of temperature rungs used.
#' @param auto_converge whether convergence should be assessed automatically
#'   every \code{converge_test} iterations, leading to termination of the
#'   burn-in phase. If \code{FALSE} then the full \code{burnin} iterations are
#'   used.
#' @param converge_test test for convergence every \code{convergence_test}
#'   iterations if \code{auto_converge} is being used.
#' @param coupling_on whether to implement Metropolis coupling.
#' @param GTI_pow power applied to thermodynamic rungs. Higher values lead to
#'   rungs clustered around zero.
#' @param beta_manual if defined, allows manual specification of thermodynamic
#'   powers used. Overrides \code{rungs} and \code{GTI_pow}.
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
                     rungs = 1,
                     auto_converge = TRUE,
                     converge_test = 1e2,
                     coupling_on = TRUE,
                     GTI_pow = 1,
                     beta_manual = NULL,
                     cluster = NULL,
                     pb_markdown = FALSE,
                     store_raw = TRUE,
                     create_maps = TRUE,
                     silent = !is.null(cluster)) {
  
  # start timer
  t0 <- Sys.time()
  
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
  assert_single_pos(GTI_pow)
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
  
  # define beta_vec manually or from rungs and GTI_pow
  beta_vec <- beta_manual
  if (is.null(beta_vec)) {
    if (rungs == 1) {
      beta_vec <- 1
    } else {
      beta_vec <- seq(0, 1, l = rungs)^GTI_pow
    }
  }
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
  }
  
  # input arguments list
  args_inputs <- list(burnin = burnin,
                      samples = samples,
                      beta_vec = beta_vec,
                      auto_converge = auto_converge,
                      converge_test = converge_test,
                      coupling_on = coupling_on,
                      pb_markdown = pb_markdown,
                      silent = silent)
  
  # extract spatial prior object
  spatial_prior <- project$parameter_sets[[s]]$spatial_prior
  spatial_prior_values <- raster::values(spatial_prior)
  spatial_prior_values[is.na(spatial_prior_values)] <- 0
  
  # initialise sources in a non-NA cell
  source_init <- raster::xyFromCell(spatial_prior, which(!is.na(raster::values(spatial_prior)))[1])
  
  # convert sigma_model and expected_popsize_model to numeric
  sigma_model_numeric <- match(project$parameter_sets[[s]]$sigma_model, c("single", "independent"))
  fixed_sigma_model <- project$parameter_sets[[s]]$sigma_prior_sd == 0
  
  expected_popsize_model_numeric <- match(project$parameter_sets[[s]]$expected_popsize_model, c("single", "independent"))
  fixed_expected_popsize_model <- project$parameter_sets[[s]]$expected_popsize_prior_sd == 0
  
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
                          source_init = source_init,
                          sigma_model_numeric = sigma_model_numeric,
                          fixed_sigma_model = fixed_sigma_model,
                          expected_popsize_model_numeric = expected_popsize_model_numeric,
                          fixed_expected_popsize_model = fixed_expected_popsize_model)
  
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
    cat("Processing results\n")
  }
  
  # loop through K
  ret <- list()
  all_converged <- TRUE
  for (i in 1:length(K)) {
    
    # create name lists
    group_names <- paste0("group", 1:K[i])
    rung_names <- paste0("rung", 1:rungs)
    
    # ---------- raw mcmc results ----------
    
    # get iteration at which each rung converged
    convergence_iteration <- output_raw[[i]]$convergence_iteration
    
    # get loglikelihood in coda::mcmc format
    loglike_burnin <- coda::mcmc(t(rcpp_to_mat(output_raw[[i]]$loglike_burnin))[1:convergence_iteration,,drop = FALSE])
    loglike_sampling <- coda::mcmc(t(rcpp_to_mat(output_raw[[i]]$loglike_sampling)))
    colnames(loglike_sampling) <- colnames(loglike_burnin) <- rung_names
    
    # get source lon lat in coda::mcmc format
    source_lon_burnin <- coda::mcmc(rcpp_to_mat(output_raw[[i]]$source_lon_burnin)[1:convergence_iteration,,drop = FALSE])
    source_lat_burnin <- coda::mcmc(rcpp_to_mat(output_raw[[i]]$source_lat_burnin)[1:convergence_iteration,,drop = FALSE])
    colnames(source_lon_burnin) <- colnames(source_lat_burnin) <- group_names
    source_lon_sampling <- coda::mcmc(rcpp_to_mat(output_raw[[i]]$source_lon_sampling))
    source_lat_sampling <- coda::mcmc(rcpp_to_mat(output_raw[[i]]$source_lat_sampling))
    colnames(source_lon_sampling) <- colnames(source_lat_sampling) <- group_names
    
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
    # get expected_popsize in coda::mcmc format
    expected_popsize_burnin <- coda::mcmc(rcpp_to_mat(output_raw[[i]]$ep_burnin)[1:convergence_iteration, ,drop = FALSE])
    expected_popsize_sampling <- coda::mcmc(rcpp_to_mat(output_raw[[i]]$ep_sampling))
    if (args_model$expected_popsize_model == "single") {
      expected_popsize_burnin <- expected_popsize_burnin[, 1, drop = FALSE]
      expected_popsize_sampling <- expected_popsize_sampling[, 1, drop = FALSE]
      colnames(expected_popsize_burnin) <- colnames(expected_popsize_sampling) <- "all_groups"
    } else {
      colnames(expected_popsize_burnin) <- colnames(expected_popsize_sampling) <- group_names
    }
    
    # ---------- summary results ----------
    # get 95% credible intervals over sampling and burnin loglikelihoods
    loglike_intervals_burnin <- as.data.frame(t(apply(loglike_burnin, 2, quantile_95)))
    loglike_intervals_sampling <- as.data.frame(t(apply(loglike_sampling, 2, quantile_95)))
    
    # get 95% credible intervals over sigma
    sigma_intervals <- as.data.frame(t(apply(sigma_sampling, 2, quantile_95)))
    
    # get 95% credible intervals over expected_popsize
    expected_popsize_intervals <- as.data.frame(t(apply(expected_popsize_sampling, 2, quantile_95)))
    
    # process Q-matrix
    qmatrix <- rcpp_to_mat(output_raw[[i]]$qmatrix)/samples
    if(project$data$data_type == "counts"){
      qmatrix[project$data$frame$counts == 0,] <- rep(NA, K[i])
    } else if(project$data$data_type == "prevalence"){
      qmatrix[project$data$frame$positive == 0,] <- rep(NA, K[i])
    }    
    colnames(qmatrix) <- group_names
    class(qmatrix) <- "rgeoprofile_qmatrix"
    
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
    for (k in 1:K[i]) {
      
      if (create_maps) {
        # get prob_surface for this K by smoothing
        prob_surface_split_mat <- kernel_smooth(source_lon_sampling[,k],
                                                source_lat_sampling[,k],
                                                breaks_lon,
                                                breaks_lat)
        prob_surface_split_mat <- prob_surface_split_mat[nrow(prob_surface_split_mat):1,]
        prob_surface_split_mat <- prob_surface_split_mat/sum(prob_surface_split_mat)
      } else {
        prob_surface_split_mat <- matrix(NA, length(breaks_lon) - 1, length(breaks_lat) - 1)
      }
      
      # add raster layer
      prob_surface_split_k <- setValues(raster_empty, prob_surface_split_mat)
      raster::values(prob_surface_split_k)[is.na(raster::values(spatial_prior))] <- NA
      prob_surface_split <- raster::addLayer(prob_surface_split, prob_surface_split_k)
      
      # add to combined surface matrix
      prob_surface_mat <- prob_surface_mat + prob_surface_split_mat/K[i]
    }
    
    # make combined raster
    prob_surface <- setValues(raster_empty, prob_surface_mat)
    values(prob_surface)[is.na(values(spatial_prior))] <- NA
    
    # produce geoprofile rasters
    geoprofile_split <- raster()
    geoprofile_mat <- 0
    for (k in 1:K[i]) {
      
      if (create_maps) {
        # make geoprofile matrix from probability surface
        geoprofile_split_mat <- rank(values(prob_surface_split[[k]]), ties.method = "first")
        geoprofile_split_mat <- 100 * (1 - geoprofile_split_mat/max(geoprofile_split_mat, na.rm = TRUE))
      } else {
        geoprofile_split_mat <- matrix(NA, length(breaks_lon) - 1, length(breaks_lat) - 1)
      }
      
      # add raster layer
      geoprofile_split_k <- setValues(raster_empty, geoprofile_split_mat)
      geoprofile_split <- addLayer(geoprofile_split, geoprofile_split_k)
    }
    
    # make combined raster
    geoprofile_mat <- rank(values(prob_surface), ties.method = "first", na.last = FALSE)
    geoprofile_mat <- 100 * (1 - geoprofile_mat/max(geoprofile_mat, na.rm = TRUE))
    geoprofile <- setValues(raster_empty, geoprofile_mat)
    values(geoprofile)[is.na(values(spatial_prior))] <- NA
    
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
        return(effectiveSize(x))
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
    source_accept_burnin <- output_raw[[i]]$source_accept_burnin/burnin
    source_accept_sampling <- output_raw[[i]]$source_accept_sampling/samples
    names(source_accept_burnin) <- names(source_accept_sampling) <- group_names
    
    sigma_accept_burnin <- output_raw[[i]]$sigma_accept_burnin/burnin
    sigma_accept_sampling <- output_raw[[i]]$sigma_accept_sampling/samples
    names(sigma_accept_burnin) <- names(sigma_accept_sampling) <- group_names
    
    expected_popsize_accept_burnin <- expected_popsize_accept_sampling <- NULL
    
    if (project$data$data_type == "prevalence" | expected_popsize_model_numeric == 2) {
      expected_popsize_accept_burnin <- output_raw[[i]]$ep_accept_burnin/burnin
      expected_popsize_accept_sampling <- output_raw[[i]]$ep_accept_sampling/samples
    }
    
    coupling_accept_burnin <- output_raw[[i]]$coupling_accept_burnin/(convergence_iteration)
    coupling_accept_sampling <- output_raw[[i]]$coupling_accept_sampling/(samples)
      
    # ---------- save arguments ----------
    
    output_args <- list(burnin = burnin,
                        samples = samples,
                        auto_converge = auto_converge,
                        converge_test = converge_test,
                        pb_markdown = pb_markdown,
                        GTI_pow,
                        rungs = rungs,
                        silent = silent)
    
    # ---------- save results ----------
    
    # add to project
    project$output$single_set[[s]]$single_K[[K[i]]] <- list()
    
    project$output$single_set[[s]]$single_K[[K[i]]]$summary <- list(loglike_intervals_burnin = loglike_intervals_burnin,
                                                                    loglike_intervals_sampling = loglike_intervals_sampling,
                                                                    prob_surface_split = prob_surface_split,
                                                                    prob_surface = prob_surface,
                                                                    geoprofile_split = geoprofile_split,
                                                                    geoprofile = geoprofile,
                                                                    qmatrix = qmatrix,
                                                                    sigma_intervals = sigma_intervals, 
                                                                    expected_popsize_intervals = expected_popsize_intervals,
                                                                    ESS = ESS,
                                                                    DIC_gelman = DIC_gelman,
                                                                    converged = converged,
                                                                    source_accept_burnin = source_accept_burnin,
                                                                    source_accept_sampling = source_accept_sampling,
                                                                    sigma_accept_burnin = sigma_accept_burnin,
                                                                    sigma_accept_sampling = sigma_accept_sampling,
                                                                    expected_popsize_accept_burnin = expected_popsize_accept_burnin,
                                                                    expected_popsize_accept_sampling = expected_popsize_accept_sampling,
                                                                    coupling_accept_burnin = coupling_accept_burnin,
                                                                    coupling_accept_sampling = coupling_accept_sampling,
                                                                    beta_vec = beta_vec,
                                                                    GTI_pow = GTI_pow)
    
    if (store_raw) {
      project$output$single_set[[s]]$single_K[[K[i]]]$raw <- list(loglike_burnin = loglike_burnin,
                                                                  source_lon_burnin = source_lon_burnin,
                                                                  source_lat_burnin = source_lat_burnin,
                                                                  sigma_burnin = sigma_burnin,
                                                                  expected_popsize_burnin = expected_popsize_burnin,
                                                                  loglike_sampling = loglike_sampling,
                                                                  source_lon_sampling = source_lon_sampling,
                                                                  source_lat_sampling = source_lat_sampling,
                                                                  sigma_sampling = sigma_sampling,
                                                                  expected_popsize_sampling = expected_popsize_sampling)
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
      names(source_lon_burnin) <- names(source_lat_burnin) <- group_names
      project$output$single_set[[s]]$single_K[[i]]$raw$source_lon_burnin <- source_lon_burnin
      project$output$single_set[[s]]$single_K[[i]]$raw$source_lat_burnin <- source_lat_burnin
      
      # reorder source_lon_sampling and source_lat_sampling
      source_lon_sampling <- x[[i]]$raw$source_lon_sampling[, best_perm_order, drop = FALSE]
      source_lat_sampling <- x[[i]]$raw$source_lat_sampling[, best_perm_order, drop = FALSE]
      names(source_lon_sampling) <- names(source_lat_sampling) <- group_names
      project$output$single_set[[s]]$single_K[[i]]$raw$source_lon_sampling <- source_lon_sampling
      project$output$single_set[[s]]$single_K[[i]]$raw$source_lat_sampling <- source_lat_sampling
      
      # reorder sigma
      sigma_burnin <- x[[i]]$raw$sigma_burnin
      sigma_sampling <- x[[i]]$raw$sigma_sampling
      if (ncol(sigma_sampling) > 1) {
        sigma_burnin <- sigma_burnin[, best_perm_order, drop = FALSE]
        sigma_sampling <- sigma_sampling[, best_perm_order, drop = FALSE]
        names(sigma_burnin) <- names(sigma_sampling) <- group_names
        project$output$single_set[[s]]$single_K[[i]]$raw$sigma_burnin <- sigma_burnin
        project$output$single_set[[s]]$single_K[[i]]$raw$sigma_sampling <- sigma_sampling
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
    source_accept_sampling <- x[[i]]$summary$source_accept_sampling[best_perm_order]
    names(source_accept_sampling) <- group_names
    project$output$single_set[[s]]$single_K[[i]]$summary$source_accept_sampling <- source_accept_sampling
    
    # reorder sigma_accept
    sigma_accept_sampling <- x[[i]]$summary$sigma_accept_sampling[best_perm_order]
    names(sigma_accept_sampling) <- group_names
    project$output$single_set[[s]]$single_K[[i]]$summary$sigma_accept_sampling <- sigma_accept_sampling
    
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
# ring-search
#' @importFrom raster extent<- extent res<- res setValues flip
#' @noRd
ring_search <- function(project, r) {
  
  # check that there is at least one positive observation
  if (sum(project$data$frame$counts) == 0 & sum(project$data$frame$positive) == 0) {
    stop("ring search not possible: no positive counts")
  }
  
  # avoid "no visible binding" error
  counts <- NULL
  
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
# repeatedly run analysis to find vector of beta values that achieves a required
# coupling acceptance rate between all rungs
#' @noRd
optimise_beta <- function(proj,
                          K = 3,
                          target_acceptance = 0.4,
                          max_iterations = 1e3,
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

  # TODO Keep or remove Bob's version of the beta function 
  # # update beta sequence
  # for (i in 1:(rungs - 2)) {
  #   if (coupling_burnin[i] < target_acceptance) {
  # 
  #     # get acceptance rate relative to target (small additions to control max
  #     # and min possible values)
  #     rel_accept <- (coupling_burnin[i] + 0.01)/(target_acceptance + 0.1)
  # 
  #     # calculate how far to adjust sequence based on rel_accept
  #     move_left <- diff(beta_vec)[i] * (1 - rel_accept)
  # 
  #     # adjust sequence
  #     beta_vec[-(1:i)] <- beta_vec[-(1:i)] - move_left
  #   }
  # }
  # 
  # # if final coupling value is less than target, add another rung
  # if (coupling_burnin[rungs - 1] < target_acceptance) {
  #   beta_vec <- c(beta_vec, 1)
  #   rungs <- length(beta_vec)
  # 
  # 
  # # ensure final value in sequence always equals 1
  # beta_vec[length(beta_vec)] <- 1

}
