
#------------------------------------------------
#' @title Create a new RgeoProfile project
#' 
#' @description Create a new RgeoProfile project.
#' 
#' @export

rgeoprofile_project <- function() {
  
  # initialise project with default values
  ret <- list(data = NULL,
              parameter_sets = NULL,
              active_set = 0,
              output = list(single_set = list())
              )
  
  # create class and return
  class(ret) <- "rgeoprofile_project"
  return(ret)
}

#------------------------------------------------
# overload print() function for rgeoprofile_project
#' @method print rgeoprofile_project
#' @export
print.rgeoprofile_project <- function(x, ...) {
  
  # print summary
  summary(x)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
#' @title Ordinary print function for class rgeoprofile_project
#' 
#' @description Calling \code{print()} on an object of class rgeoprofile_project results
#'   in custom output. This function therefore stands in for the base
#'   \code{print()} function, and is equivalent to running
#'   \code{print(unclass(x))}.
#' 
#' @param x object of class \code{rgeoprofile_project}
#' @param ... other arguments passed to \code{print()}
#' 
#' @export

print_full <- function(x, ...) {
  
  # check inputs
  assert_custom_class(x, "rgeoprofile_project")
  
  # print un-classed object
  print(unclass(x), ...)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
# overload summary() function for rgeoprofile_project
#' @method summary rgeoprofile_project
#' @export
summary.rgeoprofile_project <- function(object, ...) {
  
  # print data summary
  cat("DATA:\n")
  if (is.null(object$data)) {
    cat("   (none loaded)\n")
  } else if(object$data$data_type == "counts") {
    
    # extract data properties
    n_sentinel <- nrow(object$data$frame)
    n_positive <- sum(object$data$frame$counts > 0)
    n_obs <- sum(object$data$frame$counts)
    
    cat("   type: count\n")
    cat(sprintf("   sentinel sites = %s (%s positive, %s negative)\n", n_sentinel, n_positive, n_sentinel - n_positive))
    cat(sprintf("   total observations = %s\n", n_obs))
  } else if(object$data$data_type == "prevalence") {
  
    # extract data properties
    n_trial <- nrow(object$data$frame)
    n_tested <- sum(object$data$frame$tested)
    n_positive <- sum(object$data$frame$positive)
    
    cat("   type: prevalence\n")
    cat(sprintf("   trial sites = %s\n", n_trial))
    cat(sprintf("   total tested = %s\n", n_tested))
    cat(sprintf("   total testing positive = %s\n", n_positive))
    
  } else if(object$data$data_type == "point-pattern") {
    
    # extract data properties
    n_obs <- nrow(object$data$frame)
    
    cat("   type: point-pattern\n")
    cat(sprintf("   total observations = %s\n", n_obs))
  }
  cat("\n")
  
  # print parameter sets summary
  cat("PARAMETER SETS:\n")
  if (length(object$parameter_sets) == 0) {
    cat("   (none defined)\n")
  } else {
    
    # print names of all sets
    s <- object$active_set
    for (i in 1:length(object$parameter_sets)) {
      
      # star next to active set
      if (i == s) {
        cat(" * ")
      } else {
        cat("   ")
      }
      
      # print name of set
      cat(sprintf("SET%s: %s\n", i, object$parameter_sets[[i]]$name))
    }
    cat("\n")
    
    # get details of active set
    sentinel_radius <- expected_popsize_prior_mean <- expected_popsize_prior_sd <- NA

    if(object$data$data_type == "counts"){
      sentinel_radius <- object$parameter_sets[[s]]$sentinel_radius
    } 
        
    if(object$data$data_type == "counts" | object$data$data_type == "prevalence"){
      expected_popsize_model <- object$parameter_sets[[s]]$expected_popsize_model
      expected_popsize_prior_mean <- object$parameter_sets[[s]]$expected_popsize_prior_mean
      expected_popsize_prior_sd <- object$parameter_sets[[s]]$expected_popsize_prior_sd
    }
    
    dispersal_model <- object$parameter_sets[[s]]$dispersal_model
    spatial_prior <- object$parameter_sets[[s]]$spatial_prior
    study_area <- round(object$parameter_sets[[s]]$study_area, digits = 3)
    min_lon <- round(xmin(spatial_prior), digits = 3)
    max_lon <- round(xmax(spatial_prior), digits = 3)
    min_lat <- round(ymin(spatial_prior), digits = 3)
    max_lat <- round(ymax(spatial_prior), digits = 3)
    cells_lon <- ncol(spatial_prior)
    cells_lat <- nrow(spatial_prior)
    sigma_model <- object$parameter_sets[[s]]$sigma_model
    sigma_prior_mean <- object$parameter_sets[[s]]$sigma_prior_mean
    sigma_prior_sd <- object$parameter_sets[[s]]$sigma_prior_sd
    
    # print details of active set
    cat(sprintf("ACTIVE SET: SET%s\n", s))
    
    if(object$data$data_type == "counts"){
      cat(sprintf("   sentinel radius = %s\n", sentinel_radius))
    }    
    
    cat(sprintf("   dispersal/decay kernel = %s\n", dispersal_model))
    cat(sprintf("   spatial prior:\n"))
    cat(sprintf("      longitude range = [%s, %s]\n", min_lon, max_lon))
    cat(sprintf("      latitude range = [%s, %s]\n", min_lat, max_lat))
    cat(sprintf("      area (km^2) = %s\n", study_area))
    cat(sprintf("      cells = %s, %s (lon, lat)\n", cells_lon, cells_lat))
    cat(sprintf("   sigma prior:\n"))
    cat(sprintf("      model = %s\n", sigma_model))
    
    if (sigma_prior_sd == 0) {
      cat(sprintf("      value = %s (exact prior)\n", sigma_prior_mean))
    } else {
      cat(sprintf("      prior mean = %s\n", sigma_prior_mean))
      cat(sprintf("      prior SD = %s\n", sigma_prior_sd))
    }
    
    if(object$data$data_type == "counts" | object$data$data_type == "prevalence"){
      cat(sprintf("   expected population size prior:\n"))
      cat(sprintf("      model = %s\n", expected_popsize_model))
    
      if (expected_popsize_prior_sd == 0) {
          cat(sprintf("      value = %s (exact prior)\n", expected_popsize_prior_mean))
      } else {
          cat(sprintf("      prior mean = %s\n", expected_popsize_prior_mean))
          cat(sprintf("      prior SD = %s\n", expected_popsize_prior_sd))
      }
    } 
     
    
  }
  cat("\n")
  
}

#------------------------------------------------
#' @title Determine if object is of class rgeoprofile_project
#'
#' @description Determine if object is of class rgeoprofile_project
#'
#' @param x object of class \code{rgeoprofile_project}
#'
#' @export

is.rgeoprofile_project <- function(x) {
  inherits(x, "rgeoprofile_project")
}
