#------------------------------------------------
#' @title Bayesian Geographic Profiling from Presence-Absence Data
#'
#' @description Bayesian geographic profiling model that allows estimation of
#'   source locations given presence and absence data.
#'
#' @docType package
#' @name silverblaze
NULL

#------------------------------------------------
#' @useDynLib silverblaze, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#------------------------------------------------
# unload DLL when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("silverblaze", libpath)  # nocov
}