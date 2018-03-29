
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
