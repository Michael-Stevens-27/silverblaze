% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{dispersal_sphere}
\alias{dispersal_sphere}
\title{Draw from spherical distribution}
\usage{
dispersal_sphere(
  n,
  centre_lon,
  centre_lat,
  dispersal_model = "normal",
  scale = 1
)
}
\arguments{
\item{n}{The number of points to draw}

\item{centre_lon}{The mean longitude of the distribution}

\item{centre_lat}{The mean latitude of the distribution}

\item{dispersal_model}{The model we draw points from (normal or cauchy)}

\item{scale}{The scale parameter of the dispersal distribution defined by 
the parameter "dispersal_model"}
}
\description{
Draw from distribution converted to spherical coordinate
  system. Points are first drawn from an ordinary cartesian 2D normal
  distribution. The distances to points are then assumed to be great circle
  distances, and are combined with a random bearing from the point
  {centre_lat, centre_lon} to produce a final set of lat/lon points.
}
\examples{
dispersal_sphere(n = 100, centre_lat = 0, centre_lon = 0, dispersal_model = "normal", scale = 1)
}
