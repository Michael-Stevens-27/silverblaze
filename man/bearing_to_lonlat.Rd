% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{bearing_to_lonlat}
\alias{bearing_to_lonlat}
\title{Get spatial coordinate given an origin, a great circle distance and a
  bearing}
\usage{
bearing_to_lonlat(origin_lon, origin_lat, bearing, gc_dist)
}
\arguments{
\item{origin_lon}{The origin longitude}

\item{origin_lat}{The origin latitude}

\item{bearing}{The angle in degrees relative to due north}

\item{gc_dist}{The great circle distance in (km)}
}
\description{
Calculate destination lat/lon given an origin, a great circle
  distance of travel, and a bearing.
}
\examples{
# one degree longitude is approximately 111km at the equator. Therefore if we
# travel 111km due east from the coordinate {0,0} we can verify that we have
# moved approximately 1 degree longitude and zero degrees latitude
bearing_to_lonlat(0, 0, 90, 111)
}
