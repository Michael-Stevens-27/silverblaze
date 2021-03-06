% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{overlay_sentinels}
\alias{overlay_sentinels}
\title{Add sentinel sites to dynamic map}
\usage{
overlay_sentinels(
  myplot,
  project,
  sentinel_radius = NULL,
  fill = TRUE,
  fill_colour = c(grey(0.5), "red"),
  fill_opacity = 0.5,
  border = FALSE,
  border_colour = "black",
  border_weight = 1,
  border_opacity = 1,
  legend = FALSE,
  label = FALSE,
  label_size = 15
)
}
\arguments{
\item{myplot}{dynamic map produced by \code{plot_map()} function}

\item{project}{an RgeoProfile project, as produced by the function
\code{rgeoprofile_project()}.}

\item{sentinel_radius}{the radius of sentinel sites. Taken from the active
parameter set if unspecified.}

\item{fill}{whether to fill circles.}

\item{fill_colour}{colour of circle fill.}

\item{fill_opacity}{fill opacity.}

\item{border}{whether to add border to circles.}

\item{border_colour}{colour of circle borders.}

\item{border_weight}{thickness of circle borders.}

\item{border_opacity}{opacity of circle borders.}

\item{legend}{whether to add a legend for site count.}

\item{label}{whether to label sentinel sites with densities.}

\item{label_size}{size of the label.}
}
\description{
Add sentinel sites to dynamic map
}
\examples{
\dontshow{p <- rgeoprofile_file("tutorial1_project.rds")}
plot1 <- plot_map()
plot1 <- overlay_sentinels(plot1, project = p, fill_opacity = 0.9, fill = TRUE,
                           fill_colour = c(grey(0.7), "red"), border = c(FALSE, TRUE),
                           border_colour = "black",border_weight = 0.5)
plot1
}
