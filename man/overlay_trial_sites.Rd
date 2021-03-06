% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{overlay_trial_sites}
\alias{overlay_trial_sites}
\title{Add trial sites to dynamic map}
\usage{
overlay_trial_sites(
  myplot,
  project,
  fill = TRUE,
  fill_colour = c(grey(0.5), "red"),
  fill_opacity = 0.5,
  border = FALSE,
  border_colour = "black",
  border_weight = 1,
  border_opacity = 1,
  legend = FALSE,
  site_radius = 20,
  plot_type = "piecharts"
)
}
\arguments{
\item{myplot}{dynamic map produced by \code{plot_map()} function}

\item{project}{an RgeoProfile project, as produced by the function
\code{rgeoprofile_project()}.}

\item{fill}{whether to fill circles.}

\item{fill_colour}{colour of circle fill.}

\item{fill_opacity}{fill opacity.}

\item{border}{whether to add border to circles.}

\item{border_colour}{colour of circle borders.}

\item{border_weight}{thickness of circle borders.}

\item{border_opacity}{opacity of circle borders.}

\item{legend}{whether to add a legend for site count.}

\item{site_radius}{radius in Km shown at each site}

\item{plot_type}{plot trial sites as circles or piecharts}
}
\description{
Add trial sites to dynamic map
}
\examples{
#\dontshow{p <- rgeoprofile_file("tutorial2_project.rds")}
#plot1 <- plot_map()
#plot1 <- overlay_trial_sites(plot1, project = p, fill_opacity = 0.9, fill = TRUE,
#                           fill_colour = c(grey(0.7), "red"), border = c(FALSE, TRUE),
#                           border_colour = "black",border_weight = 0.5)
#plot1
}
