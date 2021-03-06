% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_structure}
\alias{plot_structure}
\title{Posterior allocation plot}
\usage{
plot_structure(project, K = NULL, divide_ind_on = FALSE)
}
\arguments{
\item{project}{an RgeoProfile project, as produced by the function
\code{rgeoprofile_project()}.}

\item{K}{which value of K to produce the plot for.}

\item{divide_ind_on}{whether to add dividing lines between bars.}
}
\description{
Produce posterior allocation plot of current active set.
}
\examples{
\dontshow{p <- rgeoprofile_file("tutorial1_project.rds")}
# Plot the structure for a single K value.
plot_structure(project = p, K = 2)
# Similarly, plot the allocation structure for every K.
plot_structure(project = p, divide_ind_on = TRUE)
}
