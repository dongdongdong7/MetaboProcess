% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/A_parameter_optimization.R
\name{edgeTrack}
\alias{edgeTrack}
\title{edgeTrack}
\usage{
edgeTrack(chrDf)
}
\arguments{
\item{chrDf}{A chrDf with baseline.}
}
\value{
A chrDf but start and end are edges.
}
\description{
locate the edge of a peak.
}
\examples{
chrDf_j <- ZOIList[[j]]
chrDf_j <- baselineEs(chrDf_j, loops = 8, tol_m = 30)
chrDf_j_new <- edgeTrack(chrDf = chrDf_j)
}
