% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/old_function.R
\name{getChromPeak}
\alias{getChromPeak}
\title{getChromPeak}
\usage{
getChromPeak(
  data,
  cpid,
  expandRt = 0,
  smooth = FALSE,
  method = "mean",
  size = 3,
  p = 3,
  n = p + 3 - p\%\%2,
  m = 0,
  ts = 1
)
}
\arguments{
\item{data}{A XcmsExperiment object.}

\item{cpid}{cpid.}

\item{expandRt}{expandRt.}

\item{smooth}{Whether or not smooth.}

\item{method}{Smooth method.}

\item{size}{size.}

\item{p}{p.}

\item{n}{n.}

\item{m}{m.}

\item{ts}{ts.}
}
\value{
A chrDf.
}
\description{
Get a chromPeak from XcmsExperiment.
}
\examples{
getChromPeak(data, cpid = cpid)
}
