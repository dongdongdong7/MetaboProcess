% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Paramounter_function.R
\name{paramounterPart2}
\alias{paramounterPart2}
\title{paramounterPart2}
\usage{
paramounterPart2(
  file_path,
  massSDrange = 2,
  smooth = 0,
  ppmCut = NA,
  thread = 1,
  msLevel = 1L
)
}
\arguments{
\item{file_path}{Files path of QC sample.}

\item{massSDrange}{massSDrange.}

\item{smooth}{smooth.}

\item{ppmCut}{ppmCut calculated by part1.}

\item{thread}{thread.}

\item{msLevel}{msLevel.}
}
\value{
A vector of XCMS parameters.
}
\description{
Paramounter function part2 for xcms parameters.
}
\examples{
parameters <- paramounterPart2(file_path = file_path, ppmCut = ppmCut, thread = 16)
}
