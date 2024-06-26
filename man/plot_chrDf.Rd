% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_function.R
\name{plot_chrDf}
\alias{plot_chrDf}
\title{plot_chrDf}
\usage{
plot_chrDf(chrDf, linewidth = 1, noise = NA, xlim = NA, baseline = FALSE)
}
\arguments{
\item{chrDf}{chrDf.}

\item{linewidth}{linewidth.}

\item{noise}{noise.}

\item{xlim}{xlim.}

\item{baseline}{Whether to draw the baseline, if TRUE, you need to ensure that there is a baseline in chrDf.}
}
\value{
ggplot object.
}
\description{
Plot a chromatogram using chromatographic data frame.
}
\examples{
plot_chrDf(chrDf = chrDf_i, linewidth = 2)
}
