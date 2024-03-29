% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_error_bar_fns.R
\name{PlotSpectraAsSigsWithUncertainty}
\alias{PlotSpectraAsSigsWithUncertainty}
\title{Convert spectra to signatures and then plot mean with "error" bars}
\usage{
PlotSpectraAsSigsWithUncertainty(
  spectra,
  mean.weighted = TRUE,
  conf.int = 0.95,
  num.of.bootstrap.replicates = 10^4,
  title = "Mean.as.signature"
)
}
\arguments{
\item{spectra}{An \code{\link[ICAMS]{ICAMS}} spectrum catalog.
Convert each spectrum to a signature and then compute the
mean.}

\item{mean.weighted}{Logical. Whether to weigh the samples according to the
number of mutations in them to calculate the weighted mean as the consensus
signature. Default is TRUE. If FALSE, then arithmetic mean will be
calculated.}

\item{conf.int}{A number specifying the required confidence interval. The
error bars will be plotted as bootstrap confidence interval for the mean.
If NULL, then use the maximum and minimum value of the spectra to plot
error bars.}

\item{num.of.bootstrap.replicates}{The number of bootstrap replicates to use.
Default is 10000.}

\item{title}{The name of the output signature.}
}
\value{
The mean of the spectra as a signature, the
constituent spectra as signatures, and the y positions of the
arrowheads.
}
\description{
Convert spectra to signatures and then plot mean with "error" bars
}
