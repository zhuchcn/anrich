% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrich.R
\name{plot.enrichment.ks}
\alias{plot.enrichment.ks}
\title{ecdf plot for enrichment test}
\usage{
\method{plot}{enrichment.ks}(x, level, alt, ...)
}
\arguments{
\item{x}{enrichment.ks. Returned by \code{\link{enrichement_test}}
when test was set to ks.}

\item{level}{character. The level to plot the ecdf plot. Must be in the
group.}

\item{alt}{character. Must be either greater or less.}
}
\value{
a ggplot2 object
}
\description{
This function makes a ecdf plot for enrichment test result
tested using ks.test.
}
\seealso{
\code{\link{enrichment_test}}
}
\author{
Chenghao Zhu
}
