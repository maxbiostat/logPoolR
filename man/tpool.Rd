% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pooling.r
\name{tpool}
\alias{tpool}
\title{Normalising constant of the pooled distribution.}
\usage{
tpool(
  alpha,
  D,
  trapez = TRUE,
  lwr = -Inf,
  upr = Inf,
  lwr_num = -1e+05,
  upr_num = 1e+05
)
}
\arguments{
\item{alpha}{vector of weights, which must lie in a simplex (can be a scalar).}

\item{D}{list containing the density functions.}

\item{trapez}{logical. Whether to use the trapezoid method of integration.}

\item{lwr}{lower limit of integration (exact).}

\item{upr}{upper limit of integration (exact).}

\item{lwr_num}{approximate lower limit of integration if 'trapezoid' is \code{TRUE}.}

\item{upr_num}{approximate upper limit of integration if 'trapezoid' is \code{TRUE}.}
}
\value{
numeric. The normalising constant.
}
\description{
Normalising constant of the pooled distribution.
}
