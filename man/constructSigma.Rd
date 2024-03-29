% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.r
\name{constructSigma}
\alias{constructSigma}
\title{Construct the variance-covariance matrix for a moment-matching logistic normal.}
\usage{
constructSigma(X)
}
\arguments{
\item{X}{a vector of positive  Dirichlet hyperparameters.}
}
\value{
a positive semi-definite covariance matrix.
}
\description{
Construct the variance-covariance matrix for a moment-matching logistic normal.
}
\details{
Take a Dirichlet with parameter \code{X}: this routine builds the
covariance  matrix for a logistic-normal distribution that has the smallest
Kullback-Leibler divergence to the Dirichlet. This is what is meant by
'moment-matching'.
}
