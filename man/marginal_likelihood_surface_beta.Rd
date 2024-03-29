% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marginal_likelihoods.r
\name{marginal_likelihood_surface_beta}
\alias{marginal_likelihood_surface_beta}
\title{Compute marginal likelihoods for a Beta-binomial model on a grid}
\usage{
marginal_likelihood_surface_beta(y, n, av, bv, N = 100)
}
\arguments{
\item{y}{number of successes}

\item{n}{number of trials.}

\item{av}{vector of shape1 hyperparameters.}

\item{bv}{vector of shape2 hyperparameters.}

\item{N}{number of grid points to be considered.}
}
\value{
a list containing a matrix with the marginal likelihoods computed on a N x N grid, as
well as vectors containing values of \code{a} and \code{b} in the range of \code{av} and \code{bv}.
}
\description{
Compute marginal likelihoods for a Beta-binomial model on a grid
}
