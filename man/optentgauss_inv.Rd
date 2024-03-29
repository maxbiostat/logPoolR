% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entropy_and_KL.r
\name{optentgauss_inv}
\alias{optentgauss_inv}
\title{Loss function for maximising the entropy of a Gaussian pool.}
\usage{
optentgauss_inv(alpha.inv, mp, vp)
}
\arguments{
\item{alpha.inv}{vector of weights mapped onto R^(K-1).}

\item{mp}{vector of means.}

\item{vp}{vector of variances, all larger than zero.}
}
\value{
the entropy of the pool, a scalar.
}
\description{
Loss function for maximising the entropy of a Gaussian pool.
}
