% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pooling.r
\name{dpoolnorm_unit}
\alias{dpoolnorm_unit}
\title{Normalised pooled density with support on (0, 1).}
\usage{
dpoolnorm_unit(x, D, alpha, log = FALSE)
}
\arguments{
\item{x}{numeric. values at which to evaluate the density.}

\item{D}{list containing the density functions.}

\item{alpha}{vector of weights, which must lie in a simplex (can be a scalar).}

\item{log}{logical. Whether to return the log-density.}
}
\value{
numeric. the normalised (log) density.
}
\description{
Normalised pooled density with support on (0, 1).
}
\examples{
opinions <- list(
f0 = function(x) dbeta(x = x, shape1 = 1/2, shape2 = 1/2),
f1 = function(x) dbeta(x = x, shape1 = 3, shape2 = 3),
f2 = function(x) dbeta(x = x, shape1 = 9, shape2 = 1)
)
set.seed(666)
weights.raw <- rgamma(3, 1, 1)
weights <- weights.raw/sum(weights.raw)
par.star <- pool_par(alpha = weights,
                 a = c(1/2, 3, 9),
                 b = c(1/2, 3, 1))
curve(dbeta(x, shape1 = par.star[1], shape2 = par.star[2]),
   ylab = "Density",
     lwd = 4, col = "black", lty = 2)
curve(dpoolnorm_unit(x = x,  D = opinions, alpha = weights),
     add = TRUE, lwd = 2, col = "red")
legend(x = "bottom",
      legend = c("exact solution", "dpoolnorm"),
      col = c("black", "red"),
      lty = c(2, 1),
      lwd = 2,
      bty = 'n')
      
}
