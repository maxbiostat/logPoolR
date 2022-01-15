#' Obtain the parameters of the pooled distribution (Beta and Gamma).
#'
#' @param alpha vector of weights, which must lie in a simplex (can be a scalar).
#' @param a vector of 'a' parameters, same size as 'alpha' and 'b'.
#' @param b vector of 'b' parameters.
#'
#' @return a vector with two entries, a_star and b_star.
#' @export pool_par
#'
pool_par <- function(alpha, a, b){
  if(any(alpha <0)) stop("weights must be non-negative")
  if(!all.equal(sum(alpha), 1)) stop("weights do not sum to 1")
  return(
    c(crossprod(a, alpha), crossprod(b, alpha))
  )
}

#' Obtain the parameters of the pooled distribution (Gaussian).
#'
#' @param alpha vector of weights, which must lie in a simplex (can be a scalar).
#' @param m vector of means, same dimension as 'alpha' and 'v'.
#' @param v vector of *variances*. 
#'
#' @return pooled mean and *STANDARD DEVIATION*.
#' @export
#'
pool_par_gauss <- function(alpha, m, v){
  ws <- alpha/v
  vstar <-  1/sum(ws)
  mstar <- sum(ws*m) * vstar
  c(mstar, sqrt(vstar))
}

#' Obtain the pooled density
#'
#' @param x value at which to evaluate the density.
#' @param D list containing the density functions.
#' @param alpha vector of the same size as the length of 'D' with weights.
#' @param log logical. Whether to return the log-density.
#'
#' @return the log-pooled probability density/mass at 'x'.
#' @export dpool
#'
dpool <- function(x, D, alpha, log = FALSE){
  if(any(alpha<0)) { stop ("All weights should be non-negative")}
  if(!identical(sum(alpha), 1, num.eq = TRUE)) {
    stop ("Weigths should sum up to 1")
  }
  log_dens <- function(x){
    sum(unlist(lapply(D, function(d) log(d(x)) )) * alpha)
  } 
  ans <- sapply(x,  function(x) sum(log_dens(x)))
  if(!log) ans <- exp(ans)
  return(ans)
}

#' Normalising constant of the pooled distribution.
#'
#' @param alpha vector of weights, which must lie in a simplex (can be a scalar).
#' @param D list containing the density functions.
#' @param trapez logical. Whether to use the trapezoid method of integration.
#' @param lwr lower limit of integration (exact).
#' @param upr upper limit of integration (exact).
#' @param lwr_num approximate lower limit of integration if 'trapezoid' is `TRUE`.
#' @param upr_num approximate upper limit of integration if 'trapezoid' is `TRUE`.
#'
#' @return numeric. The normalising constant.
#' @export tpool
#'
tpool <- function(alpha, D, trapez = TRUE,
                  lwr = -Inf, upr = Inf,
                  lwr_num = -1E5, upr_num =  1E5){
  if(trapez){
    the_int <- caTools::trapz(x = seq(lwr_num, upr_num, length.out = 1000L),
                              y = dpool(x = seq(lwr_num, upr_num,
                                                length.out = 1000L),
                                        D = D, alpha = alpha))
  }else{
    # proper integration
    integrand <- function(x){
      dpool(x, D = D, alpha = alpha)
    }
    the_int <- stats::integrate(integrand, lwr, upr)$value 
  }
  return(1/the_int)
}

#' Normalising constant of the pooled distributions with support on (0, Inf).
#'
#' @param alpha vector of weights, which must lie in a simplex (can be a scalar).
#' @param D list containing the density functions.
#' @param trapez logical. Whether to use the trapezoid method of integration.
#' @param lwr_num approximate lower limit of integration if 'trapezoid' is `TRUE`.
#' @param upr_num approximate upper limit of integration if 'trapezoid' is `TRUE`.
#'
#' @return numeric. The normalising constant.
#' @export tpool_positive
#'
tpool_positive <- function(alpha, D, trapez = FALSE,
                           lwr_num = 0, upr_num =  1E5){
  ct <- tpool(alpha = alpha,
                   D = D,
                   trapez = trapez,
                   lwr = 0, upr = Inf,
                   lwr_num = lwr_num, upr_num = upr_num)
  return(ct)
}

#' Normalising constant of the pooled distributions with support on (0, 1).
#'
#' @param alpha vector of weights, which must lie in a simplex (can be a scalar).
#' @param D list containing the density functions.
#' @param lwr_num approximate lower limit of integration if 'trapezoid' is `TRUE`.
#' @param upr_num approximate upper limit of integration if 'trapezoid' is `TRUE`.
#'
#' @return numeric. The normalising constant.
#' @export tpool_unit
#'
tpool_unit <- function(alpha,
                       D,
                       lwr_num = 1E-10,
                       upr_num = 1-1E-10){
  ct <- tpool(alpha = alpha, D = D,
                   trapez = FALSE,
                   lwr = 0, upr = 1,
                   lwr_num = lwr_num, upr_num = upr_num)
  return(ct)
}

#' Normalised pooled density.
#'
#' @param x numeric. values at which to evaluate the density.
#' @param D list containing the density functions.
#' @param alpha vector of weights, which must lie in a simplex (can be a scalar).
#' @param log logical. Whether to return the log-density.
#'
#' @return numeric. the normalised (log) density.
#' @export dpoolnorm
#'
dpoolnorm <- function(x, D, alpha, log = FALSE){
  ldens <- dpool(x = x, D = D, alpha = alpha, log = TRUE)
  ans <- log(tpool(alpha = alpha, D = D)) + ldens
  if(!log) ans <- exp(ans)
  return(ans)
}

#' Normalised pooled density with support on (0, Inf).
#'
#' @param x numeric. values at which to evaluate the density.
#' @param D list containing the density functions.
#' @param alpha vector of weights, which must lie in a simplex (can be a scalar).
#' @param log logical. Whether to return the log-density.
#'
#' @return numeric. the normalised (log) density.
#' @export dpoolnorm_positive
#'
dpoolnorm_positive <- function(x, D, alpha, log = FALSE){
  ldens <- dpool(x = x, D = D, alpha = alpha, log = TRUE)
  ans <- log(tpool_positive(alpha = alpha, D = D)) + ldens
  if(!log) ans <- exp(ans)
  return(ans)
}

#' Normalised pooled density with support on (0, 1).
#'
#' @param x numeric. values at which to evaluate the density.
#' @param D list containing the density functions.
#' @param alpha vector of weights, which must lie in a simplex (can be a scalar).
#' @param log logical. Whether to return the log-density.
#'
#' @return numeric. the normalised (log) density.
#' @export dpoolnorm_unit
#'
#' @examples
#' opinions <- list(
#' f0 = function(x) dbeta(x = x, shape1 = 1/2, shape2 = 1/2),
#' f1 = function(x) dbeta(x = x, shape1 = 3, shape2 = 3),
#'f2 = function(x) dbeta(x = x, shape1 = 9, shape2 = 1)
#')
#'set.seed(666)
#' weights.raw <- rgamma(3, 1, 1)
#' weights <- weights.raw/sum(weights.raw)
#' par.star <- pool_par(alpha = weights,
#'                  a = c(1/2, 3, 9),
#'                  b = c(1/2, 3, 1))
#'curve(dbeta(x, shape1 = par.star[1], shape2 = par.star[2]),
#'    ylab = "Density",
#'      lwd = 4, col = "black", lty = 2)
#' curve(dpoolnorm_unit(x = x,  D = opinions, alpha = weights),
#'      add = TRUE, lwd = 2, col = "red")
#'legend(x = "bottom",
#'       legend = c("exact solution", "dpoolnorm"),
#'       col = c("black", "red"),
#'       lty = c(2, 1),
#'       lwd = 2,
#'       bty = 'n')
#'       
dpoolnorm_unit <- function(x, D, alpha, log = FALSE){
  ldens <- dpool(x = x, D = D, alpha = alpha, log = TRUE)
  ans <- log(tpool_unit(alpha = alpha, D = D)) + ldens
  if(!log) ans <- exp(ans)
  return(ans)
}