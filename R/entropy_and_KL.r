#' Entropy of a Beta distribution
#'
#' @param a Shape parameter. Scalar larger than zero.
#' @param b Shape parameter. Scalar larger than zero.
#'
#' @return the entropy, a scalar.
#' @export entropy_beta
#'
entropy_beta <- function(a, b){
  lbeta(a, b) - (a-1)*digamma(a) - (b-1)*digamma(b) + (a+b-2)*digamma(a+b)
}

#' Loss function for maximising the entropy of a Beta pool.
#'
#' @param alpha vector of weights, which must lie in a simplex (can be a scalar).
#' @param ap vector of shape parameters, all larger than zero.
#' @param bp vector of shape parameters, all larger than zero.
#' @return the entropy of the pool, a scalar.
#' @export optentbeta
#'
optentbeta <- function(alpha, ap, bp){
  entropy_beta(a = sum(alpha*ap), b = sum(alpha*bp))
}

#' Loss function for maximising the entropy of a Beta pool.
#'
#' @param alpha.inv vector of weights mapped onto R^(K-1).
#' @param ap vector of shape parameters, all larger than zero.
#' @param bp vector of shape parameters, all larger than zero.
#' @return the negative of the entropy, a scalar.
#' @export optentbeta_inv
#'
optentbeta_inv <- function(alpha.inv, ap, bp){
  alpha <- alpha_01(alpha.inv)
  -optentbeta(alpha, ap, bp)
}

#' Compute the entropy of Beta distribution on a grid of shape parameters.
#'
#' @param av a vector of shape parameters, all greater than zero.
#' @param bv a vector of shape parameters, all greater than zero.
#' @param N Number of grid points to compute. 
#'
#' @return a list containing a matrix with the entropies computed on a N x N grid, as 
#' well as vectors containing values of `a` and `b` in the range of `av` and `bv`.
#' @export entropy_surface_beta
#'
entropy_surface_beta <- function(av, bv, N = 100){
  # Since the maximum value for astar is max(av) and
  # the same goes for bstar, we can draw a surface for a given set of a's and b's 
  # to look at the face of the entropy surface
  amax <- max(av)
  amin <- min(av)
  bmax <- max(bv)
  bmin <- min(bv)
  as <- seq(amin, amax, length.out = N)
  bs <- seq(bmin, bmax, length.out = N)
  grid <- expand.grid(as, bs)
  Es <- apply(grid, 1, function(row) entropy_beta(a = row[1], b = row[2]))
  ME <- matrix(Es, nrow = N)
  return(list(M = ME, as = as, bs = bs))
}

#' Kullback-Leibler divergence a Beta and the pool.
#'
#' @param astar shape parameter of the pool \code{pi}.
#' @param bstar shape parameter of the pool \code{pi}.
#' @param ai shape parameter of \code{f}.
#' @param bi shape parameter of \code{f}.
#' @param type  if \code{type = pf}, computes KL(pi||f), whereas 
#' \code{type = fp}, computes KL(f || pi).
#' @details  Here pi ~ Beta(astar, bstar) and f ~ Beta(ai, bi).
#'
#' @return the KL divergence. A scalar.
#' @export kl_beta
#' 
kl_beta <- function(astar, bstar, ai, bi, type = c("pf", "fp")){
  if(type == "pf"){
    a1 = astar
    b1 = bstar
    a2 = ai
    b2 = bi  
  }else{
    a1 = ai
    b1 = bi
    a2 = astar
    b2 = bstar  
  }
  res <-  lbeta(a2, b2) - lbeta(a1, b1) + (a1-a2)*digamma(a1) +
    (b1-b2)*digamma(b1)  + (a2 - a1 + b2 - b1)*digamma(a1 + b1)
  return(res)
}

#' Loss function for minimising the sum of KL divergences of a Beta pool
#'
#' @param alpha vector of weights, which must lie in a simplex (can be a scalar).
#' @param ap vector of shape parameters, all larger than zero.
#' @param bp vector of shape parameters, all larger than zero.
#' @param type if \code{type = pf}, computes KL(pi||f), whereas 
#' \code{type = fp}, computes KL(f || pi). See \code{\link[logPooR]{kl_beta}}.
#'
#' @return a vector containing the KL divergences.
#' @export optklbeta
#'
optklbeta <- function(alpha, ap, bp, type){
  K <- length(alpha)
  astar <- sum(alpha*ap)
  bstar <- sum(alpha*bp)
  ds <- rep(NA, K) # the distances from each f_i to \pi
  for (i in 1:K){ ds[i] <-  kl_beta(astar = astar, bstar = bstar,
                                    ai = ap[i], bi = bp[i], type = type)} 
  return(ds)
}

#' Loss function for maximising the entropy of a Beta pool.
#'
#' @param alpha.inv vector of weights mapped onto R^(K-1).
#' @param ap vector of shape parameters, all larger than zero.
#' @param bp vector of shape parameters, all larger than zero.
#' @param type if \code{type = pf}, computes KL(pi||f), whereas 
#' \code{type = fp}, computes KL(f || pi). See \code{\link[logPooR]{kl_beta}}.
#'
#' @return the sum of KL divergences, a scalar.
#' @export optklbeta_inv
#'
optklbeta_inv <- function(alpha.inv, ap, bp, type = "pf"){
  alpha <- alpha_01(alpha.inv)
  sum(optklbeta(alpha, ap, bp, type))
}

#' Entropy of a Gamma distribution
#'
#' @param a Shape parameter. Scalar larger than zero.
#' @param b Rate parameter. Scalar larger than zero.
#'
#' @return the entropy, a scalar.
#' @export entropy_gamma
#'
entropy_gamma <- function(a, b){
  a - log(b) + log(gamma(a)) + (1-a)*digamma(a)
}

entropy_surface_gamma <- function(av, bv, N = 100){
  amax <- max(av)
  amin <- min(av)
  bmax <- max(bv)
  bmin <- min(bv)
  as <- seq(amin, amax, length.out = N)
  bs <- seq(bmin, bmax, length.out = N)
  grid <- expand.grid(as, bs)
  Es <- apply(grid, 1, function(row) entropy_gamma(a = row[1], b = row[2]))
  ME <- matrix(Es, ncol = N)
  return(list(M = ME, as = as, bs = bs))
}

#' Kullback-Leibler divergence a Gamma and the pool.
#'
#' @param astar shape parameter of the pool \code{pi}.
#' @param bstar rate parameter of the pool \code{pi}.
#' @param ai shape parameter of \code{f}.
#' @param bi rate parameter of \code{f}.
#' @param type  if \code{type = pf}, computes KL(pi||f), whereas 
#' \code{type = fp}, computes KL(f || pi).
#' @details  Here pi ~ Gamma(astar, bstar) and f ~ Gamma(ai, bi)
#'
#' @return the KL divergence. A scalar.
#' @export kl_gamma
#'
kl_gamma <- function(astar, bstar, ai, bi, type = c("pf", "fp")){
  if(type == "pf"){
    a0 = astar
    b0 = bstar
    a1 = ai
    b1 = bi  
  }else{
    a0 = ai
    b0 = bi
    a1 = astar
    b1 = bstar  
  }
  ans <- (a0-a1)*digamma(a0) - log(gamma(a0)) + log(gamma(a1)) +
    a1*(log(b0/b1)) + a0*((b1-b0)/b0)
  return(ans)
}

#' Loss function for minimising the sum of KL divergences of a Gamma pool.
#'
#' @param alpha vector of weights, which must lie in a simplex (can be a scalar).
#' @param ap vector of shape parameters, all larger than zero.
#' @param bp vector of rate parameters, all larger than zero.
#' @param type if \code{type = pf}, computes KL(pi||f), whereas 
#' \code{type = fp}, computes KL(f || pi). See  \code{\link[logPooR]{kl_gamma}}.
#'
#' @return a vector containing the KL divergences.
#' @export optklgamma
#'
optklgamma <- function(alpha, ap, bp, type){
  K <- length(alpha)
  astar <- sum(alpha*ap)
  bstar <- sum(alpha*bp)
  ds <- rep(NA, K) # the distances from each f_i to \pi
  for (i in 1:K){ ds[i] <-  kl_gamma(astar = astar, bstar = bstar,
                                     ai = ap[i], bi = bp[i], type = type)} 
  return(ds)
}

#' Loss function for minimising the sum of KL divergences of a Gamma pool.
#'
#' @param alpha.inv vector of weights mapped onto R^(K-1).
#' @param ap vector of shape parameters, all larger than zero.
#' @param bp vector of rate parameters, all larger than zero.
#' @param type if \code{type = pf}, computes KL(pi||f), whereas 
#' \code{type = fp}, computes KL(f || pi). See  \code{\link[logPooR]{kl_gamma}}.
#'
#' @return the sum of KL divergences, a scalar.
#' @export optklgamma_inv
#'
optklgamma_inv <- function(alpha.inv, ap, bp, type = "pf"){
  alpha <- alpha_01(alpha.inv)
  sum(optklgamma(alpha, ap, bp, type))
}

entropy_gauss <- function(m, v){
  .5*log(2*pi*exp(1)*v)
}

#' Loss function for maximising the entropy of a Gaussian pool.
#'
#' @param alpha vector of weights, which must lie in a simplex (can be a scalar).
#' @param mp vector of means.
#' @param vp vector of variances, all larger than zero.
#' @return the entropy of the pool, a scalar.
#' @export optentgauss
#'
optentgauss <- function(alpha, mp, vp){
  ### WARNING: there's no need to optimise the entropy, as there is a "closed-form" solution, namely
  ### picking the distribution with  the largest variance, but we'll indulge for the sake of completeness.
  ws <- alpha/vp
  mstar <- sum(ws*mp)/sum(ws)
  vstar <-  1/sum(ws)
  entropy_gauss(m = mstar, v = vstar)
}

#' Loss function for maximising the entropy of a Gaussian pool.
#'
#' @param alpha.inv vector of weights mapped onto R^(K-1).
#' @param mp vector of means.
#' @param vp vector of variances, all larger than zero.
#' @return the entropy of the pool, a scalar.
#' @export optentgauss_inv
#'
optentgauss_inv <- function(alpha.inv, mp, vp){
  alpha <- alpha_01(alpha.inv)
  -optentgauss(alpha, mp, vp)
}

#' Kullback-Leibler divergence a Gaussian and the pool.
#'
#' @param mstar mean parameter of the pool \code{pi}.
#' @param vstar variance of the pool \code{pi}.
#' @param mi mean parameter of \code{f}.
#' @param vi variance parameter of \code{f}.
#' @param type  if \code{type = pf}, computes KL(pi||f), whereas 
#' \code{type = fp}, computes KL(f || pi).
#' @details  Here pi ~ Normal(mstar, vstar) and f ~ Normal(mi, vi) and parametrisation
#' is mean and variance.
#'
#' @return the KL divergence. A scalar.
#' @export kl_gauss
#'
kl_gauss <- function(mstar, vstar, mi, vi, type = c("pf", "fp")){
  ## WARNING: parametrisation is MEAN and VARIANCE!
  if(type == "pf"){
    m0 = mstar
    v0 = vstar
    m1 = mi
    v1 = vi  
  }else{
    m0 = mi
    v0 = vi
    m1 = mstar
    v1 = vstar  
  }
  ans <- log(sqrt(v1)) - log( sqrt(v0))  + (v0 + (m0-m1)^2)/(2 * v1) - 1/2
  return(ans)
}

#' Loss function for minimising the sum of KL divergences of a Gaussian pool.
#'
#' @param alpha vector of weights, which must lie in a simplex (can be a scalar).
#' @param mp vector of means.
#' @param vp vector of variances, all larger than zero.
#' @param type if \code{type = pf}, computes KL(pi||f), whereas 
#' \code{type = fp}, computes KL(f || pi). See  \code{\link[logPooR]{kl_gauss}}.
#'
#' @return a vector containing the KL divergences.
#' @export optklgauss
#'
optklgauss <- function(alpha, mp, vp, type){
  K <- length(alpha)
  ws <- alpha/vp
  mstar <- sum(ws*mp)/sum(ws)
  vstar <-  1/sum(ws)
  ds <- rep(NA, K) # the distances from each f_i to \pi
  for (i in 1:K){ ds[i] <-  kl_gauss(mstar = mstar, vstar = vstar,
                                     mi = mp[i], vi = vp[i], type = type)} 
  return(ds)
}

#' Loss function for minimising the sum of KL divergences of a Gaussian pool.
#'
#' @param alpha.inv vector of weights mapped onto R^(K-1).
#' @param mp vector of means.
#' @param vp vector of variances, all larger than zero.
#' @param type if \code{type = pf}, computes KL(pi||f), whereas 
#' \code{type = fp}, computes KL(f || pi). See  \code{\link[logPooR]{kl_gauss}}.
#'
#' @return a vector containing the KL divergences.
#' @export optklgauss_inv
#'
optklgauss_inv <- function(alpha.inv, mp, vp, type = "pf"){
  alpha <- alpha_01(alpha.inv)
  sum(optklgauss(alpha, mp, vp, type))
}

#' Compute the entropy of a general pool of distributions
#'
#' @param D list containing the density functions.
#' @param alpha vector of weights, which must lie in a simplex (can be a scalar).
#'
#' @return the entropy, a scalar.
#' @export entropy_pool
#'
entropy_pool <- function(D, alpha){
  expectlog <- function(x) {-log(x) * dpoolnorm_positive(x, D, alpha)}
  # using dpoolnorm_positive DELIBERATELY introduces a bug
  ent <- stats::integrate(expectlog, 0, Inf)$value
  return(ent)
}
