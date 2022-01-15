ml_beta <- function(yi, ni, a, b, log = FALSE){ # Equation (9) in Raftery et al (2007) 
  ans <- lgamma(ni + 1) - lgamma(ni - yi + 1) - lgamma(yi + 1)  +
    lgamma(a + b)-lgamma(a + b + ni) +
    ( lgamma(a+yi)-lgamma(a) ) + ( lgamma(b + ni - yi) -lgamma(b))
  if(!log) ans <- exp(ans)
  return(ans)
}

marginal_likelihood_surface_beta <- function(y, n, av, bv, N = 100){
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
  MLs <- apply(grid, 1,
               function(row) ml_beta(yi = y, ni = n, a = row[1], b = row[2]))
  ML <- matrix(MLs, nrow = N)
  return(list(M = ML, as = as, bs = bs))
}

normal_mean_marg_like <- function(y, sigma,  m, v, log = TRUE){
  s2 <- sum(y^2)
  xbar <- mean(y)
  n <- length(y)
  sigmasq <- sigma^2
  ###
  lt1 <- log(sigma) - (n/2 * log(2*pi*sigmasq) + 1/2 * log(n*v + sigmasq) )
  lt2 <- -s2/(2*sigmasq) - m^2/(2*v)
  lt3 <- ( (v*n^2*xbar^2)/sigmasq + (sigmasq*m^2)/v + 2*n*xbar*m)/(2*(n*v + sigmasq))
  ###
  ans <- lt1 + lt2 + lt3
  if(!log) ans <- exp(ans)
  return(ans)
}

get_normal_weights <- function(obs, sigma, ms, vs){
  K <- length(ms)
  mls <- rep(NA, K)
  for (k in 1:K){
    mls[k] <- normal_mean_marg_like(y = obs, sigma = sigma,
                                    m = ms[k], v = vs[k])
  }
  logS <-  matrixStats::logSumExp(mls)
  return(
    list(
      mls = mls,
      weights = exp(mls-logS)
    )
  )
}