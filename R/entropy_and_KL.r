entropy_beta <- function(a, b){
  lbeta(a, b) - (a-1)*digamma(a) - (b-1)*digamma(b) + (a+b-2)*digamma(a+b)
}

optentbeta <- function(alpha, ap, bp){
  entropy_beta(a = sum(alpha*ap), b = sum(alpha*bp))
}

optentbeta_inv <- function(alpha.inv, ap, bp){
  alpha <- alpha_01(alpha.inv)
  -optentbeta(alpha, ap, bp)
}

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

kl_beta <- function(astar, bstar, ai, bi, type = c("pf", "fp")){
  ## if type = pf, computes KL(pi||f),
  ## if type = fp, computes KL(f || pi)
  # where pi ~ Beta(astar, bstar) and f ~ B(ai, bi)
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

optklbeta <- function(alpha, ap, bp, type){
  K <- length(alpha)
  astar <- sum(alpha*ap)
  bstar <- sum(alpha*bp)
  ds <- rep(NA, K) # the distances from each f_i to \pi
  for (i in 1:K){ ds[i] <-  kl_beta(astar = astar, bstar = bstar,
                                    ai = ap[i], bi = bp[i], type = type)} 
  return(ds)
}

optklbeta_inv <- function(alpha.inv, ap, bp, type = "pf"){
  alpha <- alpha_01(alpha.inv)
  sum(optklbeta(alpha, ap, bp, type))
}

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

kl_gamma <- function(astar, bstar, ai, bi, type = c("pf", "fp")){
  ## KL(f||g), where f~Gamma(a0, b0) and g ~ Gamma(a1, b1)
  ## KL(a, b0, a1, b1)
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

optklgamma <- function(alpha, ap, bp, type){
  K <- length(alpha)
  astar <- sum(alpha*ap)
  bstar <- sum(alpha*bp)
  ds <- rep(NA, K) # the distances from each f_i to \pi
  for (i in 1:K){ ds[i] <-  kl_gamma(astar = astar, bstar = bstar,
                                     ai = ap[i], bi = bp[i], type = type)} 
  return(ds)
}

optklgamma_inv <- function(alpha.inv, ap, bp, type = "pf"){
  alpha <- alpha_01(alpha.inv)
  sum(optklgamma(alpha, ap, bp, type))
}

entropy_gauss <- function(m, v){
  .5*log(2*pi*exp(1)*v)
}
### WARNING: there's no need to optimise the entropy, as there is a "closed-form" solution, namely
### picking the distribution with  the largest variance, but we'll indulge for the sake of completeness.
optentgauss <- function(alpha, mp, vp){
  ws <- alpha/vp
  mstar <- sum(ws*mp)/sum(ws)
  vstar <-  1/sum(ws)
  entropy_gauss(m = mstar, v = vstar)
}

optentgauss_inv <- function(alpha.inv, mp, vp){
  alpha <- alpha_01(alpha.inv)
  -optentgauss(alpha, mp, vp)
}

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

optklgauss_inv <- function(alpha.inv, mp, vp, type = "pf"){
  alpha <- alpha_01(alpha.inv)
  sum(optklgauss(alpha, mp, vp, type))
}

entropy_pool <- function(D, alphas){
  expectlog <- function(x) {-log(x) * dpoolnorm_positive(x, D, alphas)}
  # using dpoolnorm_positive DELIBERATELY introduces a bug
  ent <- stats::integrate(expectlog, 0, Inf)$value
  return(ent)
}