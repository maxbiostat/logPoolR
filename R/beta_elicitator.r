elicit_beta_mean_cv <- function(m0, v0 = NULL, cv = 1){
  if(!is.null(v0)){
    a <- -(m0*v0 + m0^3 - m0^2)/v0
    b <- ((m0-1)*v0 + m0^3 - 2*m0^2 + m0)/v0
  }else{
    a <- -(m0*(cv*m0)^2 + m0^3 - m0^2)/(cv*m0)^2
    b <- ((m0-1)*(cv*m0)^2 + m0^3 - 2*m0^2 + m0)/(cv*m0)^2
  }
  if(a < 0 || b <0) warning("Warning: at least one of the obtained parameters is not valid")
  return(list(a = a, b = b))
}

elicit_beta_median_iq <- function(m, d, q = .90){
  u <- m + d
  loss <- function(x){
    a <- x[1]
    b <- x[2]
    m.hat <- stats::qbeta(.5, shape1 = a, shape2 = b)
    u.hat <- stats::qbeta(q, shape1 = a, shape2 = b)
    error <-  (m.hat - m)^2 + (u.hat-u)^2
    return(error)  
  }
  opt <- suppressWarnings( stats::optim(loss, par = c(1, 1), lower = c(1E-3, 1E-3)) )
  a <- opt$par[1]
  b <- opt$par[2]
  if(a < 0 || b < 0) warning("Warning: at least one of the obtained parameters is not valid")
  return(list(a = a, b = b))
}

elicit_beta_quantiles <- function(l, u, alpha = .95){
  q0 <- (1-alpha)/2
  q1 <- (1 + alpha)/2
  loss <- function(x){
    a <- x[1]
    b <- x[2]
    l.hat <- stats::qbeta(q0, shape1 = a, shape2 = b)
    u.hat <- stats::qbeta(q1, shape1 = a, shape2 = b)
    # error <- (l.hat - l)^2 + (u.hat-u)^2
    error <- abs(l.hat - l) + abs(u.hat-u) ## L1 norm: better?
    return(error)  
  }
  opt <- suppressWarnings( stats::optim(loss, par = c(1, 1), lower = c(1E-3, 1E-3)) )
  a <- opt$par[1]
  b <- opt$par[2]
  if(a < 0 || b < 0){
    warning("Warning: at least one of the obtained parameters is not valid")
  } 
  return(list(a = a, b = b))
}