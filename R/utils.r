renorm <- function(x) x + abs(min(x))

invlogit <- function (y){
  return(
    1/(1 + exp(-y))
  )
}

logit <- function(x){
  return(
    log(x) - log1p(-x)
  )
}

#' Maps a (n+1) open simplex to R^n.
#'
#' @param alpha vector of weights that live on an open simplex.
#'
#' @return a projection of 'alpha' onto R^n.
#' @export alpha_real 
#'
alpha_real <- function(alpha){
  p <- length(alpha)
  if(p == 1){
    ans <- logit(alpha)
  }else{
    ans <- log( alpha[-p] / alpha[p]  )
  }
  return(ans)
}

#' Maps from R^n to the open simplex.
#'
#' @param alpha.inv a vector in R^n.
#'
#' @return a vector on the (n+1) open simplex.
#' @export alpha_01
#' @references See \url{https://cran.r-project.org/web/packages/SALTSampler/vignettes/SALTSampler.html}
alpha_01 <- function(alpha.inv){
  K <- length(alpha.inv) + 1
  z <- rep(NA, K-1)
  alphas <- rep(0, K)
  for(k in 1:(K-1)){
    z[k] <- invlogit(alpha.inv[k] + log( 1/(K-k) ))
    alphas[k] <- (1 - sum(alphas[1:(k-1)])) * z[k]
  }
  alphas[K] <- 1-sum(alphas[-K])
  return(alphas)
}

#' Construct the variance-covariance matrix for a moment-matching logistic normal.
#'
#' @param X a vector of positive  Dirichlet hyperparameters.
#'
#' @return a positive semi-definite covariance matrix.
#' @export constructSigma
#' @details  Take a Dirichlet with parameter \code{X}: this routine builds the
#'  covariance  matrix for a logistic-normal distribution that has the smallest
#'   Kullback-Leibler divergence to the Dirichlet. This is what is meant by
#'    'moment-matching'.
constructSigma <- function(X){
  K <- length(X)
  Sigma <- matrix(rep(trigamma(X[K]), K^2), ncol = K, nrow = K)
  diag(Sigma) <- trigamma(X) + trigamma(X[K])
  return(Sigma)
}

#' Generate samples from a multivariate logistic normal.
#'
#' @param N  number of draws.
#' @param m is a vector of means. If a scalar is given, it is recycled \code{K} times.
#' @param Sigma the variance-covariance matrix.
#'
#' @return
#' @export rlogisticnorm
#'
rlogisticnorm <- function(N, m, Sigma){
  raw.draws <-  mvtnorm::rmvnorm(n = N, mean = m, sigma = Sigma)
  logisticNorm <- function(x){
    D <- length(x)
    y <- rep(NA, D)
    for(d in 1:(D-1)) y[d] <- exp(x[d])/(1 + sum(exp(x[-D])) )
    y[D] <- 1/(1 + sum(exp(x[-D])) )
    return(y)
  } 
  draws <- apply(raw.draws, 1, logisticNorm)
  return(t(draws))
}

#' Return mean and 100*alpha% quantiles.
#'
#' @param x vector of values to be summarised.
#' @param alpha probability level. A number between 0 and 1.
#' @param na.rm logical.  If \code{TRUE} \code{NA}s are excluded.
#'
#' @return a data.frame with mean and upper and lower quantiles.
#' @export mean_ci
#'
mean_ci <- function(x, alpha = .95, na.rm = TRUE){
  qs <- stats::quantile(x,
                 probs = c((1 - alpha), (1 + alpha) )/2,
                 na.rm = na.rm)
  return(data.frame(
    mean = mean(x, na.rm = na.rm),
    lwr = as.numeric(qs[1]),
    upr = as.numeric(qs[2])
  ))
}

get_ratio <- function(x) {
  ## outputs the ratio of the largest and second largest values in a vector x
  y <- x[order(x, decreasing = TRUE)]
  y[1]/y[2]
}

get_ratio_correctExpert <- function(x, p = 3) {
  ## outputs the ratio of the largest and second largest values in a vector x
  y <- x[-p][order(x[-p], decreasing = TRUE)]
  x[p]/y[1]
}

stat_beta <- function(a, alpha = .95){
  qs <- stats::qbeta(c((1-alpha)/2, (1+alpha)/2), a[1], a[2])
return(data.frame(
  mean = a[1] / (a[1] + a[2]),
  lwr = as.numeric(qs[1]),
  upr = as.numeric(qs[2])
))
}

beta_mode <- function(a, b) (a-1)/(a +b -2)
beta_mean <- function(a, b) a/ (a + b)
beta_sd <- function(a, b) sqrt( (a*b)/( (a + b + 1) * (a + b)^2 ))
beta_var <- function(a, b) (beta_sd(a, b) )^2

fbeta <- function(x, par){
  stats::dbeta(x, par[1], par[2])
}

stat_gamma <- function(a, alpha = .95){
  qs <- stats::qgamma( c((1-alpha)/2, (1 + alpha)/2), a[1], a[2])
  return(data.frame(
    mean = a[1]/a[2],
    lwr = as.numeric(qs[1]),
    upr = as.numeric(qs[2])
  ))
}

fgamma <- function(x, par){
  stats::dgamma(x, par[1], par[2])
}

stat_gauss <- function(p, alpha = .95){
  qs <- stats::qnorm( c((1-alpha)/2, (1+alpha)/2), p[1], p[2])
  return(data.frame(
    mean = p[1],
    lwr = as.numeric(qs[1]),
    upr = as.numeric(qs[2])
  ))
}

