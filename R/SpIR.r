makedf <- function(mc.samples){
  kde <- stats::density(mc.samples)
  return(Vectorize(stats::approxfun(kde)))
}

# Logarithmic pooling (LP) via Sampling-importance-resampling (SpIR)
LP_SpIR <- function(k, l, Model, rq1, dq2, dL1, dL2, alpha, cores = 4){
  ## currently tuned to perform LP on two dimensions
  theta.samp <- rq1(n = k)
  phi.transf <- unlist(
    parallel::mclapply(1:nrow(theta.samp), function(i) Model(theta.samp[i, ]), mc.cores = cores) 
  )
  q1star <- makedf(phi.transf)
  get_Weight <- function(theta, phi, log = FALSE){
    # (dq2(phi)/q1star(phi))^{1-alpha} * dL1(theta) * dL2(phi)
    ans <- (1-alpha)*{log(dq2(phi))-log(q1star(phi)) } + log(dL1(theta)) + log(dL2(phi))  
    if(!log) ans <- exp(ans)
    return(ans)
  }
  get_Weight <- compiler::cmpfun(get_Weight)
  ws <- unlist(
    parallel::mclapply(seq_len(nrow(theta.samp)), function(i) {
      get_Weight(theta = theta.samp[i, ], phi = phi.transf[i])
    }, mc.cores = cores)
  ) 
  ws[which(ws == -Inf)] <- 0 ## giving  weight = 0 to "impossible" values
  ws[which(ws == "NaN")] <- 0 ## giving  weight = 0 to "weird" values ;-)
  resamp.Pos <-  sample(seq_len(nrow(theta.samp)), size = l,
                        replace = TRUE, prob = ws/sum(ws))
  return(list(
    theta.resamp = theta.samp[resamp.Pos, ],
    phi.resamp = phi.transf[resamp.Pos])
  )
}

LPSpIR_varying_alpha <- function(k, l, Model, rq1, dq2, dL1, dL2, cores = 4){
  ## currently tuned to perform LP on two dimensions
  theta.samp <- rq1(n = k)
  phi.transf <- unlist(
    parallel::mclapply(1:nrow(theta.samp), 
                       function(i) Model(theta.samp[i, 1:2]), mc.cores = cores) 
  )
  q1star <- makedf(phi.transf)
  corrected_q1star <-  function(x){
    res <- q1star(x)
    res[is.na(res)] <- 0
    return(res)
  } 
  getKa <- function(a){
    tpool(alpha = c(a, 1-a),
          D = list(
            f0 = function(x) corrected_q1star(x),
            f1 = function(x) stats::dnorm(x, mean = 7800, sd = 1300) )
    )
  }
  getKa <- compiler::cmpfun(getKa)
  #
  get_Weight <- function(theta, phi, log = FALSE){
    # getKa(theta[3]) * (dq2(phi)/q1star(phi))^{1-theta[3]} * dL1(theta) * dL2(phi)
    ans <- log(getKa(theta[3])) + (1-theta[3])*{log(dq2(phi)) - log(q1star(phi))} + log(dL1(theta)) + log(dL2(phi))
    if(!log) ans <- exp(ans)
    return(ans)
  }
  get_Weight <- compiler::cmpfun(get_Weight)
  ws <- unlist(
    parallel::mclapply(seq_len(nrow(theta.samp)), function(i) {
      get_Weight(theta = theta.samp[i, ], phi = phi.transf[i])
    }, mc.cores = cores)
  ) 
  ws[which(ws == -Inf)] <- 0 ## giving  weight = 0 to "impossible" values
  ws[which(ws == "NaN")] <- 0 ## giving  weight = 0 to "weird" values ;-)
  resamp.Pos <-  sample(seq_len(nrow(theta.samp)), size = l,
                        replace = TRUE, prob = ws/sum(ws))
  return(list(
    theta.resamp = theta.samp[resamp.Pos, ],
    phi.resamp = phi.transf[resamp.Pos])
  )
}