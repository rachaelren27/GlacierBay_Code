log_grad <- function(beta, X, y, n, p, mu.beta, Sigma.beta.inv){
  exp_X_beta <- exp(X %*% beta)
  
  logit <- exp_X_beta / (1 + exp_X_beta)
  
  term1 <- t(y - logit) %*% X
  
  return(term1 - t(beta - mu.beta) %*% Sigma.beta.inv)
}

SVGD_update <- function(beta.particles, X, y, mu.beta, Sigma.beta.inv, epsilon){
  n.particles <- ncol(beta.particles)
  n <- length(y)
  p <- ncol(X)
  
  pairwise.dist <- as.matrix(dist(t(beta.particles)))
  h <- (median(pairwise.dist)^2) / log(n.particles)
  kernel.matrix <- exp(-pairwise.dist^2 / h)
  
  log.grads <- apply(beta.particles, 2, function(beta) {
    log_grad(beta, X, y, n, p, mu.beta, Sigma.beta.inv)
  })
  
  beta.updates <- matrix(0, nrow = p, ncol = n.particles)
  
  for (i in 1:n.particles) {
    beta.i <- beta.particles[, i]
    
    kernel.grads <- sweep(beta.particles, 2, beta.i, `-`) * (-2 / h) * kernel.matrix[i, ]
    
    beta.updates[, i] <- rowSums(
      sweep(log.grads, 2, kernel.matrix[i, ], `*`) + kernel.grads
    )
  }
  
  beta.particles <- beta.particles + (epsilon / n.particles) * beta.updates
  
  return(beta.particles)
}

SVGD <- function(beta.init, epsilon, n.iter, X, y, mu.beta, Sigma.beta.inv){
  SVGD.out <- list()
  SVGD.out[[1]] <- beta.init
  for(k in 2:n.iter){
    beta.particles <- SVGD.out[[k-1]]
    SVGD.out[[k]] <- SVGD_update(beta.particles, X.pg, y.binary, mu.beta,
                                 Sigma.beta.inv, epsilon)
    if(k %% 100 == 0){
      print(k)
    }
  }
  
  return(SVGD.out)
}