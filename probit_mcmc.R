probit_mcmc <- function(y, X,
                        mu.beta, sigma.beta,
                        n.mcmc){

  n <- length(y)
  p <- length(mu.beta)
  
  # initialization
  beta <- mu.beta
  
  beta.save <- matrix(nrow = n.mcmc, ncol = p)
  
  A <- crossprod(X,X) + diag(p)
  A.inv <- solve(A)
  mu.sigma.beta <- crossprod(mu.beta, sigma.beta)

  for(k in 1:n.mcmc){
    # update z
    z.means <- X %*% beta
    
    z <- sapply(1:length(y), function(i) {
      if (y[i] == 1) {
        rtruncnorm(1, a = 0, b = Inf, mean = z.means[i])
      } else {
        rtruncnorm(1, a = -Inf, b = 0, mean = z.means[i])
      }
    })
    
    # update beta
    b <- crossprod(z, X) + mu.sigma.beta
    beta <- t(rmvn(1, mu = tcrossprod(A,b), sigma = A.inv))
    
    beta.save[k,] <- beta
    
    if (k %% 100 == 0) {cat(k, " ")}
  
  }
  
}