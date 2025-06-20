spp.stg3.mcmc.nb <- function(out){
  
  spp.loglik <- function(X.beta.sum, lam.int, theta, ds, n){
    theta.lam.int <- theta*lam.int
    n*log(theta) + X.beta.sum - theta.lam.int - lfactorial(n)
  }
  
  ###
  ###  Setup Variables
  ###
  
  n <- out$n
  n.mcmc <- out$n.mcmc
  X.full <- out$X.full
  p <- dim(X.full)[2]
  ds <- out$ds
  lam.int.save <- out$lam.int.save
  X.beta.sum.save <- out$X.beta.sum.save
  mu.beta <- out$mu.beta
  sigma.beta <- out$sigma.beta
  accept <- 0
  
  # beta prior hyperparameters
  mu.0 <- rep(0,p)
  sig.0 <- 10000*diag(p)
  
  # gamma prior hyperparameters
  a <- 0.0000000001
  b <- 0.0000000001
  
  beta.0.save <- rep(0, n.mcmc)
  beta.save <- matrix(0, p, n.mcmc)
  
  # starting values
  beta <- c(out$beta.save[,n.mcmc])
  lam.int <- lam.int.save[n.mcmc]
  X.beta.sum <- X.beta.sum.save[[n.mcmc]]
  
  for(k in 1:n.mcmc){
    if(k %% 1000 == 0){cat(k, " ")}
    
    ###
    ###  Propose theta
    ###
    theta <- rgamma(1, a + n, b + lam.int)
    
    ###
    ###  Propose beta
    ###
    
    idx.star <- sample(1:n.mcmc,1)
    beta.star <- c(out$beta.save[,idx.star])
    lam.int.star <- lam.int.save[idx.star]
    X.beta.sum.star <- X.beta.sum.save[[idx.star]]
    
    mh.1 <- spp.loglik(X.beta.sum.star, lam.int.star, theta, ds, n) + 
            sum(mvnfast::dmvn(t(beta.star), mu.0, sig.0, log=TRUE)) + 
            sum(mvnfast::dmvn(t(beta), mu.beta, sigma.beta, log=TRUE))
    
    mh.2 <- spp.loglik(X.beta.sum, lam.int, theta, ds, n) + 
            sum(mvnfast::dmvn(t(beta), mu.0, sig.0, log=TRUE)) + 
            sum(mvnfast::dmvn(t(beta.star), mu.beta, sigma.beta, log=TRUE))
    
    if(exp(mh.1 - mh.2) > runif(1)){
      beta <- beta.star
      lam.int <- lam.int.star
      X.beta.sum <- X.beta.sum.star
      accept <- accept + 1
    }
    
    ###
    ###  Save Samples
    ###
    
    beta.0.save[k] <- log(theta)
    beta.save[,k] <- beta
    
  };cat("\n")
  
  print(accept/n.mcmc)
  
  ###
  ###  Write Output 
  ###
  
  list(beta.save=beta.save, beta.0.save=beta.0.save, n.mcmc=n.mcmc)
  
}
