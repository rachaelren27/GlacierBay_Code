spp.comp.mcmc <- function(s.mat,X,X.full,ds,n.mcmc,theta.tune,beta.tune){
  
  ###
  ###  Subroutine 
  ###
  
  spp.loglik <- function(theta,beta,X,X.full,ds,n){
    llam <- log(theta) + X%*%beta
    lam.int <- sum(theta*exp(log(ds) + X.full%*%beta))
    sum(llam) - lam.int - lfactorial(n)
  }
  
  ###
  ###  Set up variables
  ###
  
  n <- nrow(s.mat)
  p <- ncol(X)
  
  beta.save <- matrix(0,p,n.mcmc)
  # beta.0.save=rep(0,n.mcmc)
  theta.save <- rep(0, n.mcmc)
  # D <- rep(0, n.mcmc)
  
  ###
  ###  Priors and Starting Values
  ###
  
  mu.00 <- 0
  sig.00 <- 10
  a <- 0.0001
  b <- 0.0001
  mu.0 <- rep(0,p)
  sig.0 <- 100*diag(p)
  
  theta <- 1
  beta <- rep(0,p)
  
  # beta.0.tune=.1
  # beta.tune=.1
  
  ###
  ###  MCMC Loop 
  ###
  
  accept.theta <- 0
  accept.beta <- 0
  
  for(k in 1:n.mcmc){
    if(k%%1000 == 0){cat(k," ")}
    
    ###
    ###  Sample beta 
    ###
    
    beta.star <- t(mvnfast::rmvn(1, beta, beta.tune*diag(p)))
    
    mh.1 <- spp.loglik(theta, beta.star, X, X.full, ds, n) +
         sum(mvnfast::dmvn(t(beta.star), mu.0, sig.0, log=TRUE))
    
    mh.2 <- spp.loglik(theta, beta, X, X.full, ds, n) +
         sum(mvnfast::dmvn(t(beta), mu.0, sig.0, log=TRUE))
    
    if(exp(mh.1 - mh.2) > runif(1)){
      beta <- beta.star
      accept.beta <- accept.beta + 1
    }
    
    # beta.star=rnorm(p,beta,beta.tune)
    # mh.1=spp.loglik(theta,beta.star,X,X.full,ds,n)+sum(dnorm(beta.star,mu.0,sig.0,log=TRUE))
    # mh.2=spp.loglik(theta,beta,X,X.full,ds,n)+sum(dnorm(beta,mu.0,sig.0,log=TRUE))
    # if(exp(mh.1-mh.2)>runif(1)){
    #   beta=beta.star 
    # }
    
    ###
    ###  Sample beta.0
    ###
    
    # beta.0.star=rnorm(1,beta.0,beta.0.tune)
    # mh.1=spp.loglik(beta.0.star,beta,X,X.full,ds,n)
    # + dnorm(beta.0.star, mu.00, sig.00, log = TRUE)
    # mh.2=spp.loglik(beta.0,beta,X,X.full,ds,n)
    # + dnorm(beta.0, mu.00, sig.00, log = TRUE)
    # if(exp(mh.1-mh.2)>runif(1)){
    #   beta.0 <- beta.0.star
    # }
    
    theta.star <- rlnorm(1, log(theta), theta.tune)
    mh.1 <- spp.loglik(theta.star, beta, X, X.full, ds, n)
      + dgamma(theta.star, a, b, log = TRUE)
      + dlnorm(log(theta), theta.tune, log = TRUE)
    mh.2 <- spp.loglik(theta, beta, X, X.full, ds, n)
      + dgamma(theta, a, b, log = TRUE)
      + dlnorm(log(theta.star), theta.tune, log = TRUE)
    if(exp(mh.1 - mh.2) > runif(1)){
      theta <- theta.star
      accept.theta <- accept.theta + 1
    }
    
    ###
    ###  Save Samples
    ###
    
    # beta.0.save[k]=beta.0 
    theta.save[k] <- theta
    beta.save[,k] <- beta 
    # D[k] = -2*spp.loglik(theta, beta, X, X.full, ds, n)
    
  };cat("\n")
  
  ###
  ###  Write Output
  ###
  
  print(paste("theta acceptance ratio:", accept.theta/n.mcmc))
  print(paste("beta acceptance ratio:", accept.beta/n.mcmc))
  
  list(beta.save=beta.save, beta.0.save=log(theta.save), n.mcmc=n.mcmc)
  # beta.0.save = log(theta.save)/beta.0.save
  
}