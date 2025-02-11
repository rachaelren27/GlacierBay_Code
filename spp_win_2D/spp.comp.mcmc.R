spp.comp.mcmc <- function(s.mat,X,X.full,ds,n.mcmc,theta.tune,beta.tune){
  # theta.tune/beta.0.tune
  
  #
  #  1-D SPP w/ complete data (s and n) likelihood and windowing
  #
  
  ###
  ###  Subroutine 
  ###
  
  spp.loglik <- function(theta,beta,X,X.full,ds,n){ # pass in X.full, ds
    llam=log(theta)+X%*%beta
    lam.int=sum(theta*exp(log(ds)+X.full%*%beta))
    sum(llam)-lam.int-lfactorial(n) # check if we can delete lfactorial(n)
  }
  
  ###
  ###  Set up variables
  ###
  
  n=nrow(s.mat)
  p=dim(X)[2]
  
  beta.save=matrix(0,p,n.mcmc)
  # beta.0.save=rep(0,n.mcmc)
  theta.save <- rep(0, n.mcmc)
  # D <- rep(0, n.mcmc)
  
  ###
  ###  Priors and Starting Values
  ###
  
  mu.00=0
  sig.00=10
  a <- 0.0001
  b <- 0.0001
  mu.0=rep(0,p)
  sig.0=100*diag(p)
  
  theta <- exp(4)
  beta=rep(0,p)
  
  # beta.0.tune=.1
  # beta.tune=.1
  
  ###
  ###  MCMC Loop 
  ###
  
  for(k in 1:n.mcmc){
    if(k%%1000==0){cat(k," ")}
    
    ###
    ###  Sample beta 
    ###
    
    beta.star=t(mvnfast::rmvn(1,beta,beta.tune*diag(p)))
    
    mh.1=spp.loglik(theta,beta.star,X,X.full,ds,n) +
         sum(mvnfast::dmvn(t(beta.star),mu.0,sig.0,log=TRUE))
    
    mh.2=spp.loglik(theta,beta,X,X.full,ds,n) +
         sum(mvnfast::dmvn(t(beta),mu.0,sig.0,log=TRUE))
    
    if(exp(mh.1-mh.2)>runif(1)){
      beta=beta.star
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
    mh.1=spp.loglik(theta.star,beta,X,X.full,ds,n)
    + dgamma(theta.star, a, b, log = TRUE)
    + dlnorm(log(theta), theta.tune, log = TRUE)
    mh.2=spp.loglik(theta,beta,X,X.full,ds,n)
    + dgamma(theta, a, b, log = TRUE)
    + dlnorm(log(theta.star), theta.tune, log = TRUE)
    if(exp(mh.1-mh.2)>runif(1)){
      theta <- theta.star
    }
    
    ###
    ###  Save Samples
    ###
    
    # beta.0.save[k]=beta.0 
    theta.save[k] = theta
    beta.save[,k] = beta 
    # D[k] = -2*spp.loglik(theta, beta, X, X.full, ds, n)
    
  };cat("\n")
  
  ###
  ###  Write Output
  ###
  
  list(beta.save=beta.save,beta.0.save=log(theta.save),n.mcmc=n.mcmc)
  # beta.0.save = log(theta.save)/beta.0.save
  
}