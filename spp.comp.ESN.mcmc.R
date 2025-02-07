spp.comp.ESN.mcmc <- function(s.mat,X.full,full.win.idx,obs.win.idx,q,ds,n.mcmc,beta.tune){
  
  require(stats)
  
  #
  #  1-D SPP w/ complete data (s and n) likelihood and windowing
  #
  
  ###
  ###  Subroutine 
  ###
  
  spp.loglik <- function(beta,W.obs.beta,W.win.full.beta,ds,n){ # pass in X.full, ds
    llam=W.obs.beta
    lam.int=sum(exp(log(ds)+W.win.full.beta))
    sum(llam)-lam.int-lfactorial(n) # check if we can delete lfactorial(n)
  }
  
  gelu <- function(z){	
    z*pnorm(z)
  }
  
  ###
  ###  Priors and Starting Values
  ###
  
  mu.00=0
  sig.00=10
  mu.0=rep(0,q+1)
  sig.0=rep(100,q+1)
  
  ###
  ###  Set up variables
  ###
  
  n=nrow(s.mat)
  p=dim(X.full)[2]
  X.win.full=X.full[full.win.idx,]
  X.obs=X.full[obs.win.idx,]
  
  beta.save=matrix(0,q+1,n.mcmc)
  
  A=matrix(rnorm(q*p),p,q)
  
  W.full=gelu(X.full%*%A)
  W.full=prcomp(W.full)$x
  W.obs=W.full[obs.win.idx,]
  
  # add intercept
  q=q+1
  beta=rep(0,q)
  W.full=cbind(rep(1,nrow(W.full)),W.full)
  
  W.full.beta=W.full%*%beta
  W.win.full.beta=W.full.beta[full.win.idx,]
  W.obs.beta=W.full.beta[obs.win.idx,]
  
  ###
  ###  MCMC Loop 
  ###
  
  for(k in 1:n.mcmc){
    if(k%%1000==0){cat(k," ")}
    
    ###
    ###  Sample beta 
    ###z
    
    beta.star=rnorm(q,beta,beta.tune)
    mh.1=spp.loglik(beta.star,W.obs.beta,W.win.full.beta,ds,n)+sum(dnorm(beta.star,mu.0,sig.0,log=TRUE))
    mh.2=spp.loglik(beta,W.obs.beta,W.win.full.beta,ds,n)+sum(dnorm(beta,mu.0,sig.0,log=TRUE))
    if(exp(mh.1-mh.2)>runif(1)){
      beta=beta.star
    }
 
    ###
    ###  Save Samples
    ###
    
    beta.save[,k] = beta 
    W.full.beta = W.full%*%beta
    W.win.full.beta = W.full.beta[full.win.idx,]
    W.obs.beta = W.full.beta[obs.win.idx,]
    
  };cat("\n")
  
  ###
  ###  Write Output
  ###
  
  list(beta.save=beta.save,n.mcmc=n.mcmc,W.full=W.full,W.obs=W.obs)

}