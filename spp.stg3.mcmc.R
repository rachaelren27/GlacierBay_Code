spp.stg3.mcmc <- function(out){
  
  #
  #  Second stage algorithm for fitting SPP model w/ unknown n before observed
  #
  
  ###
  ###  Setup Variables
  ###
  
  n=out$n
  n.mcmc=out$n.mcmc
  X.full=out$X.full
  p=dim(X.full)[2]
  ds=out$ds
  lam.int.save <- out$lam.int.save
  accept <- 0
  
  #gamma prior hyperparameters
  a <- 0.01
  b <- 0.01
  
  beta.0.save=rep(0,n.mcmc)
  beta.save=matrix(0,p,n.mcmc)
  
  # starting values
  beta=c(out$beta.save[,n.mcmc])
  lam.int <- lam.int.save[n.mcmc]
  theta <- rgamma(1, a + n, b + lam.int)
  
  for(k in 1:n.mcmc){
    if(k%%1000==0){cat(k," ")}
    
    ###
    ###  Propose theta
    ###
    theta.star <- rgamma(1, a + n, b + lam.int.save[k])
    
    ###
    ###  Propose beta
    ###
    
    idx.star=sample(1:n.mcmc,1)
    beta.star=c(out$beta.save[,idx.star])
    lam.int.star <- lam.int.save[idx.star]
    theta.lam.int.star=theta.star*lam.int.star 
    
    theta.lam.int <- theta*lam.int
      
    mh.1=dpois(n,theta.lam.int.star,log=TRUE)
    mh.2=dpois(n,theta.lam.int,log=TRUE)
    
    if(exp(mh.1-mh.2)>runif(1)){
      theta=theta.star
      beta=beta.star
      lam.int=lam.int.star
      accept = accept + 1
    }
    
    ###
    ###  Save Samples
    ###
    
    beta.0.save[k]=log(theta)
    beta.save[,k]=beta
    
  };cat("\n")
  
  print(accept/n.mcmc)
  
  ###
  ###  Write Output 
  ###
  
  list(beta.save=beta.save,beta.0.save=beta.0.save,n.mcmc=n.mcmc)
  
}
