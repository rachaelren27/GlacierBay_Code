spp.cond.mcmc <- function(s.mat,X,X.full,ds,n.mcmc){

#
#  1-D SPP while conditioning on n 
#

###
###  Subroutine 
###

spp.loglik <- function(beta,X,X.full,ds,n){
  llam=X%*%beta
  lam.int=sum(exp(log(ds)+X.full%*%beta))
  # list(loglik = sum(llam)-n*log(lam.int), lam.int = lam.int)
  sum(llam)-n*log(lam.int)
}

###
###  Set up variables
###

n=nrow(s.mat)
p=dim(X)[2]

beta.save=matrix(0,p,n.mcmc)
# beta.0.save=rep(0,n.mcmc) # these are MC draws from the prior
lam.int.save <- rep(0,n.mcmc)

###
###  Priors and Starting Values
###

mu.00=0
sig.00=10
mu.0=rep(0,p)
sig.0=rep(100,p)

beta=rep(0,p)
# spp.loglik.beta <- spp.loglik(beta,X,X.full,ds,n)

beta.tune=.1

###
###  MCMC Loop 
###

for(k in 1:n.mcmc){
  if(k%%1000==0){cat(k," ")}

  ###
  ###  Sample beta 
  ###
  

  beta.star=rnorm(p,beta,beta.tune)
  # spp.loglik.beta.star <- spp.loglik(beta.star,X,X.full,ds,n)
  
  mh.1=spp.loglik(beta.star,X,X.full,ds,n) + sum(dnorm(beta.star,mu.0,sig.0,log=TRUE))
  mh.2=spp.loglik(beta,X,X.full,ds,n) + sum(dnorm(beta,mu.0,sig.0,log=TRUE))
  if(exp(mh.1-mh.2)>runif(1)){
    beta=beta.star 
    # spp.loglik.beta <- spp.loglik.beta.star
  }

  ###
  ###  Save Samples
  ###
  
  beta.0.save[k]=rnorm(1,mu.00,sig.00)
  beta.save[,k]=beta 
  # lam.int.save[k] <- spp.loglik.beta$lam.int

};cat("\n")

###
###  Write Output
###

list(beta.save=beta.save,n.mcmc=n.mcmc,n=n,ds=ds,X.full=X.full)
# omit beta.0.save, add lam.int.save
}
