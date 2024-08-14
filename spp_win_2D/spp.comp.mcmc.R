spp.comp.mcmc <- function(s.mat,X,X.full,ds,area,n.mcmc){

#
#  1-D SPP w/ complete data (s and n) likelihood and windowing
#

###
###  Subroutine 
###

spp.loglik <- function(beta.0,beta){ # pass in X.full, ds
  llam=beta.0+X%*%beta
  lam.int=sum(exp(log(ds)+beta.0+X.full%*%beta))
  sum(llam)-lam.int-lfactorial(n) # check if we can delete lfactorial(n)
}

###
###  Set up variables
###

n=nrow(s.mat)
p=dim(X)[2]

beta.save=matrix(0,p,n.mcmc)
beta.0.save=rep(0,n.mcmc)

###
###  Priors and Starting Values
###

mu.00=0
sig.00=10
mu.0=rep(0,p)
sig.0=rep(100,p)

beta.0=log(n)-log(area)
beta=rep(0,p)

beta.0.tune=.1
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
  mh.1=spp.loglik(beta.0,beta.star)+sum(dnorm(beta.star,mu.0,sig.0,log=TRUE))
  mh.2=spp.loglik(beta.0,beta)+sum(dnorm(beta,mu.0,sig.0,log=TRUE))
  if(exp(mh.1-mh.2)>runif(1)){
    beta=beta.star 
  }

  ###
  ###  Sample beta.0
  ###

  beta.0.star=rnorm(1,beta.0,beta.0.tune)
  mh.1=spp.loglik(beta.0.star,beta)+dnorm(beta.0.star,mu.00,sig.00,log=TRUE)
  mh.2=spp.loglik(beta.0,beta)+dnorm(beta.0,mu.00,sig.00,log=TRUE)
  if(exp(mh.1-mh.2)>runif(1)){
    beta.0=beta.0.star 
  }

  ###
  ###  Save Samples
  ###

  beta.0.save[k]=beta.0 
  beta.save[,k]=beta 

};cat("\n")

###
###  Write Output
###

list(beta.save=beta.save,beta.0.save=beta.0.save,n.mcmc=n.mcmc)

}
