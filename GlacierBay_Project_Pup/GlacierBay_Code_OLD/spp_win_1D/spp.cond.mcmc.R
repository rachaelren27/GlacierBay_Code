spp.cond.mcmc <- function(s.mat,X,X.full,ds,n.mcmc){

#
#  1-D SPP while conditioning on n 
#

###
###  Subroutine 
###

spp.loglik <- function(beta){
  llam=X%*%beta
  lam.int=sum(exp(log(ds)+X.full%*%beta))
  sum(llam)-n*log(lam.int)
}

###
###  Set up variables
###

n=nrow(s.mat)
p=dim(X)[2]

beta.save=matrix(0,p,n.mcmc)
beta.0.save=rep(0,n.mcmc) # these are MC draws from the prior

###
###  Priors and Starting Values
###

mu.00=0
sig.00=10
mu.0=rep(0,p)
sig.0=rep(100,p)

beta=rep(0,p)

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
  mh.1=spp.loglik(beta.star)+sum(dnorm(beta.star,mu.0,sig.0,log=TRUE))
  mh.2=spp.loglik(beta)+sum(dnorm(beta,mu.0,sig.0,log=TRUE))
  if(exp(mh.1-mh.2)>runif(1)){
    beta=beta.star 
  }

  ###
  ###  Save Samples
  ###

  beta.0.save[k]=rnorm(1,mu.00,sig.00)
  beta.save[,k]=beta 

};cat("\n")

###
###  Write Output
###

list(beta.save=beta.save,beta.0.save=beta.0.save,n.mcmc=n.mcmc,n=n,ds=ds,X.full=X.full)

}
