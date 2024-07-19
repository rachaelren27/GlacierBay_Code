spp.stg2.mcmc <- function(out){

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

beta.0.save=rep(0,n.mcmc)
beta.save=matrix(0,p,n.mcmc)

###
###  Calculate k=0 
###

beta.0=out$beta.0.save[n.mcmc]
beta=c(out$beta.save[,n.mcmc])
lam.int=sum(exp(log(ds)+beta.0+X.full%*%beta))

###
###  MCMC Loop 
###

for(k in 1:n.mcmc){
  if(k%%1000==0){cat(k," ")}

  ###
  ###  Update beta.0 and beta 
  ###
  
  idx.star=sample(1:n.mcmc,1)
  beta.0.star=out$beta.0.save[idx.star]
  beta.star=c(out$beta.save[,idx.star])
  lam.int.star=sum(exp(log(ds)+beta.0.star+X.full%*%beta.star)) # can be done before in parallel

  mh.1=dpois(n,lam.int.star,log=TRUE)
  mh.2=dpois(n,lam.int,log=TRUE)

  if(exp(mh.1-mh.2)>runif(1)){
    beta.0=beta.0.star
    beta=beta.star
    lam.int=lam.int.star
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
