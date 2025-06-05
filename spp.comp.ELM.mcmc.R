spp.comp.ELM.mcmc <- function(s.mat, X.full, full.win.idx, obs.win.idx, ds, n.mcmc,
                              theta.tune, beta.tune, q, lambda){
  
  require(stats)
  require(VGAM)
  
  # subroutines
  
  spp.loglik <- function(theta, W.obs.beta, W.win.full.beta,ds){
    llam <- log(theta) + W.obs.beta
    lam.int <- sum(theta*exp(log(ds) + W.win.full.beta))
    sum(llam) - lam.int
  }
  
  gelu <- function(z){	
    z*pnorm(z)
  }
  
  # set up variables
  
  n <- nrow(s.mat)
  p <- dim(X.full)[2]
  
  theta.save <- rep(0, n.mcmc)
  beta.save <- matrix(0, q, n.mcmc)
  
  # set basis functions
  A <- matrix(rnorm(q*p), p, q)
   
  W.full <- gelu(X.full%*%A)
  W.obs <- W.full[obs.win.idx,]
  
  theta <- exp(10)
  beta <- rep(0, q)
  
  W.full.beta <- W.full%*%beta
  W.win.full.beta <- W.full.beta[full.win.idx,]
  W.obs.beta <- W.full.beta[obs.win.idx,]
  
  # priors/starting values
  
  # mu.00=0
  # b.00=1/10
  a <- 0.0000001
  b <- 0.0000001
  # mu.00=rep(0,q)
  # sigma.00=1*diag(q)
  mu.00=rep(0,q)
  b.00=1/(lambda*n)
  
  # mcmc
  
  for(k in 1:n.mcmc){
    if(k %% 1000 == 0){cat(k, " ")}
    
    beta.star <- t(mvnfast::rmvn(1,beta,beta.tune*diag(q)))
    W.full.beta.star <- W.full%*%beta.star
    W.win.full.beta.star <- W.full.beta.star[full.win.idx,]
    W.obs.beta.star <- W.full.beta.star[obs.win.idx,]
    
    mh.1 = spp.loglik(theta,W.obs.beta.star,W.win.full.beta.star,ds) +
           sum(dlaplace(beta.star,mu.00,b.00,log=TRUE))
    mh.2 = spp.loglik(theta,W.obs.beta,W.win.full.beta,ds) +
           sum(dlaplace(beta,mu.00,b.00,log=TRUE))
    
    if(exp(mh.1-mh.2)>runif(1)){
      beta <- beta.star
      W.full.beta <- W.full%*%beta
      W.win.full.beta <- W.full.beta[full.win.idx,]
      W.obs.beta <- W.full.beta[obs.win.idx,]
    }
    
    # sample theta
    
    theta.star <- rlnorm(1, log(theta), theta.tune)
    mh.1=spp.loglik(theta.star, W.obs.beta, W.win.full.beta,ds) + 
         dgamma(theta.star, a, b, log = TRUE) +
         dlnorm(theta, log(theta), theta.tune, log = TRUE)
    mh.2=spp.loglik(theta,W.obs.beta,W.win.full.beta,ds) + 
         dgamma(theta, a, b, log = TRUE) +
         dlnorm(theta.star, log(theta.star), theta.tune, log = TRUE)
    if(exp(mh.1-mh.2)>runif(1)){
      theta <- theta.star
    }
 
    # sample theta
    
    theta.save[k] <- theta
    beta.save[,k] <- beta 
    
  };cat("\n")
  
  list(beta.0.save=log(theta.save), beta.save=beta.save, n.mcmc=n.mcmc,
       W.full=W.full, W.obs=W.obs, A=A)

}