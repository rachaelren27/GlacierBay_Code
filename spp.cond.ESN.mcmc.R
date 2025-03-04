spp.logit.bayesreg.ELM <- function(s.mat, X.obs.aug, X.full.aug, X.win.aug,
                                   win.idx, obs.idx, ds, n.mcmc, q){
  
  require(stats)
  require(VGAM)
  
  # subroutines
  
  spp.loglik <- function(W.obs.beta, W.win.full.beta, ds, n){
    lam.int=sum(exp(log(ds) + W.win.full.beta))
    sum(W.obs.beta) - n*log(lam.int)
  }
  
  gelu <- function(z){	
    z*pnorm(z)
  }
  
  # set up variables
  
  n <- nrow(s.mat)
  p <- dim(X.full)[2]
  
  theta.save <- rep(0, n.mcmc)
  beta.save <- matrix(0, q, n.mcmc)
  
  n.sim <- 100
  A.array <- array(0, c(p, q, n.sim))
  W.array <- array(0, c(nrow(X.obs.aug), q, n.sim))
  aic.vec <- rep(0, n.sim)
  beta.mat <- matrix(0, q+1, n.mcmc)
    
  for(l in 1:n.sim){
    A.array[,,l] <- matrix(rnorm(q*p), p, q)
    W.array[,,l] <- gelu(X.obs.aug%*%A.array[,,l])
    tmp.lm <- glm(y.obs.binary ~ W.array[,,l], family = binomial(link = "logit"))
    aic.vec[l] <- AIC(tmp.lm)
    beta.mat[,l] <- coef(tmp.lm)
  }
  best.idx <- (1:n.sim)[aic.vec == min(aic.vec)]
  cat("best AIC:", aic.vec[best.idx], "\n")
    
  A <- A.array[,, best.idx]
  W <- W.array[,, best.idx]
  beta0 <- beta.mat[1, best.idx]
  beta <- beta.mat[-1, best.idx]
  W.beta <- W%*%beta
  
  # set basis functions
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
    
    # sample beta
    
    # beta.star=t(mvnfast::rmvn(1,beta,beta.tune*diag(q)))
    # 
    # mh.1 = spp.loglik(theta,beta.star,W.obs.beta,W.win.full.beta,ds,n) +
    #        sum(mvnfast::dmvn(t(beta.star),mu.00,sigma.00,log=TRUE))
    # mh.2 = spp.loglik(theta,beta,W.obs.beta,W.win.full.beta,ds,n) +
    #        sum(mvnfast::dmvn(t(beta),mu.00,sigma.00,log=TRUE))
    
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