spp.logit.bayesreg.ELM <- function(X.obs.aug, X.full.aug, y.obs.binary, n, q, n.mcmc,
                                   n.burn, n.train){
  
  require(stats)
  require(VGAM)
  
  gelu <- function(z){	
    z*pnorm(z)
  }
  
  # set up variables
  n <- n
  p <- dim(X.full)[2]
  
  # pre-train ELM
  n.train <- 100
  A.array <- array(0, c(p, q, n.train))
  W.array <- array(0, c(nrow(X.obs.aug), q, n.train))
  aic.vec <- rep(0, n.train)
  beta.mat <- matrix(0, q+1, n.train)
  
  for(l in 1:n.train){
    A.array[,,l] <- matrix(rnorm(q*p), p, q)
    W.array[,,l] <- gelu(X.obs.aug%*%A.array[,,l])
    tmp.lm <- glm(as.factor(y.obs.binary) ~ W.array[,,l], family = binomial(link = "logit"))
    aic.vec[l] <- AIC(tmp.lm)
    beta.mat[,l] <- coef(tmp.lm)
  }
  best.idx <- (1:n.train)[aic.vec == min(aic.vec)]
  cat("best AIC:", aic.vec[best.idx], "\n")
  
  A <- A.array[,, best.idx]
  W.obs <- W.array[,, best.idx]
  W.full <- gelu(X.full.aug%*%A)
  
  # fit logistic regression
  logit.obs.df <- data.frame(y = as.factor(y.obs.binary), W.obs)
  out.bern.cond <- bayesreg(y ~ X1 + X2 + X3 + X4 + X5, data = logit.obs.df,
                            model = "logistic",
                            n.samples = n.mcmc, burnin = n.burn)
  
  beta.save <- out.bern.cond$beta
  
  list(beta.save=beta.save, W.full=W.full, W.obs=W.obs, A=A)
  
}