spp.logit.bayesreg.ELM <- function(X.obs.aug, X.full.aug, y.obs.binary, n, q, n.train){
  
  require(bayesreg)
  
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
    W.array[,,l] <- scale(gelu(X.obs.aug%*%A.array[,,l]), center = FALSE)
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
  out.bern.cond <- glm(y ~ W.obs, data = logit.obs.df,
                           family = binomial(link="logit"))
  
  beta.save <- out.bern.cond$beta
  
  beta.glm <- coef(out.bern.cond)[-1]
  vcov.glm <- vcov(out.bern.cond)[-1,-1]
  
  list(beta.glm=beta.glm, vcov.glm=vcov.glm, W.full=W.full, W.obs=W.obs, A=A)
  
}