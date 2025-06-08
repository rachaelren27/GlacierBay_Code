spp.logit.ELM <- function(X.obs.aug, X.full.aug, y.obs.binary, q.vec){
  
  gelu <- function(z){	
    z*pnorm(z)
  }
  
  # pre-train ELM
  n.sim <- length(q.vec)
  p <- ncol(X.full.aug)
  A.save <- list()
  W.save <- list()
  beta.save <- list()
  out.glm <- list()
  aic.vec <- rep(0, n.sim)
    
  for(l in 1:n.sim){
    q <- q.vec[l]
    A.save[[l]] <- matrix(rnorm(q*p), p, q)
    W.save[[l]] <- scale(gelu(X.obs.aug%*%A.save[[l]]))
    out.glm[[l]] <- glm(y.obs.binary ~ W.save[[l]],
                        family = binomial(link = "logit"))
    aic.vec[l] <- AIC(out.glm[[l]])
  }
  best.idx <- (1:n.sim)[aic.vec == min(aic.vec)]
  cat("best AIC:", aic.vec[best.idx], "\n")
    
  A <- A.save[[best.idx]]
  W.obs <- W.save[[best.idx]]
  W.full <- scale(gelu(X.full.aug%*%A))
  
  out.glm <- glm(y.obs.binary ~ W.obs, family = binomial(link = "logit"))
  
  beta.glm <- coef(out.glm)[-1]
  vcov.glm <- vcov(out.glm)[-1,-1]
  
  list(beta.glm=beta.glm, vcov.glm=vcov.glm, W.full=W.full, W.obs=W.obs, A=A)
  
}