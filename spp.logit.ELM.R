spp.logit.ELM <- function(X.obs.aug, X.full, y.obs.binary, q.vec){
  
  # sigmoid <- function(x){	
  #   1/(1 + exp(-x))
  # }
  
  gelu <- function(x){
    x*pnorm(x)
  }
  
  ## pre-train ELM
  n.sim <- length(q.vec)
  p <- ncol(X.full.aug)
  A.save <- list()
  W.save <- list()
  glm.list <- list()
  aic.vec <- c()
    
  for(l in 1:n.sim){
    q <- q.vec[l]
    A.save[[l]] <- matrix(rnorm(q*p), p, q)
    W.save[[l]] <- gelu(X.obs.aug%*%A.save[[l]])
    glm.list[[l]] <- glm(y.obs.binary ~ W.save[[l]],
                        family = binomial(link = "logit"))
    aic.vec[l] <- AIC(glm.list[[l]])
  }
  best.idx <- (1:n.sim)[aic.vec == min(aic.vec)]
  cat("best AIC:", aic.vec[best.idx], "\n")
    
  A <- A.save[[best.idx]]
  W.obs <- W.save[[best.idx]]
  W.full <- gelu(X.full%*%A)
  out.glm <- glm.list[[best.idx]]
  
  beta.glm <- coef(out.glm)[-1]
  vcov.glm <- vcov(out.glm)[-1,-1]
  
  list(beta.glm=beta.glm, vcov.glm=vcov.glm, W.full=W.full, W.obs=W.obs, A=A)
  
}