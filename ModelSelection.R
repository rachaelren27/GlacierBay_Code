setwd("/Users/rlr3795/Desktop/GlacierBay_Project")
library(parallel)
library(doParallel)
library(foreach)

calc_DIC <- function(beta.post,X.obs,X.full,ds,n){
  spp.loglik <- function(theta,beta,X.obs,X.full,ds,n){ # pass in X.full, ds
    llam=log(theta)+X.obs%*%beta
    lam.int=sum(theta*exp(log(ds)+X.full%*%beta))
    sum(llam)-lam.int
  }
  
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)

  D <- foreach(k = 1:nrow(beta.post), .combine = c) %dopar% {
    D <- -2*spp.loglik(exp(beta.post[k,1]),beta.post[k,-1],X.obs,X.full,ds,n)

    return(D)
  }

  stopCluster(cl)
  
  theta.mean <- exp(mean(beta.post[,1]))
  beta.post <- as.matrix(beta.post[,-1])
  beta.mean <- apply(beta.post, 2, mean)
  D.hat <- -2*spp.loglik(theta.mean, beta.mean, X.obs, X.full, ds, n)
  D.bar <- mean(D)
  p_D <- D.bar - D.hat
  
  list(DIC = 2*D.bar - D.hat, p_D = p_D)
}

source(here("GlacierBay_Code", "spp_win_2D", "spp.comp.mcmc.R"))

n.mcmc <- 100000
n.burn <- 0.1*n.mcmc

n.models <- 6
X.win.full.list <- list()
X.obs.list <- list()

# set design matrices
X.win.full.list[[1]] <- as.matrix(X.win.full[,1]) # intercept + ice
X.obs.list[[1]] <- as.matrix(X.obs[,1])

X.win.full.list[[2]] <- as.matrix(X.win.full[,2]) # intercept + bath
X.obs.list[[2]] <- as.matrix(X.obs[,2])

X.win.full.list[[3]] <- as.matrix(X.win.full[,3]) # intercept + glacier
X.obs.list[[3]] <- as.matrix(X.obs[,3])

X.win.full.list[[4]] <- X.win.full[,1:2] # intercept + ice + bath
X.obs.list[[4]] <- X.obs[,1:2]

X.win.full.list[[5]] <- X.win.full[,c(1,3)] # intercept + ice + glacier
X.obs.list[[5]] <- X.obs[,c(1,3)]

X.win.full.list[[6]] <- X.win.full[,2:3] # intercept + bath + glacier
X.obs.list[[6]] <- X.obs[,2:3]

# fit models
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

beta.post.list <- foreach(k = 1:n.models) %dopar% {
  spp.comp.mcmc(seal.mat, X.obs.list[[k]], X.win.full.list[[k]],
                                  ds, n.mcmc, 0.1, 0.1)
}

stopCluster(cl)

# calculate DIC
DIC <- list()
for(i in 1:n.models){
  beta.out <- beta.post.list[[i]]
  beta.post <- t(rbind(beta.out$beta.0.save, beta.out$beta.save))[-(1:n.burn),]
  D <- beta.out$D
  
  DIC[[i]] <- calc_DIC(beta.post, X.obs.list[[i]], X.win.full.list[[i]], ds, n)
}

# --- Model 7: intercept + ice + bathymetry + glacier distance -----------------
out.comp.full=spp.comp.mcmc(seal.mat,X.obs,X.win.full,ds,n.mcmc,0.1,0.1)

# discard burn-in
beta.save.full.lik <- out.comp.full$beta.save[,-(1:n.burn)]
beta.0.save.full.lik <- out.comp.full$beta.0.save[-(1:n.burn)]

beta.post <- t(rbind(beta.0.save.full.lik, beta.save.full.lik))

# trace plots
layout(matrix(1:2,2,1))
plot(beta.0.save.full.lik,type="l")
matplot(t(beta.save.full.lik),lty=1,type="l", col = c("black", "red"))

m7.DIC <- calc_DIC(beta.post, X.obs, X.win.full, ds, n)


