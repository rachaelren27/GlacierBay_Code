setwd("/Users/rlr3795/Desktop/GlacierBay_Project")

library(tidyverse)
library(raster)
library(spatstat)
library(tictoc)
library(viridis)
# library(rstanarm)
# library(pgdraw)
library(mvnfast)
library(coda)
library(foreach)
library(doParallel)
# library(vioplot)

load(here("GlacierBay_Code","spp_win_2D", "sim2.RData"))
set.seed(1234)

# --- Simulate 2D data ---------------------------------------------------------
# set domain
x.domain <- c(0,1.05)
y.domain <- c(0,1.05)

# define the coordinates for window squares
win.length <- 0.2
gap <- 0.05
domain.length <- x.domain[2] - x.domain[1]
coords <- expand.grid(x = seq(gap, domain.length - win.length - gap, by = win.length + gap), 
                      y = seq(gap, domain.length - win.length - gap, by = win.length + gap))

# create individual squares
squares <- lapply(1:nrow(coords), function(i) {
  x0 <- coords$x[i]
  y0 <- coords$y[i]
  owin(xrange = c(x0, x0 + win.length), yrange = c(y0, y0 + win.length))
})

# combine squares into single window
combined.window <- do.call(union.owin, squares)

# calculate areas
tot.area <- (x.domain[2] - x.domain[1])*(y.domain[2] - y.domain[1])
tot.win.area <- area.owin(combined.window)
n.win <- 16
win.area <- tot.win.area/n.win
tot.nonwin.area <- tot.area - tot.win.area

# # plot the window
domain <- owin(xrange = c(0,domain.length), yrange = c(0,domain.length))

# plot(domain)
# plot(combined.window, add = TRUE)

# get quadrature grid
x.m <- 800
y.m <- 800
m <- x.m*y.m

x.full <- seq(x.domain[1], x.domain[2], length.out = x.m)
y.full <- x.full

s.full <- expand.grid(x = x.full, y = y.full)

# plot(domain)
# plot(combined.window, add = TRUE)
# points(s.full, pch = 19, cex = 0.05)

ds <- (x.full[2] - x.full[1])^2

# set X matrix
X.full <- matrix(0,m,2)
x1 <- s.full[,1]
x2 <- s.full[,2]
X.full[,1] <- x1
X.full[,2] <- x2

# set beta
beta <- c(2,1)
beta.0 <- 5
lam.full <- exp(beta.0+X.full%*%beta)
lam.max <- max(lam.full)

full.df <- as.data.frame(cbind(s.full, x1, x2, lam.full))
# # plot covariates and lambda
# ggplot(data = full.df, aes(x = x, y = y, col = x)) + 
#   geom_point(size = 0.5)
# ggplot(data = full.df, aes(x = x, y = y, col = y)) + 
#   geom_point(size = 0.5)
# 
# ggplot() +
#   geom_tile(data = full.df, aes(x = x, y = y, fill = lam.full)) + 
#   labs(fill = "lambda")
# 
# # create full raster
full.df <- full.df %>% rename(z = lam.full)
full.raster <- rasterFromXYZ(full.df)
# plot(full.raster, color = viridis(100))

# simulate observed points
M <- rpois(1, lam.max) 
x.superpop <- runif(M, x.domain[1], x.domain[2])
y.superpop <- runif(M, y.domain[1], y.domain[2])
s.superpop <- cbind(x.superpop, y.superpop)
X.superpop <- cbind(x.superpop, y.superpop)
lam.superpop <- exp(beta.0 + X.superpop%*%beta)

obs.idx <- rbinom(M, 1,lam.superpop/lam.max) == 1
s.obs <- s.superpop[obs.idx,] # total observed points 
X.obs <- X.superpop[obs.idx,] 
lam.obs <- lam.superpop[obs.idx]
N <- nrow(s.obs) 

# # plot superpop lambda
# superpop.df <- as.data.frame(cbind(x.superpop, y.superpop, lam.superpop))
# ggplot(data = superpop.df, aes(x = x.superpop, y = y.superpop, col = lam.superpop)) + 
#   geom_point(size = 0.5)


# --- Get windowed data --------------------------------------------------------
obs.win <- inside.owin(s.obs[,1], s.obs[,2], combined.window)
obs.win.idx <- (1:N)[obs.win]
n <- length(obs.win.idx)

full.win <- inside.owin(s.full[,1], s.full[,2], combined.window)
full.win.idx <- (1:m)[full.win]

s.win=s.obs[obs.win.idx,]
X.win=X.obs[obs.win.idx,]

X.win.full=X.full[full.win.idx,]
X.nowin.full=X.full[-full.win.idx,]

# plot windowed data
obs.df <- as.data.frame(s.obs)
ggplot(data = obs.df, aes(x = x.superpop, y = y.superpop,
                          col = factor(obs.win))) +
  geom_point()

# plot(domain)
# plot(combined.window, add = TRUE)
# points(s.obs[,1], s.obs[,2], col = factor(obs.win), pch = 19, cex = 0.5)
# 
pdf("sim2_data.pdf")
domain.sf <- st_as_sfc(as.polygonal(domain)) %>% st_sf()
combined.window.sf <- st_as_sfc(as.polygonal(combined.window)) %>% st_sf()
ggplot() +
  geom_sf(data = domain.sf) +
  geom_sf(data = combined.window.sf) +
  geom_point(aes(x = s.win[,1], y = s.win[,2]), col = "red", size = 0.5) +
  # geom_point(aes(x = s.obs[-obs.win.idx,1], y = s.obs[-obs.win.idx,2]), size = 0.5) +
  theme(axis.title = element_blank())
dev.off()

# --- Fit SPP w/ complete likelihood -------------------------------------------
n.mcmc <- 100000
source(here("GlacierBay_Code", "spp_win_2D", "spp.comp.mcmc.R"))
tic()
out.comp.full=spp.comp.mcmc(s.win, X.win, X.win.full, ds, n.mcmc, 0.1, 0.01)
toc() # 385 sec

# discard burn-in
n.burn <- 0.1*n.mcmc
beta.0.save <- out.comp.full$beta.0.save[-(1:n.burn)]
beta.save <- out.comp.full$beta.save[,-(1:n.burn)]

effectiveSize(beta.0.save)
effectiveSize(beta.save[1,])
effectiveSize(beta.save[2,])

# trace plot
layout(matrix(1:2,2,1))
plot(beta.0.save,type="l")
abline(h=beta.0,col=rgb(0,1,0,.8),lty=2)
matplot(t(beta.save),lty=1,type="l")
abline(h=beta,col=rgb(0,1,0,.8),lty=2)

# posterior summary
beta.save.full <- t(rbind(beta.0.save, beta.save))
apply(beta.save.full,2,mean) 
apply(beta.save.full,2,sd) 
apply(beta.save.full,2,quantile,c(0.025,.975))

# posterior for N
N.comp.save <- rep(0, n.mcmc - n.burn)

for(k in 1:(n.mcmc - n.burn)){
  if(k%%10000==0){cat(k," ")}
  beta.0.tmp=beta.0.save[k]
  beta.tmp=beta.save[,k]
  lam.nowin.int=sum(exp(log(ds)+beta.0.tmp+X.nowin.full%*%beta.tmp))
  N.comp.save[k]=n+rpois(1,lam.nowin.int)
};cat("\n")

hist(N.comp.save,breaks=50,prob=TRUE,main="",xlab="N")
abline(v=N,col=rgb(0,1,0,.8),lty=2,lwd=2)


# --- Prepare Logistic Regression X matrix -------------------------------------
# sample background points
n.bg <- 50000
bg.pts <- rpoint(n.bg, win = combined.window)

n.mcmc <- 100000
n.burn <- 0.1*n.mcmc

# plot(domain)
# plot(combined.window, add = TRUE)
# points(bg.pts$x, bg.pts$y)

# prepare X matrix
X.bg <- cbind(bg.pts$x, bg.pts$y)
X.bern <- rbind(X.win, X.bg)
y.bern <- rep(0, n + n.bg)
y.bern[1:n] <- 1
bern.df <- data.frame(y = y.bern, x1 = X.bern[,1], x2 = X.bern[,2])

# --- Bayesian GLM (logistic) --------------------------------------------------
# stage 1
out.bayes.glm <- stan_glm(y ~ x1 + x2, family = binomial(link="logit"), 
                          data = bern.df, iter = n.mcmc, warmup = n.burn,
                          chains = 1) # 75 sec
beta.save <- t(as.matrix(out.bayes.glm)[,-1])

# stage 2
lam.int.save <- c()

cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

tic()
lam.int.save <- foreach(k = 1:ncol(beta.save), .combine = c) %dopar% {
  lam.int <- sum(exp(log(ds) + X.win.full %*% beta.save[,k]))
  return(lam.int) 
}
toc()

stopCluster(cl) # 16.7 sec

out.bayes.glm2 <- list(beta.save = beta.save, 
                    n.mcmc = n.mcmc - n.burn, n = n, ds = ds, X.full = X.win.full,
                    lam.int.save = lam.int.save)

# stage 3
source(here("GlacierBay_Code", "spp.stg3.mcmc.R"))
tic()
out.bayes.glm3 <- spp.stg3.mcmc(out.bayes.glm2)
toc()

beta.save <- out.bayes.glm3$beta.save
beta.0.save <- out.bayes.glm3$beta.0.save

effectiveSize(beta.0.save)
effectiveSize(beta.save[1,])
effectiveSize(beta.save[2,])

layout(matrix(1:2,2,1))
plot(beta.0.save,type="l")
abline(h=beta.0,col=rgb(0,1,0,.8),lty=2)
matplot(t(beta.save),lty=1,type="l")
abline(h=beta,col=rgb(0,1,0,.8),lty=2)


# --- Fit SPP w/ cond. likelihood (Polya-Gamma) --------------------------------
source(here("GlacierBay_Code", "Polya_Gamma.R"))
X.pg <- cbind(rep(1, nrow(X.bern)), X.bern) # add intercept
p <- ncol(X.pg)
mu.beta <- rep(0, p)
sigma.beta <- diag(100, p)

tic()
beta.save.pg <- polya_gamma(y.bern, X.pg,
                            mu.beta, sigma.beta, n.mcmc)
toc() # 1039 sec

beta.save <- beta.save.pg$beta[-1,]

# stage 2
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

tic()
lam.int.save <- foreach(k = 1:ncol(beta.save), .combine = c) %dopar% {
  lam.int <- sum(exp(log(ds) + X.win.full%*%beta.save[,k]))
  return(lam.int)
}
toc() # 12 sec

stopCluster(cl)

out.pg2 <- list(beta.save = beta.save,
                 n.mcmc = n.mcmc, n = n, ds = ds, X.full = X.win.full,
                 lam.int.save = lam.int.save)

# stage 3
source(here("GlacierBay_Code", "spp.stg3.mcmc.R"))
out.pg3 <- spp.stg3.mcmc(out.pg2)

beta.save <- out.pg3$beta.save
beta.0.save <- out.pg3$beta.0.save

effectiveSize(beta.0.save)
effectiveSize(beta.save[1,])
effectiveSize(beta.save[2,])

layout(matrix(1:2,2,1))
plot(beta.0.save,type="l")
abline(h=beta.0,col=rgb(0,1,0,.8),lty=2)
matplot(t(beta.save),lty=1,type="l")
abline(h=beta,col=rgb(0,1,0,.8),lty=2)


# --- Non-Bayesian GLM (exact, logistic) ---------------------------------------
## stage 1
out.glm <- glm(y ~ x1 + x2, family = binomial(link="logit"), data = bern.df)
beta.glm <- coef(out.glm)[-1]
vcov.glm <- vcov(out.glm)[-1,-1]
# sample from glm estimated density
beta.save <- t(mvnfast::rmvn(n.mcmc, mu = beta.glm, sigma = vcov.glm))

## stage 2
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

tic()
lam.int.save <- foreach(k = 1:ncol(beta.save), .combine = c) %dopar% {
  lam.int <- sum(exp(log(ds) + X.win.full%*%beta.save[,k]))
  return(lam.int)
}

X.beta.sum.save <- foreach(k = 1:ncol(beta.save)) %dopar% {
  X.beta.sum <- sum(X.win%*%beta.save[,k])
  return(X.beta.sum)
}
toc() # 31 sec

stopCluster(cl) 

out.glm2 <- list(beta.save = beta.save, mu.beta = beta.glm, sigma.beta = diag(se.glm^2),
                 n.mcmc = n.mcmc, n = n, ds = ds, X.full = X.win.full,
                 X.beta.sum.save = X.beta.sum.save, lam.int.save = lam.int.save)

## stage 3
source(here("GlacierBay_Code", "spp_win_2D", "spp.stg3.mcmc.nb.R"))
out.glm3 <- spp.stg3.mcmc.nb(out.glm2)

beta.save <- out.glm3$beta.save
beta.0.save <- out.glm3$beta.0.save

effectiveSize(beta.0.save)
effectiveSize(beta.save[1,])
effectiveSize(beta.save[2,])

layout(matrix(1:2,2,1))
plot(beta.0.save,type="l")
abline(h=beta.0,col=rgb(0,1,0,.8),lty=2)
matplot(t(beta.save),lty=1,type="l")
abline(h=beta,col=rgb(0,1,0,.8),lty=2)


# --- Non-Bayesian GLM (approx., logistic) -------------------------------------
# stage 1
out.glm <- glm(y ~ x1 + x2, family=binomial(link="logit"), data=bern.df)
beta.glm <- coef(out.glm)[-1]
se.glm <- summary(out.glm)$coefficients[-1, 2]
# sample from glm estimated density
beta.save <- t(mvnfast::rmvn(n.mcmc, mu = beta.glm, sigma = diag(se.glm^2)))

# stage 2
lam.int.save <- c()

tic()
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

lam.int.save <- foreach(k = 1:ncol(beta.save), .combine = c) %dopar% {
  lam.int <- sum(exp(log(ds) + X.win.full%*%beta.save[,k]))
  return(lam.int)
}

stopCluster(cl)
toc()

out.glm.approx2 <- list(beta.save = beta.save, 
                       n.mcmc = n.mcmc, n = n, ds = ds, X.full = X.win.full,
                       lam.int.save = lam.int.save)

# stage 3
source(here("GlacierBay_Code", "spp.stg3.mcmc.R"))
out.glm.approx3 <- spp.stg3.mcmc(out.glm.approx2)

beta.save <- out.glm.approx3$beta.save
beta.0.save <- out.glm.approx3$beta.0.save

effectiveSize(beta.0.save)
effectiveSize(beta.save[1,])
effectiveSize(beta.save[2,])

layout(matrix(1:2,2,1))
plot(beta.0.save,type="l")
abline(h=beta.0,col=rgb(0,1,0,.8),lty=2)
matplot(t(beta.save),lty=1,type="l")
abline(h=beta,col=rgb(0,1,0,.8),lty=2)

# ## compare Bayesian vs non-Bayesian stage 1 samples
# layout(matrix(1:2,1,2))
# plot(density(beta.save.bayes[,2], n = 1000), col = "red", 
#      lwd = 2, ylim = c(0,1))
# lines(density(beta.save[,1], n = 1000), col = "green",
#       lwd = 2)
# abline(v = beta[1], lwd = 2, lty = 2)
# plot(density(beta.save.bayes[,3], n = 1000), col = "red",
#      lwd = 2, ylim = c(0,1))
# lines(density(beta.save[,2], n=1000), col = "green",
#       lwd = 2)
# abline(v = beta[2], lwd = 2, lty = 2)
# 
# # compare joint posterior cross-sections
# plot(x = beta.save.bayes[,2], y = beta.save.bayes[,3],
#      xlab = TeX('$\\beta_1$'), ylab = TeX('$\\beta_2$'))
# points(x = beta.save[,1], y = beta.save[,2],
#        col = "red")


# --- Non-Bayesian GLM (exact, Poisson) ----------------------------------------
obs.cell <- cellFromXY(full.raster, s.win)
row.counts <- table(factor(obs.cell, levels = 1:nrow(s.full)))
X.pois <- cbind(X.full, row.counts)[full.win.idx,]
pois.df <- data.frame(x1 = X.pois[,1], x2 = -X.pois[,2], y = X.pois[,3])

out.pois.glm <- glm(y ~ x1 + x2, data = pois.df, family = poisson(link = "log"))
beta.glm <- coef(out.pois.glm)[-1]
se.glm <- summary(out.pois.glm)$coefficients[-1, 2]
# sample from glm estimated density
beta.save <- t(mvnfast::rmvn(n.mcmc, mu = beta.glm, sigma = diag(se.glm^2)))

## stage 2
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

tic()
lam.int.save <- foreach(k = 1:ncol(beta.save), .combine = c) %dopar% {
  lam.int <- sum(exp(log(ds) + X.win.full%*%beta.save[,k]))
  return(lam.int)
}

X.beta.sum.save <- foreach(k = 1:ncol(beta.save)) %dopar% {
  X.beta.sum <- sum(X.win%*%beta.save[,k])
  return(X.beta.sum)
}
toc() # 31 sec

stopCluster(cl) 

out.pois.glm2 <- list(beta.save = beta.save, mu.beta = beta.glm, sigma.beta = diag(se.glm^2),
                 n.mcmc = n.mcmc, n = n, ds = ds, X.full = X.win.full,
                 X.beta.sum.save = X.beta.sum.save, lam.int.save = lam.int.save)

## stage 3
source(here("GlacierBay_Code", "spp_win_2D", "spp.stg3.mcmc.nb.R"))
out.pois.glm3 <- spp.stg3.mcmc.nb(out.pois.glm2)

beta.save <- out.pois.glm3$beta.save
beta.0.save <- out.pois.glm3$beta.0.save

effectiveSize(beta.0.save)
effectiveSize(beta.save[1,])
effectiveSize(beta.save[2,])

layout(matrix(1:2,2,1))
plot(beta.0.save,type="l")
abline(h=beta.0,col=rgb(0,1,0,.8),lty=2)
matplot(t(beta.save),lty=1,type="l")
abline(h=beta,col=rgb(0,1,0,.8),lty=2)


# --- Beta Posterior Summary ---------------------------------------------------
# posterior summary
beta.save.full <- t(rbind(beta.0.save, beta.save))
apply(beta.save.full,2,mean) 
apply(beta.save.full,2,sd) 
apply(beta.save.full,2,quantile,c(0.025,.975))


# --- N Posterior Predictive ---------------------------------------------------
# posterior for N
N.comp.save <- rep(0, n.mcmc - n.burn)

for(k in 1:(n.mcmc - n.burn)){
  if(k%%10000 == 0){cat(k," ")}
  beta.0.tmp <- beta.0.save[k]
  beta.tmp <- beta.save[,k]
  lam.nowin.int <- sum(exp(log(ds) + beta.0.tmp + X.nowin.full%*%beta.tmp))
  N.comp.save[k] <- n + rpois(1, lam.nowin.int)
};cat("\n")

mean(N.comp.save)
sd(N.comp.save)

hist(N.comp.save,breaks=50,prob=TRUE,main="",xlab="N")
abline(v=N,col=rgb(0,1,0,.8),lty=2,lwd=2)

# --- Compare marginal posteriors ----------------------------------------------

layout(matrix(1:3,1,3))
# hist(out.comp.full$beta.0.save[-(1:n.burn)],prob=TRUE,breaks=60,main="",xlab=bquote(beta[0]),
#      ylim = c(0,1), cex.axis=1.5, cex.lab=1.3)
plot(density(out.comp.full$beta.0.save[-(1:n.burn)], n = 1000), col = "red", 
     lwd = 2, ylim = c(0,10))
lines(density(out.glm3$beta.0.save[-(1:n.burn)], n = 1000), col = "green", 
      lwd = 2)
lines(density(out.glm.approx3$beta.0.save[-(1:n.burn)], n = 1000), col = "blue",
      lwd = 2)
lines(density(out.bayes.glm3$beta.0.save[-(1:n.burn)], n = 1000), col = "orange",
      lwd = 2)
# lines(density(out.pois.glm3$beta.0.save[-(1:n.burn)], n = 1000), col = "purple", 
#       lwd = 2)
lines(density(out.pg3$beta.0.save[-(1:n.burn)], n = 1000), col = "pink", 
      lwd = 2)
abline(v = beta.0, lwd = 2, lty = 2)
# hist(out.comp.full$beta.save[1,-(1:n.burn)],prob=TRUE,breaks=60,main="",xlab=bquote(beta[1]),
#      ylim = c(0,1), cex.axis=1.5, cex.lab=1.3)
plot(density(out.comp.full$beta.save[1,-(1:n.burn)], n = 1000), col = "red", 
     lwd = 2, ylim = c(0,10))
lines(density(out.glm3$beta.save[1,-(1:n.burn)], n = 1000), col = "green",
      lwd = 2)
lines(density(out.pg3$beta.save[1,-(1:n.burn)], n = 1000), col = "pink",
      lwd = 2)
lines(density(out.glm.approx3$beta.save[1,-(1:n.burn)], n = 1000), col = "blue",
      lwd = 2)
lines(density(out.bayes.glm3$beta.save[1,-(1:n.burn)], n = 1000), col = "orange",
      lwd = 2)
# lines(density(out.pois.glm3$beta.save[1,-(1:n.burn)], n = 1000), col = "purple",
#       lwd = 2)
abline(v = beta[1], lwd = 2, lty = 2)
# hist(out.comp.full$beta.save[2,-(1:n.burn)],prob=TRUE,breaks=60,main="",xlab=bquote(beta[2]),
#      ylim = c(0,1), cex.axis=1.5, cex.lab=1.3)
plot(density(out.comp.full$beta.save[2,-(1:n.burn)], n = 1000), col = "red",
     lwd = 2, ylim = c(0,10))
lines(density(out.glm3$beta.save[2, -(1:n.burn)], n=1000), col = "green",
      lwd = 2)
lines(density(out.glm.approx3$beta.save[2, -(1:n.burn)], n=1000), col = "blue",
      lwd = 2)
lines(density(out.bayes.glm3$beta.save[2,-(1:n.burn)], n = 1000), col = "orange",
      lwd = 2)
# lines(density(out.pois.glm3$beta.save[2,-(1:n.burn)], n = 1000), col = "purple",
#       lwd = 2)
lines(density(out.pg3$beta.save[2,-(1:n.burn)], n = 1000), col = "pink",
      lwd = 2)
abline(v = beta[2], lwd = 2, lty = 2)
# dev.off()

n_keep <- n.mcmc - n.burn
df <- data.frame(
  value = c(
    # beta[0]
    out.comp.full$beta.0.save[-(1:n.burn)],
    out.pg3$beta.0.save[-(1:n.burn)],
    out.bayes.glm3$beta.0.save,
    out.glm3$beta.0.save[-(1:n.burn)],
    out.glm.approx3$beta.0.save[-(1:n.burn)],
    
    # beta[1]
    out.comp.full$beta.save[1, -(1:n.burn)],
    out.pg3$beta.save[1, -(1:n.burn)],
    out.bayes.glm3$beta.save[1,],
    out.glm3$beta.save[1, -(1:n.burn)],
    out.glm.approx3$beta.save[1, -(1:n.burn)],
    
    # beta[2]
    out.comp.full$beta.save[2, -(1:n.burn)],
    out.pg3$beta.save[2, -(1:n.burn)],
    out.bayes.glm3$beta.save[2,],
    out.glm3$beta.save[2, -(1:n.burn)],
    out.glm.approx3$beta.save[2, -(1:n.burn)]
  ),
  
  method = rep(
    c("complete", "PG", "stan_glm", "NB-GLM-E", "NB-GLM-A"),
    each = n_keep, times = 3
  ),
  
  coefficient = rep(
    c("beta[0]", "beta[1]", "beta[2]"),
    each = 5 * n_keep
  )
)

df$method <- factor(df$method, levels = c("complete", "PG", "stan_glm", "NB-GLM-E", "NB-GLM-A"))
df$coefficient <- factor(df$coefficient, levels = c("beta[0]", "beta[1]", "beta[2]"))

pdf("compare_marginals.pdf")
ggplot(df, aes(x = coefficient, y = value, fill = method)) +
  geom_violin(position = position_dodge(0.95), trim = FALSE, color = "black") +
  scale_x_discrete(labels = c(TeX("$\\beta_0$"), TeX("$\\beta_1$"), 
                              TeX("$\\beta_2$"))) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Coefficient", y = "Value", fill = "Method"
  )
dev.off()