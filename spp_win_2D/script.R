setwd("/Users/rlr3795/Desktop/GlacierBay_Project")

library(tidyverse)
library(raster)
library(spatstat)
library(tictoc)
library(here)

# --- Simulate 2D data ---------------------------------------------------------
set.seed(1234)

# set domain
x.domain <- c(0,0.9)
y.domain <- c(0,0.9)

# define the coordinates for window squares
win.length <- 0.1
gap <- 0.1
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

# plot the window
domain <- owin(xrange = c(0,domain.length), yrange = c(0,domain.length))

plot(domain)
plot(combined.window, add = TRUE)

# n.win <- 2
# tot.area <- (x.domain[2] - x.domain[1])*(y.domain[2] - y.domain[1])
# tot.win.area <- tot.area/2
# win.area <- tot.win.area/n.win
# tot.nonwin.area <- tot.area - tot.win.area
# nonwin.area <- tot.nonwin.area/n.win

# get quadrature grid
x.m <- 100
y.m <- 100
m <- x.m*y.m

x.full <- seq(x.domain[1], x.domain[2], length.out = x.m)
y.full <- x.full

s.full <- expand.grid(x = x.full, y = y.full)

plot(domain)
plot(combined.window, add = TRUE)
points(s.full, pch = 19, cex = 0.05)

ds <- (x.full[2] - x.full[1])^2

# set X matrix
X.full <- matrix(0,m,2)
x1 <- s.full[,1]
x2 <- s.full[,2]
X.full[,1] <- x1
X.full[,2] <- x2

# set beta
beta <- c(2,1)
beta.0 <- 4
lam.full <- exp(beta.0+X.full%*%beta)
lam.max <- max(lam.full)

full.df <- as.data.frame(cbind(s.full, x1, x2, lam.full))
# plot covariates and lambda
ggplot(data = full.df, aes(x = x, y = y, col = x)) + 
  geom_point(size = 0.5)
ggplot(data = full.df, aes(x = x, y = y, col = y)) + 
  geom_point(size = 0.5)

ggplot(data = full.df, aes(x = x, y = y, col = lam.full)) + 
  geom_point(size = 0.5)

# create full raster
full.df <- full.df %>% rename(z = lam.full)
full.raster <- rasterFromXYZ(full.df)
plot(full.raster)

# simulate observed points
M=rpois(1, lam.max) 
x.superpop <- runif(M, x.domain[1], x.domain[2])
y.superpop <- runif(M, y.domain[1], y.domain[2])
s.superpop <- cbind(x.superpop, y.superpop)
X.superpop <- cbind(x.superpop, y.superpop)
lam.superpop=exp(beta.0+X.superpop%*%beta)

obs.idx=rbinom(M,1,lam.superpop/lam.max)==1
s.obs=s.superpop[obs.idx,] # total observed points 
X.obs=X.superpop[obs.idx,] 
lam.obs <- lam.superpop[obs.idx]
N=nrow(s.obs) # 218

# plot superpop lambda
superpop.df <- as.data.frame(cbind(x.superpop, y.superpop, lam.superpop))
ggplot(data = superpop.df, aes(x = x.superpop, y = y.superpop, col = lam.superpop)) + 
  geom_point(size = 0.5)

# --- Get windowed data --------------------------------------------------------
obs.win <- inside.owin(s.obs[,1], s.obs[,2], combined.window)
obs.win.idx <- (1:N)[obs.win]
n=length(obs.win.idx) # 44

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

plot(domain)
plot(combined.window, add = TRUE)
points(s.obs[,1], s.obs[,2], col = factor(obs.win), pch = 19, cex = 0.5)

# --- Fit SPP w/ complete likelihood -------------------------------------------
n.mcmc=100000
source(here("GlacierBay_Code", "spp_win_2D", "spp.comp.mcmc.R"))
tic()
out.comp.full=spp.comp.mcmc(s.win,X.win,X.win.full,ds,win.area,n.mcmc,0.1,0.1)
toc() # 6.568 sec

# discard burn-in
n.burn <- 0.1*n.mcmc
beta.0.save <- out.comp.full$beta.0.save[-(1:n.burn)]
beta.save <- out.comp.full$beta.save[,-(1:n.burn)]

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

# --- Fit SPP w/ conditional likelihood ----------------------------------------
n.mcmc=100000
source(here("GlacierBay_Code", "spp_win_2D", "spp.cond.mcmc.R"))
out.cond.full <- spp.cond.mcmc(s.win,X.win,X.win.full,ds,n.mcmc)

layout(matrix(1:2,2,1))
plot(out.cond.full$beta.0.save,type="l")
abline(h=beta.0,col=rgb(0,1,0,.8),lty=2)
matplot(t(out.cond.full$beta.save),lty=1,type="l")
abline(h=beta,col=rgb(0,1,0,.8),lty=2)

effectiveSize(out.cond.full$beta.save[1,]) # 407
effectiveSize(out.cond.full$beta.save[2,]) # 463

# --- Fit SPP w/ cond. likelihood Bernoulli GLM --------------------------------
# sample background points
n.bg <- 10000
bg.pts <- rpoint(n.bg, win = combined.window)

plot(domain)
plot(combined.window, add = TRUE)
points(bg.pts$x, bg.pts$y)

# prepare X matrix
X.bg <- cbind(bg.pts$x, bg.pts$y)
X.bern <- rbind(X.win, X.bg)
y.bern <- rep(0, n + n.bg)
y.bern[1:n] <- 1

bern.rsf.df <- data.frame(y = y.bern, x1 = X.bern[,1], x2 = X.bern[,2])
out.bern.cond <- stan_glm(y ~ x1 + x2, family=binomial(link="logit"), data=bern.rsf.df,
                          iter = 100000, chains = 1)

# --- Fit SPP w/ cond. likelihood (Polya-Gamma 1st stage) ----------------------
source(here("GlacierBay_Code", "Polya_Gamma.R"))
X.pg <- cbind(rep(1, nrow(X.bern)), X.bern)
p <- ncol(X.pg)
mu.beta <- rep(0, p)
sigma.beta <- diag(2.25, p)
w <- c()
for(i in 1:length(y.bern)){
  w[i] <- ifelse(y.bern[i] == 0, 10000, 1)
}
tic()
beta.save.pg <- polya_gamma(y.bern, X.pg,
                            mu.beta, sigma.beta, 100000)
toc() # 236 sec

beta.0.save <- beta.save.pg$beta[1,]
beta.save <- beta.save.pg$beta[2:3,]

# # check weights
# y_i.hat <- c()
# for(i in 1:n.mcmc){
#   y_i.hat[i] <- ((tot.win.area*exp(beta.0.save[i] + t(X.pg[i,])%*%beta.save[,i]))/(100000*1000))/
#                 ((1 + tot.win.area*exp(beta.0.save[i] + t(X.pg[i,])%*%beta.save[,i]))/(100000*1000))
# }

# plot(beta.save.pg$beta[1,], type = "l")

plot(beta.save.pg$beta[2,], type = "l")
abline(h=2,col=rgb(0,1,0,.8),lty=2)
effectiveSize(beta.save.pg$beta[2,]) # 2111

plot(beta.save.pg$beta[3,], type = "l")
abline(h=1,col=rgb(0,1,0,.8),lty=2)
effectiveSize(beta.save.pg$beta[3,]) # 2326

matplot(t(beta.save.pg$beta[-1,]),lty=1,type="l")
abline(h=beta,col=rgb(0,1,0,.8),lty=2)

# prepare for second stage
# beta.0.precise <- rnorm(n.mcmc, mean = 3.7, sd = 1)
out.cond.pg <- list(beta.save = beta.save.pg$beta[-1,], beta.0.save = rnorm(n.mcmc, 0, 10),
                    n.mcmc = 100000, n = n, ds = ds, X.full = X.win.full)

# --- 2nd stage: compute lambda integrals --------------------------------------
beta.save <- beta.save.pg$beta[2:3,]
# theta.save <- rep(0,n.mcmc)
lam.int.save <- c()

for(k in 1:n.mcmc){
  if(k%%1000==0){cat(k," ")}
  lam.int.save[k] <- sum(exp(log(ds)+X.win.full%*%beta.save[,k]))
  # theta.save[k] <- rgamma(1, 0.01 + n, rate = 0.01 + lam.int)
};cat("\n")

# beta.0.save <- log(theta.save)

# plot(beta.0.save, type ="l")

# # compare marginal posterior
# hist(out.comp.full$beta.0.save[-(1:1000)],prob=TRUE,breaks=60,main="",xlab=bquote(beta[0]),
#      ylim = c(0,1))
# lines(density(beta.0.save[-(1:1000)],n=1000, adjust = 2),col="red",lwd=2)
# lines(density(out.cond.pg2$beta.0.save[-(1:1000)],n=1000, adjust = 2),
#       col="green",lwd=2)

out.cond.pg <- list(beta.save = beta.save.pg$beta[-1,], 
                    # beta.0.save = log(rgamma(n.mcmc, n + 0.001, 1.001)),
                    n.mcmc = n.mcmc, n = n, ds = ds, X.full = X.win.full,
                    lam.int.save = lam.int.save)

# --- 3rd stage MCMC -----------------------------------------------------------
source(here("GlacierBay_Code", "spp.stg3.mcmc.R"))
# tic()
out.cond.pg3 <- spp.stg3.mcmc(out.cond.pg)
# toc()

# discard burn-in
beta.save <- out.cond.pg3$beta.save[,-(1:n.burn)]
beta.0.save <- out.cond.pg3$beta.0.save[-(1:n.burn)]

layout(matrix(1:2,2,1))
plot(beta.0.save,type="l")
abline(h=beta.0,col=rgb(0,1,0,.8),lty=2)
matplot(t(beta.save),lty=1,type="l")
abline(h=beta,col=rgb(0,1,0,.8),lty=2)

# --- Fit SPP using cond. output with 2nd stage MCMC ---------------------------
source(here("GlacierBay_Code", "spp_win_2D", "spp.stg2.mcmc.R"))
out.cond.2.full=spp.stg2.mcmc(out.cond.full)
# acceptance rate: 0.018

layout(matrix(1:2,2,1))
plot(out.cond.2.full$beta.0.save,type="l")
abline(h=beta.0,col=rgb(0,1,0,.8),lty=2)
matplot(t(out.cond.2.full$beta.save),lty=1,type="l")
abline(h=beta,col=rgb(0,1,0,.8),lty=2)

effectiveSize(out.cond.2.full$beta.0.save) # 998
effectiveSize(out.cond.2.full$beta.save[1,]) # 974
effectiveSize(out.cond.2.full$beta.save[2,]) # 976

# --- Fit SPP using cond. likelihood (polya-gamma stage 2) ---------------------
source(here("spp_win_2D", "spp.stg2.mcmc.R"))

# using num quad results
tic()
out.cond.pg2 = spp.stg2.mcmc(out.cond.pg)
toc()
# acceptance rate: 0.017

# discard burn-in
beta.save <- out.cond.pg2$beta.save[,-(1:n.burn)]
beta.0.save <- out.cond.pg2$beta.0.save[-(1:n.burn)]

# trace plots
layout(matrix(1:2,2,1))
plot(beta.0.save,type="l")
abline(h=beta.0,col=rgb(0,1,0,.8),lty=2)
matplot(t(beta.save),lty=1,type="l")
abline(h=beta,col=rgb(0,1,0,.8),lty=2)

# posterior summary
beta.save.full <- t(rbind(beta.0.save, beta.save))
vioplot(data.frame(beta.save.full),
        names=expression(beta[0],beta[1],beta[2]),
        ylim = c(-10,5))
abline(h = 0, lty = 2)

apply(beta.save.full,2,mean) 
apply(beta.save.full,2,sd) 
apply(beta.save.full,2,quantile,c(0.025,.975))

effectiveSize(beta.0.save) # 777
effectiveSize(beta.save[1,]) # 735
effectiveSize(beta.save[2,]) # 744

# --- Compare Marginal Posteriors ----------------------------------------------
# pdf("comp_marginals.pdf", width = 14)

out.comp.full=spp.comp.mcmc(s.win,X.win,X.win.full,ds,win.area,n.mcmc,0.1,0.1)

layout(matrix(1:3,1,3))
hist(out.comp.full$beta.0.save[-(1:n.burn)],prob=TRUE,breaks=60,main="",xlab=bquote(beta[0]),
     ylim = c(0,1), cex.axis=1.5, cex.lab=1.3)
lines(density(out.cond.2.full$beta.0.save,n=1000),col="red",lwd=2)
lines(density(out.cond.pg3$beta.0.save[-(1:n.burn)],n=1000),col="green",lwd=2)
abline(v = 4, lwd = 2, lty = 2)
hist(out.comp.full$beta.save[1,-(1:n.burn)],prob=TRUE,breaks=60,main="",xlab=bquote(beta[1]),
     ylim = c(0,1), cex.axis=1.5, cex.lab=1.3)
lines(density(out.cond.2.full$beta.save[1,],n=1000),col="red",lwd=2)
lines(density(out.cond.pg3$beta.save[1,-(1:n.burn)],n=1000),col="green",lwd=2)
abline(v = 2, lwd = 2, lty = 2)
hist(out.comp.full$beta.save[2,-(1:n.burn)],prob=TRUE,breaks=60,main="",xlab=bquote(beta[2]),
     ylim = c(0,1), cex.axis=1.5, cex.lab=1.3)
lines(density(out.cond.2.full$beta.save[2,],n=1000),col="red",lwd=2)
lines(density(out.cond.pg3$beta.save[2,-(1:n.burn)],n=1000),col="green",lwd=2)
abline(v = 1, lwd = 2, lty = 2)
# dev.off()

# --- Posterior for N ----------------------------------------------------------
N.save=rep(0,n.mcmc)

for(k in 1:n.mcmc){
  if(k%%10000==0){cat(k," ")}
  beta.0.tmp=out.cond.pg3$beta.0.save[k]
  beta.tmp=out.cond.pg3$beta.save[,k]
  lam.nowin.int=sum(exp(log(ds)+beta.0.tmp+X.nowin.full%*%beta.tmp))
  N.save[k]=n+rpois(1,lam.nowin.int)
};cat("\n")

par(mfrow = c(1,1))

plot(N.save,type="l")
abline(h=N,col=rgb(0,1,0,.8),lty=2,lwd=2)

hist(N.save[-(1:1000)],breaks=50,prob=TRUE,main="",xlab="N")
abline(v=N,lty=2,lwd=2)

# --- PPD of lambda full area --------------------------------------------------
# idx.sm=seq(1,m,2)
# m.sm=length(idx.sm)
# s.sm=s.full[idx.sm,]
# X.sm=X.full[idx.sm,]
lam.save=matrix(0,m,n.mcmc)
for(k in 1:n.mcmc){
  if(k%%10000==0){cat(k," ")}
  beta.0.tmp=out.cond.2.full$beta.0.save[k]
  beta.tmp=out.cond.2.full$beta.save[,k]
  lam.save[,k]=exp(beta.0.tmp+X.full%*%beta.tmp)
};cat("\n")
lam.mn=apply(lam.save,1,mean)
lam.u=apply(lam.save,1,quantile,.975)
lam.l=apply(lam.save,1,quantile,.025)

# plot(s.sm,exp(beta.0+X.sm%*%beta),xlab="location",ylab=bquote(lambda),col=3,type="l")
# polygon(c(s.sm,rev(s.sm)),c(lam.u,rev(lam.l)),col=rgb(0,0,0,.2),border=NA)
# lines(s.sm,lam.mn,col=1,lwd=2)

lam.ppd.df <- as.data.frame(cbind(s.full, lam.mn))
ggplot(data = lam.ppd.df, aes(x = x, y = y, col = lam.mn)) + 
  geom_point(size = 0.5)

# sample to data level
