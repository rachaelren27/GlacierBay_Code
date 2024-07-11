####
####  Load Libraries
####

library(raster)
library(MASS)
library(spatstat)
library(vioplot)
library(rjags)
library(mvnfast)
library(polyCub)
library(tidyverse)
library(geometry)
load("mtnlion.RData")
names(mtnlion.df)=c("x","y","time")
n=dim(mtnlion.df)[1]


# --- Bernoulli Regression Approximation Setup ---------------------------------

# Set spatial support
hr.tmp=draw.contour(mtnlion.df[,1:2],bw=NULL,alpha=.99) # read-in
n.contour=length(hr.tmp)
mat=matrix(0, nrow = n.contour, ncol = 1)
for(i in 1:n.contour) {
  mat[i]=length(hr.tmp[[i]]$x)
}
hr=owin(poly=list(list(x=rev(hr.tmp[[which.max(mat)]]$x[-1]),
                       y=rev(hr.tmp[[which.max(mat)]]$y[-1]))))
pts=inside.owin(mtnlion.df$x,mtnlion.df$y,hr)

# Obtain data in Bernoulli form
pts=data.frame(x=mtnlion.df$x,y=mtnlion.df$y,inout=pts)
pts=subset(pts,pts$inout=="TRUE") # set observed locations to 1
pts.ppp=ppp(x=pts$x,y=pts$y,window=hr)
set.seed(101)
n.bg=1000
bg.ppp=rpoint(n.bg,win=hr) # obtain background sample

# plot comparison
layout(matrix(1:2,1,2,byrow=TRUE))
plot(pts.ppp,type="n",cex=1.5,pch=19,col=rgb(0,0,0,alpha=.4),cex.axis=1.2,
     cex.lab=1.5,cex.main=1.5,xlab="",ylab="",asp=TRUE,main="telemetry data")
points(pts.ppp,pch=19,col=rgb(0,0,0,alpha=.6))
plot(bg.ppp,type="n",cex=1.5,pch=19,col=rgb(0,0,0,alpha=.4),cex.axis=1.2,
     cex.lab=1.5,cex.main=1.5,xlab="",ylab="",asp=TRUE,main="background sample")
points(bg.ppp,pch=19,col=rgb(0,0,0,alpha=.6))


# --- Poisson Regression Approximation Setup -----------------------------------
grid.hr.TF=inside.owin(grid.centers[,1],grid.centers[,2],hr)
elev.NA.rast=elevation.rast
values(elev.NA.rast)[!grid.hr.TF]=NA

mtnlion.count.rast=elev.NA.rast
values(mtnlion.count.rast)[grid.hr.TF]=0
mtnlion.cells=cellFromXY(elev.NA.rast,mtnlion.df[,1:2]) # get cell index from mountain lion locs
mtnlion.table=table(mtnlion.cells)
mtnlion.cell.idx=as.numeric(attr(table(mtnlion.cells),"dimnames")$mtnlion.cells)
values(mtnlion.count.rast)[mtnlion.cell.idx]=mtnlion.table

# plot comparison
layout(matrix(1:2,1,2))
plot(pts.ppp,type="n",cex=1.5,pch=19,col=rgb(0,0,0,alpha=.4),cex.axis=1.2,
     cex.lab=1.5,cex.main=1.5,xlab="",ylab="",asp=TRUE,main="telemetry data")
points(pts.ppp,pch=19,col=rgb(0,0,0,alpha=.6))
lines(hr$bdry[[1]],lwd=3,col=rgb(0,0,0,alpha=.4))
plot(pts.ppp,type="n",cex=1.5,pch=19,col=0,cex.axis=1.2,cex.lab=1.5,
     cex.main=1.5,xlab="",ylab="",asp=TRUE,main="gridded counts")
image(mtnlion.count.rast,col=gray.colors.rev(100),add=TRUE)
lines(hr$bdry[[1]],lwd=3,col=rgb(0,0,0,alpha=.4))

# --- Prepare and Plot Covariates ----------------------------------------------
elev.NA.rast=elevation.rast
slope.NA.rast=slope.rast
exposure.NA.rast=exposure.rast
values(elev.NA.rast)[!grid.hr.TF]=NA
values(slope.NA.rast)[!grid.hr.TF]=NA
values(exposure.NA.rast)[!grid.hr.TF]=NA
values(elev.NA.rast)[grid.hr.TF]=scale(values(elev.NA.rast)[grid.hr.TF])
values(slope.NA.rast)[grid.hr.TF]=scale(values(slope.NA.rast)[grid.hr.TF])
values(exposure.NA.rast)[grid.hr.TF]=scale(values(exposure.NA.rast)[grid.hr.TF])

layout(matrix(1:3,1,3))
plot(pts.ppp,type="n",cex=1.5,pch=19,col=0,cex.axis=1.2,cex.lab=1.5,cex.main=1.5,
     xlab="",ylab="",asp=TRUE,main="elevation")
image(elev.NA.rast,col=gray.colors.rev(100),add=TRUE)
lines(hr$bdry[[1]],lwd=3,col=rgb(0,0,0,alpha=.4))
plot(pts.ppp,type="n",cex=1.5,pch=19,col=0,cex.axis=1.2,cex.lab=1.5,cex.main=1.5,
     xlab="",ylab="",asp=TRUE,main="slope")
image(slope.NA.rast,col=gray.colors.rev(100),add=TRUE)
lines(hr$bdry[[1]],lwd=3,col=rgb(0,0,0,alpha=.4))
plot(pts.ppp,type="n",cex=1.5,pch=19,col=0,cex.axis=1.2,cex.lab=1.5,cex.main=1.5,
     xlab="",ylab="",asp=TRUE,main="exposure")
image(exposure.NA.rast,col=gray.colors.rev(100),add=TRUE)
lines(hr$bdry[[1]],lwd=3,col=rgb(0,0,0,alpha=.4))


# --- Fit Models using GLM function --------------------------------------------
# prepare covariates for Poisson regression
mtnlion.count.rast=elev.NA.rast
values(mtnlion.count.rast)[grid.hr.TF]=0
mtnlion.cells=cellFromXY(elev.NA.rast,mtnlion.df[,1:2])
mtnlion.table=table(mtnlion.cells)
mtnlion.cell.idx=as.numeric(attr(table(mtnlion.cells),"dimnames")$mtnlion.cells)
values(mtnlion.count.rast)[mtnlion.cell.idx]=mtnlion.table
rsf.1.df=data.frame(y=values(mtnlion.count.rast),elev=values(elev.NA.rast),
                    slope=values(slope.NA.rast),exposure=values(exposure.NA.rast))
cor(rsf.1.df[!is.na(rsf.1.df$y),]) # check for colinearity

# prepare covariates for Bernoulli regression
bg.cells=cellFromXY(elev.NA.rast,cbind(bg.ppp$x,bg.ppp$y))
y.binary=rep(0,n+n.bg)
elev.binary=rep(0,n+n.bg)
slope.binary=rep(0,n+n.bg)
exposure.binary=rep(0,n+n.bg)
y.binary[1:n]=1
elev.binary[1:n]=values(elev.NA.rast)[mtnlion.cells]
slope.binary[1:n]=values(slope.NA.rast)[mtnlion.cells]
exposure.binary[1:n]=values(exposure.NA.rast)[mtnlion.cells]
elev.binary[-(1:n)]=values(elev.NA.rast)[bg.cells]
slope.binary[-(1:n)]=values(slope.NA.rast)[bg.cells]
exposure.binary[-(1:n)]=values(exposure.NA.rast)[bg.cells]
rsf.2.df=data.frame(y=y.binary,elev=elev.binary,slope=slope.binary,exposure=exposure.binary)

# fit Poisson GLM
summary(glm(y~elev+slope+exposure,family=poisson(link="log"),data=rsf.1.df))

# fit Bernoulli GLM
summary(glm(y~elev+slope+exposure,family=binomial(link="logit"),data=rsf.2.df))


# --- Fit Bernoulli Model using JAGS -------------------------------------------
# data preparation
rsf.noNA.df=rsf.2.df[!apply(is.na(rsf.2.df),1,any),]  # omit points w/ NA values for covariates

y=rsf.noNA.df$y
X=model.matrix(~elev+slope+exposure,data=rsf.noNA.df)
m=length(y)

# model specification
m.jags <-"
  model{
    for(i in 1:m){
      y[i] ~ dbern(p[i])
      logit(p[i]) <- b0+b1*X[i,2]+b2*X[i,3]+b3*X[i,4]
    }
    b0 ~ dnorm(mu,tau)
    b1 ~ dnorm(mu,tau)
    b2 ~ dnorm(mu,tau)
    b3 ~ dnorm(mu,tau)
}
"

mod<-textConnection(m.jags)

# set hyperpriors
mu=0
tau=1/2.25
m.out<-jags.model(mod,data=list('y'=y,'m'=m,'X'=X,'mu'=mu,'tau'=tau),n.chains=1,
                  n.adapt=0) # build model and algorithm

n.mcmc=50000
n.burn=round(.2*n.mcmc)
update(m.out,n.burn) # perform burn-in
m.samples=jags.samples(m.out,c('b0','b1','b2','b3'),n.mcmc) # fit model post-burn-in

beta.post.mat=cbind(m.samples$b0[1,,1],m.samples$b1[1,,1],m.samples$b2[1,,1],m.samples$b3[1,,1])
matplot(beta.post.mat,type="l",lty=1,ylab=bquote(beta))

# plot marginal posterior distributions
vioplot(data.frame(beta.post.mat),names=expression(beta[0],beta[1],beta[2],beta[3]))
abline(h=0,col=8)

# posterior summary
apply(beta.post.mat,2,mean) # marginal posterior means for beta
apply(beta.post.mat,2,sd) # marginal posterior std devs for beta
apply(beta.post.mat,2,quantile,c(0.025,.975)) # marginal posterior 95% CI for beta

# compare to GLM function
summary(glm(y~elev+slope+exposure,family=binomial(link="logit"),data=rsf.2.df))
# very similar!


# --- Fit Poisson Model using JAGS ---------------------------------------------
# data preparation
rsf.count.noNA.df=rsf.1.df[!apply(is.na(rsf.1.df),1,any),]  # omit points w/ NA values for covariates

y=rsf.count.noNA.df$y
X=model.matrix(~elev+slope+exposure,data=rsf.count.noNA.df)
m=length(y)

# model specification
m.jags <-"
  model{
    for(i in 1:m){
      y[i] ~ dpois(lambda[i])
      log(lambda[i]) <- b0+b1*X[i,2]+b2*X[i,3]+b3*X[i,4]
    }
    b0 ~ dnorm(mu,tau)
    b1 ~ dnorm(mu,tau)
    b2 ~ dnorm(mu,tau)
    b3 ~ dnorm(mu,tau)
}
"

mod<-textConnection(m.jags)

mu=0
tau=1/100
m.out<-jags.model(mod,data=list('y'=y,'m'=m,'X'=X,'mu'=mu,'tau'=tau),n.chains=1,
                  n.adapt=0) # build model and algorithm

n.mcmc=50000
n.burn=round(.2*n.mcmc)
update(m.out,n.burn) # perform burn-in
m.samples=jags.samples(m.out,c('b0','b1','b2','b3'),n.mcmc) # fit model post-burn-in

beta.post.mat=cbind(m.samples$b0[1,,1],m.samples$b1[1,,1],m.samples$b2[1,,1],m.samples$b3[1,,1])
matplot(beta.post.mat,type="l",lty=1,ylab=bquote(beta))

# plot marginal posterior distributions
vioplot(data.frame(beta.post.mat),names=expression(beta[0],beta[1],beta[2],beta[3]))
abline(h=0,col=8)

# posterior summary
apply(beta.post.mat,2,mean) # marginal posterior means for beta
apply(beta.post.mat,2,sd) # marginal posterior std devs for beta
apply(beta.post.mat,2,quantile,c(0.025,.975)) # marginal posterior 95% CI for beta

## comparison to MLE
summary(glm(y~elev+slope+exposure,family=poisson(link="log"),data=rsf.1.df))
# std dev estimates are smaller
# means are significantly different
# testing hypotheses at alpha = 0.05 levels still give us the same outcomes

## comparison to Bayesian Bernoulli GLM
# std dev estimates are smaller
# means are significantly different
# testing hypotheses at alpha = 0.05 levels still give us the same outcomes


# --- Numerical Quadrature -----------------------------------------------------
# x is a px1 vector of covariates at prespecified location
# beta is a px1 vector of coefficients
exp.rsf <- function(x, beta, log = FALSE){
  out <- dot(x,beta)
  if(log == FALSE){
    out <- exp(out)
  }
  return(out)
}

# X.grid is an n.grid x p matrix of covariates 
# obs is an n.obs x 1 vector of indices that map observed loc to corresponding row in X.grid
# beta is a px1 vector of coefficients
SPP.joint.log.lik.numquad <- function(X.grid, obs, beta, cell.size){
  n.obs <- length(obs)
  X.grid.obs <- X.grid[obs,]
  num <- sum(apply(X.grid.obs, 1, exp.rsf, beta = beta, log = TRUE)) 
  denom <- n.obs*log(cell.size) + n.obs*log(sum(apply(X.grid, 1, exp.rsf, beta = beta)))
  return(num - denom)
}
# multiply by cell size
# take log for stability

mcmc <- function(n.mcmc, mu.beta, sigma.beta, beta.0, X.grid, obs, cell.size){
  p <- length(beta.0)
  beta.save <- matrix(nrow = n.mcmc, ncol = p)
  beta.save[1,] <- beta.0
  sd.tune <- 0.01
  beta <- beta.0
  
  for(k in 2:n.mcmc){
    beta.prop <- as.vector(rmvn(1, mu = beta, sigma = sd.tune*diag(p)))
    mh <- (SPP.joint.log.lik.numquad(X.grid, obs, beta.prop, cell.size) + 
             dmvn(beta.prop, mu.beta, sigma.beta, log = TRUE)) - 
      (SPP.joint.log.lik.numquad(X.grid, obs, beta, cell.size) + 
         dmvn(beta, mu.beta, sigma.beta, log = TRUE))
    if(runif(1) < exp(mh)){
      beta <- beta.prop
    } 
    beta.save[k,] <- beta
    if(k %% 100 == 0){
      print(k)
    }
  }
  
  return(beta.save)
}

pois.glm.out <- glm(y~elev+slope+exposure,family=poisson(link="log"),data=rsf.1.df)
beta.0 <- pois.glm.out$coefficients[-1] # throw away intercept
p <- length(beta.0)

# prepare X.grid and obs
load("mtnlion.RData")
values(elevation.rast)=scale(values(elevation.rast))
values(slope.rast)=scale(values(slope.rast))
values(exposure.rast)=scale(values(exposure.rast))
X.grid <- cbind(values(elevation.rast), values(slope.rast), values(exposure.rast))

# sum(is.na(values(elevation.rast))) # 0
# sum(is.na(values(slope.rast))) # 148
# sum(is.na(values(exposure.rast))) # 148

obs <- cellFromXY(elevation.rast,mtnlion.df[,1:2])
# X.grid <- cbind(X.grid, ifelse(1:nrow(X.grid) %in% obs, 1, 0))
# X.grid <- na.omit(X.grid)
# obs <- which(X.grid[,4] == 1)
# X.grid <- X.grid[,-4]
row.counts <- table(factor(obs, levels = 1:nrow(X.grid)))
X.grid <- cbind(X.grid, row.counts)
X.grid <- na.omit(X.grid)
obs <- c()
for(i in 1:nrow(X.grid)){
  if(X.grid[i,4] != 0){
    obs <- c(obs, rep(i, times = X.grid[i,4]))
  }
}
X.grid <- X.grid[,-4]

# run mcmc
n.mcmc <- 5000
cell.res <- res(elevation.rast)
cell.size <- cell.res[1] * cell.res[2]
beta.samp <- mcmc(n.mcmc, rep(0, p), 10*diag(p), beta.0, X.grid, obs, cell.size)
beta.samp <- beta.samp[-(1:n.mcmc*0.2),] # discard burn-in

# trace plots
plot(beta.samp[,1], type = "l", col = "red", ylim = c(-1, 1))
lines(beta.samp[,2], col = "blue")
lines(beta.samp[,3], col = "green")

# posterior summary
vioplot(data.frame(beta.samp),names=expression(beta[1],beta[2],beta[3]),
        ylim = c(-1,1))
abline(h=0,col=8)

# posterior summary
apply(beta.samp,2,mean) # marginal posterior means for beta
apply(beta.samp,2,sd) # marginal posterior std devs for beta
apply(beta.samp,2,quantile,c(0.025,.975))

# obtain posterior predictive?
pred.ind <- 100
pred.loc <- xyFromCell(elevation.rast, cell = 100)

SPP.lik.numquad <- function(X.grid, pred.ind, beta){
  X.grid.pred <- X.grid[pred.ind,]
  num <- exp.rsf(X.grid.pred, beta)
  denom <- sum(apply(X.grid, 1, exp.rsf, beta = beta))
  return(num/denom)
}

s.pred.lik <- apply(beta.samp, 1, SPP.lik.numquad, X.grid = X.grid, pred.ind = pred.ind)
hist(s.pred.lik, freq = FALSE)
mean(s.pred.lik)
# want to sample out to the data level
