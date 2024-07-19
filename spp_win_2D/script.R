library(tidyverse)
library(raster)

# --- Simulate 2D data ---------------------------------------------------------
set.seed(1234)

# create domain
x.domain <- c(0,1)
y.domain <- c(0,1)

n.win <- 2
tot.area <- (x.domain[2] - x.domain[1])*(y.domain[2] - y.domain[1])
tot.win.area <- tot.area/2
win.area <- tot.win.area/n.win
tot.nonwin.area <- tot.area - tot.win.area
nonwin.area <- tot.nonwin.area/n.win

# get quadrature grid
x.m <- 100
y.m <- 100
m <- x.m*y.m

x.full <- seq(x.domain[1], x.domain[2], length.out = x.m)
y.full <- x.full

s.full <- expand.grid(x = x.full, y = y.full)
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

# plot covariates and lambda
full.df <- as.data.frame(cbind(s.full, x1, x2, lam.full))
ggplot(data = full.df, aes(x = x, y = y, col = x)) + 
  geom_point(size = 0.5)
ggplot(data = full.df, aes(x = x, y = y, col = y)) + 
  geom_point(size = 0.5)

ggplot(data = full.df, aes(x = x, y = y, col = lam.full)) + 
  geom_point(size = 0.5)

# # create full raster
# full.df <- full.df %>% rename(z = lam.full)
# full.raster <- rasterFromXYZ(full.df)

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
N=nrow(s.obs)

# plot superpop lambda
superpop.df <- as.data.frame(cbind(x.superpop, y.superpop, lam.superpop))
ggplot(data = superpop.df, aes(x = x.superpop, y = y.superpop, col = lam.superpop)) + 
  geom_point(size = 0.5)

# --- Get windowed data --------------------------------------------------------
obs.win <- ((s.obs[,1] < 0.5) & (s.obs[,2] > 0.5)) | 
  ((s.obs[,1] > 0.5) & (s.obs[,2] < 0.5))
obs.win.idx <- (1:N)[obs.win]

full.win <- ((s.full[,1] < 0.5) & (s.full[,2] > 0.5)) | 
  ((s.full[,1] > 0.5) & (s.full[,2] < 0.5))
full.win.idx <- (1:m)[full.win]
n=length(obs.win.idx)

s.win=s.obs[obs.win.idx,]
X.win=X.obs[obs.win.idx,]

X.win.full=X.full[full.win.idx,]
X.nowin.full=X.full[-full.win.idx,]

# plot windowed data
obs.df <- as.data.frame(s.obs)
ggplot(data = obs.df, aes(x = x.superpop, y = y.superpop,
                          col = factor(obs.win))) + 
  geom_point()

# --- Fit Poisson GLM ----------------------------------------------------------


# --- Fit SPP w/ complete likelihood -------------------------------------------
n.mcmc=100000
source("spp.comp.mcmc.R")
out.comp.full=spp.comp.mcmc(s.win,X.win,X.win.full,ds,win.area,n.mcmc)

layout(matrix(1:2,2,1))
plot(out.comp.full$beta.0.save,type="l")
abline(h=beta.0,col=rgb(0,1,0,.8),lty=2)
matplot(t(out.comp.full$beta.save),lty=1,type="l")
abline(h=beta,col=rgb(0,1,0,.8),lty=2)

# --- Fit SPP w/ conditional likelihood ----------------------------------------
n.mcmc=100000
source("spp.cond.mcmc.R")
out.cond.full=spp.cond.mcmc(s.win,X.win,X.win.full,ds,n.mcmc)

layout(matrix(1:2,2,1))
plot(out.cond.full$beta.0.save,type="l")
abline(h=beta.0,col=rgb(0,1,0,.8),lty=2)
matplot(t(out.cond.full$beta.save),lty=1,type="l")
abline(h=beta,col=rgb(0,1,0,.8),lty=2)

# --- Fit SPP uisng cond. output with 2nd stage MCMC ---------------------------
source("spp.stg2.mcmc.R")
out.cond.2.full=spp.stg2.mcmc(out.cond.full)

layout(matrix(1:2,2,1))
plot(out.cond.2.full$beta.0.save,type="l")
abline(h=beta.0,col=rgb(0,1,0,.8),lty=2)
matplot(t(out.cond.2.full$beta.save),lty=1,type="l")
abline(h=beta,col=rgb(0,1,0,.8),lty=2)

