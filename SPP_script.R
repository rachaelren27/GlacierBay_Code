setwd("/Users/rlr3795/Desktop/GlacierBay_Project")

library(sf)
library(here)
library(tidyverse)
library(terra)
library(raster)
library(hilbertSimilarity)
library(geosphere)
library(tictoc)
library(vioplot)
library(spatstat)
library(rstanarm)
library(reshape2)
library(BayesLogit)
library(mvnfast)

# --- Read in NPS data ---------------------------------------------------------
path <- here("NPS_data", "HARBORSEAL_2007", "seal_locations_final",
             "nonpup_locs")
seal.locs.20070813 <- st_read(dsn = path, layer = "JHI_20070813_nonpup_locs")

path <- here("NPS_data", "HARBORSEAL_2007", "footprints")
footprint.20070813 <- st_read(dsn = path, layer = "JHI_20070813_footprint")

path <- here("NPS_data", "HARBORSEAL_2007", "survey_polygons")
survey.poly.20070813 <- st_read(dsn = path, layer = "JHI_20070813_surveyboundary")

# convert CRS
survey.poly <- st_transform(survey.poly.20070813$geometry, 
                            CRS("+proj=longlat +datum=WGS84"))

seal.locs <- st_transform(seal.locs.20070813$geometry,
                          CRS("+proj=longlat +datum=WGS84"))

footprint <- st_transform(footprint.20070813$geometry, 
                          CRS("+proj=longlat +datum=WGS84"))

# crop footprint
footprint <- st_intersection(footprint, survey.poly)

# prepare windows
survey.poly.mat <- survey.poly[[1]][[1]]
survey.win <- owin(poly = data.frame(x=rev(survey.poly.mat[,1]),
                                     y=rev(survey.poly.mat[,2])))

footprints <- lapply(1:length(footprint), function(i) {
  footprint.mat <- footprint[[i]][[1]]
  if(class(footprint.mat)[1] == "list"){
    footprint.mat <- footprint.mat[[1]]
  }
  owin(poly = data.frame(x=rev(footprint.mat[,1]),
                         y=rev(footprint.mat[,2])))
})
footprint.win <- do.call(union.owin, footprints)

# plot
ggplot() + 
  geom_sf(data = survey.poly) + 
  geom_sf(data = footprint) + 
  geom_sf(data = seal.locs, size = 0.5)

# --- Read in covariates -------------------------------------------------------
# read in bathymetry
bath.rast <- raster(here("covariates", "bathymetry.tiff"))

# crop using survey boundary
bath.rast.survey <- raster::crop(bath.rast, extent(survey.poly.mat))
bath.rast.survey <- raster::mask(bath.rast.survey, as(survey.poly, 'Spatial'))
# length: 1619618

plot(bath.rast.survey)
plot(survey.poly, add = TRUE)

# calculate distance from southern boundary (glacier)
ggplot() + 
  geom_sf(data = survey.poly) + 
  geom_point(aes(x = -137.1311, y = 58.84288), color = "red") # westmost point

ggplot() + 
  geom_line(data = glacier.poly, aes(x = V1, y = V2)) + 
  geom_point(aes(x = glacier.poly[101,1], y = glacier.poly[101, 2]), color = "red")

survey.poly.df <- as.data.frame(survey.poly.mat)
glacier.poly <- survey.poly.df %>% filter(V2 < 58.84288)
glacier.poly <- as.matrix(glacier.poly[-(1:100),]) # found index 101 using localMinima

ggplot() + 
  geom_sf(data = survey.poly) + 
  geom_point(data = glacier.poly, aes(x = V1, y = V2), color = "red")

seal.mat <- as.matrix(st_coordinates(seal.locs))
bath.rast <- na.omit(values(bath.rast.survey))

seal.glac.dist <- dist2Line(seal.mat, glacier.poly) # in meters

bath.survey.idx <- which(!is.na(values(bath.rast.survey)))
full.coord <- xyFromCell(bath.rast.survey, bath.survey.idx)

full.glac.dist <- dist2Line(full.coord, glacier.poly) # takes a while

cor(na.omit(values(bath.rast.survey)), full.glac.dist[,1]) # -0.183

# --- Calculate areas ----------------------------------------------------------
tot.area <- area.owin(survey.win)
tot.win.area <- area.owin(footprint.win)
tot.nonwin.area <- tot.area - tot.win.area
n.win <- length(footprint)
ex.win <- owin(poly = data.frame(x = rev(footprint[[1]][[1]][,1]),
                                 y = rev(footprint[[1]][[1]][,2])))
win.area <- area.owin(ex.win) # approx. bc windows not equally sized

ds <- res(bath.rast.survey)[1]*res(bath.rast.survey)[2]
  
# --- Set X matrices -----------------------------------------------------------
glac.dist <- full.glac.dist[,1]

seal.full.idx <- cellFromXY(bath.rast.survey, seal.mat)
row.counts <- table(factor(seal.full.idx, levels = 1:length(bath.rast.survey)))
bath.full <- cbind(values(bath.rast.survey), row.counts)
bath <- na.omit(bath.full)
X.full <- cbind(bath, glac.dist)
seal.idx <- c()
for(i in 1:nrow(X.full)){
  if(X.full[i,2] != 0){
    seal.idx <- c(seal.idx, rep(i, times = X.full[i,2]))
  }
}
X.full <- scale(X.full[,-2])

win.idx <- which(inside.owin(full.coord[,1], full.coord[,2], footprint.win))
X.win.full <- X.full[win.idx,]

X.obs <- X.full[seal.idx,]
n <- length(seal.idx)

# --- Fit SPP w/ Complete Likelihood -------------------------------------------
n.mcmc=10000
source(here("GlacierBay_Code", "spp_win_2D", "spp.comp.mcmc.R"))
tic()
out.comp.full=spp.comp.mcmc(seal.mat,X.obs,X.win.full,ds,win.area,n.mcmc)
toc() # 543.252 sec elapsed (~9 min)

# discard burn-in
n.burn <- 0.1*n.mcmc
beta.save <- out.comp.full$beta.save[,-(1:n.burn)]
beta.0.save <- out.comp.full$beta.0.save[-(1:n.burn)]

# trace plots
layout(matrix(1:2,2,1))
plot(beta.0.save,type="l")
matplot(t(beta.save),lty=1,type="l", col = c("black", "red"))

# posterior summary
beta.save.full <- t(rbind(beta.0.save, beta.save))
vioplot(data.frame(beta.save.full),
        names=expression(beta[0],beta[1],beta[2]),
        ylim = c(-10,5))
abline(h = 0, lty = 2)

apply(beta.save.full,2,mean) 
apply(beta.save.full,2,sd) 
apply(beta.save.full,2,quantile,c(0.025,.975))

# --- Fit SPP w/ cond. likelihood (num quad stage 1) ---------------------------
source(here("GlacierBay_Code", "spp_win_2D", "spp.cond.mcmc.R"))
tic()
out.cond.full=spp.cond.mcmc(seal.mat,X.obs,X.win.full,ds,n.mcmc)
toc() # 282.763 sec (~4.7 min)

# discard burn-in
beta.save <- out.cond.full$beta.save[,-(1:n.burn)]
beta.0.save <- out.cond.full$beta.0.save[-(1:n.burn)]

# trace plots
layout(matrix(1:2,2,1))
plot(out.cond.full$beta.0.save,type="l")
matplot(t(out.cond.full$beta.save),lty=1,type="l")

# posterior summary
beta.save.full <- t(rbind(beta.0.save, beta.save))
vioplot(data.frame(beta.save.full),
        names=expression(beta[0],beta[1],beta[2]),
        ylim = c(-10,5))
abline(h = 0, lty = 2)

apply(beta.save.full,2,mean) 
apply(beta.save.full,2,sd) 
apply(beta.save.full,2,quantile,c(0.025,.975))

# --- Fit SPPS using cond. likelihood (Bernoulli GLM stage 1) ------------------

# obtain background sample
  n.bg <- 50000
  bg.pts <- rpoint(n.bg, win = footprint.win)
  
  ggplot() + 
    geom_sf(data = survey.poly) + 
    geom_sf(data = footprint) + 
    geom_point(aes(x = bg.pts$x, y = bg.pts$y), size = 0.3) + 
    geom_sf(data = seal.locs, size = 0.3, col = "red")
  
  # prepare covariates for background sample 
  bg.mat <- cbind(bg.pts$x, bg.pts$y)
  
  bg.full.idx <- cellFromXY(bath.rast.survey, bg.mat)
  row.counts <- table(factor(bg.full.idx, levels = 1:length(bath.rast.survey)))
  bath.full <- cbind(values(bath.rast.survey), row.counts)
  bath <- na.omit(bath.full)
  X.full <- cbind(bath, glac.dist)
  bg.idx <- c()
  for(i in 1:nrow(X.full)){
    if(X.full[i,2] != 0){
      bg.idx <- c(bg.idx, rep(i, times = X.full[i,2]))
    }
  }
  # 10000 -> 9802 background points (some correspond to NA in bath raster)
  # 50000 -> 49008 bg pts
  X.full <- scale(X.full[,-2]) 

X.obs <- X.full[c(seal.idx, bg.idx),] 

y.binary <- rep(0, n + length(bg.idx))
y.binary[1:n] <- 1

bern.rsf.df <- data.frame(y = y.binary, bath = X.obs[,1], glac.dist = X.obs[,2])
tic()
out.bern.cond <- stan_glm(y ~ bath + glac.dist, family=binomial(link="logit"), data=bern.rsf.df,
                 iter = 100000, chains = 1)
toc()
# 9802 bg pts: 376.625 sec (~6.3 min)
# 49008 bg pts: 2079 sec (~34.7 min)

# plot comparison (Bernoulli GLM vs num quad)
beta.save.full.cond1 <- cbind(beta.save.full, rep(1, nrow(beta.save.full)))
beta.save.full.stan1 <- cbind(as.matrix(out.bern.cond), rep(0, nrow(as.matrix(out.bern.cond))))
beta.save.full.stage1 <- as.data.frame(rbind(beta.save.full.cond1, beta.save.full.stan1)[,-1])

beta.save.stage1.long <- melt(beta.save.full.stage1, id.vars = "V3")

pdf("stage1_compare.pdf")
ggplot(beta.save.stage1.long, aes(x = variable, y = value, fill = as.factor(V3))) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(x = "Coefficients",
       y = "Values",
       fill = "Model") + 
  scale_x_discrete(labels = c("bathymetry", "glacier distance")) + 
  scale_fill_manual(values = c("#00BFC4", "#F8766D"), labels = c("stan_glm", "num quad"))
dev.off()

# --- Fit SPP using cond. likelihood (Polya-gamma stage 1) ---------------------
# obtain background sample
n.bg <- 10000
bg.pts <- rpoint(n.bg, win = footprint.win)

ggplot() + 
  geom_sf(data = survey.poly) + 
  geom_sf(data = footprint) + 
  geom_point(aes(x = bg.pts$x, y = bg.pts$y), size = 0.3) + 
  geom_sf(data = seal.locs, size = 0.3, col = "red")

# prepare covariates for background sample 
bg.mat <- cbind(bg.pts$x, bg.pts$y)

bg.full.idx <- cellFromXY(bath.rast.survey, bg.mat)
row.counts <- table(factor(bg.full.idx, levels = 1:length(bath.rast.survey)))
bath.full <- cbind(values(bath.rast.survey), row.counts)
bath <- na.omit(bath.full)
X.full <- cbind(bath, glac.dist)
bg.idx <- c()
for(i in 1:nrow(X.full)){
  if(X.full[i,2] != 0){
    bg.idx <- c(bg.idx, rep(i, times = X.full[i,2]))
  }
}
# 10000 -> 9802 background points (some correspond to NA in bath raster)
# 50000 -> 49008 bg pts
X.full <- scale(X.full[,-2]) 

X.obs <- X.full[c(seal.idx, bg.idx),] 

y.binary <- rep(0, n + length(bg.idx))
y.binary[1:n] <- 1

source(here("GlacierBay_Code", "Polya_Gamma.R"))
p <- ncol(X.obs)
mu.beta <- rep(0, p)
sigma.beta <- diag(100, p)
tic()
beta.save.pg <- polya_gamma(y.binary, X.obs, mu.beta, sigma.beta, 10000)
toc() # 2.4 hours

# posterior summary
plot(beta.save.pg$beta[1,], type = "l")
plot(beta.save.pg$beta[2,], type = "l")

# --- Fit SPP using cond. output (num quad stage 2) ----------------------------
source(here("GlacierBay_Code", "spp_win_2D", "spp.stg2.mcmc.R"))
tic()
out.cond.2.full=spp.stg2.mcmc(out.cond.full)
toc() # 138.445 sec (~2.3 min)

# discard burn-in
beta.save <- out.cond.2.full$beta.save[,-(1:n.burn)]
beta.0.save <- out.cond.2.full$beta.0.save[-(1:n.burn)]

# trace plots
layout(matrix(1:2,2,1))
plot(out.cond.2.full$beta.0.save,type="l")
matplot(t(out.cond.2.full$beta.save),lty=1,type="l")

# posterior summary
beta.save.full <- t(rbind(beta.0.save, beta.save))
vioplot(data.frame(beta.save.full),
        names=expression(beta[0],beta[1],beta[2]),
        ylim = c(-10,5))
abline(h = 0, lty = 2)

apply(beta.save.full,2,mean) 
apply(beta.save.full,2,sd) 
apply(beta.save.full,2,quantile,c(0.025,.975))

# --- Compare Marginal Posteriors ----------------------------------------------
# layout(matrix(1:3,1,3))
hist(out.comp.full$beta.0.save,prob=TRUE,breaks=60,main="",xlab=bquote(beta[0]),
     ylim = c(0,1))
lines(density(out.cond.full$beta.0.save,n=1000),col=2,lwd=1)
lines(density(out.cond.2.full$beta.0.save,n=1000,adj=2),col=3,lwd=1)
hist(out.comp.full$beta.save[1,],prob=TRUE,breaks=60,main="",xlab=bquote(beta[1]),
     ylim = c(0,5))
lines(density(out.cond.full$beta.save[1,],n=1000),col=2,lwd=1)
lines(density(out.cond.2.full$beta.save[1,],n=1000,adj=2),col=3,lwd=1)
hist(out.comp.full$beta.save[2,],prob=TRUE,breaks=60,main="",xlab=bquote(beta[2]),
     ylim = c(0,1.5))
lines(density(out.cond.full$beta.save[2,],n=1000),col=2,lwd=1)
lines(density(out.cond.2.full$beta.save[2,],n=1000,adj=2),col=3,lwd=1)

# --- Posterior for N ----------------------------------------------------------
N.save=rep(0,n.mcmc)

X.nowin.full <- X.full[-win.idx,]
n <- nrow(seal.mat)

tic()
for(k in 1:n.mcmc){
  if(k%%10000==0){cat(k," ")}
  beta.0.tmp=out.cond.2.full$beta.0.save[k]
  beta.tmp=out.cond.2.full$beta.save[,k]
  lam.nowin.int=sum(exp(log(ds)+beta.0.tmp+X.nowin.full%*%beta.tmp)) # can parallelize
  N.save[k]=n+rpois(1,lam.nowin.int)
};cat("\n")
toc()

par(mfrow = c(1,1))

# discard burn-in
N.save <- N.save[-(1:n.burn)]

plot(N.save,type="l")
hist(N.save,breaks=50,prob=TRUE,main="",xlab="N")

# posterior summary
mean(N.save)
sd(N.save)
quantile(N.save, c(0.025, 0.975))

