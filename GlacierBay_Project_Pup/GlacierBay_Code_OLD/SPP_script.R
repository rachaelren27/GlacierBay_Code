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
library(pgdraw)
library(coda)
library(viridis)

set.seed(1234)
load(here("SPP_0618.RData"))

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
  owin(poly = data.frame(x=footprint.mat[,1],
                         y=footprint.mat[,2]))
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

# read in ice
ice.rast <- raster(here("GlacierBay_Code", "GlacierBay_Project_Pup", "covariates",
                        "ice_20070618_FRK_clipped.tiff"))

ice.rast <- raster::crop(ice.rast, extent(survey.poly.mat))
ice.rast <- raster::mask(ice.rast, as(survey.poly, 'Spatial'))

plot(ice.rast)

# # calculate distance from southern boundary (glacier)
# ggplot() + 
#   geom_sf(data = survey.poly) + 
#   geom_point(aes(x = -137.1311, y = 58.84288), color = "red") # westmost point
# 
# ggplot() + 
#   geom_line(data = glacier.poly, aes(x = V1, y = V2)) + 
#   geom_point(aes(x = glacier.poly[101,1], y = glacier.poly[101, 2]), color = "red")
# 
# survey.poly.df <- as.data.frame(survey.poly.mat)
# glacier.poly <- survey.poly.df %>% filter(V2 < 58.84288)
# glacier.poly <- as.matrix(glacier.poly[-(1:100),]) # found index 101 using localMinima
# 
# ggplot() + 
#   geom_sf(data = survey.poly) + 
#   geom_point(data = glacier.poly, aes(x = V1, y = V2), color = "red")
# 
# seal.mat <- as.matrix(st_coordinates(seal.locs))
# bath.rast <- na.omit(values(bath.rast.survey))
# 
# seal.glac.dist <- dist2Line(seal.mat, glacier.poly) # in meters
# 
# bath.survey.idx <- which(!is.na(values(bath.rast.survey)))
# full.coord <- xyFromCell(bath.rast.survey, bath.survey.idx)
# 
# full.glac.dist <- dist2Line(full.coord, glacier.poly) # takes a while
# 
# glac.dist.df <- data.frame(x = full.coord[,1], y = full.coord[,2],
#                            z = full.glac.dist[,1])
# glac.dist.rast <- rasterFromXYZ(glac.dist.df)
# writeRaster(glac.dist.rast, filename = "glacier_dist.tiff", format = "GTiff")


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
# glac.dist <- full.glac.dist[,1]
# 
# seal.full.idx <- cellFromXY(bath.rast.survey, seal.mat)
# row.counts <- table(factor(seal.full.idx, levels = 1:length(bath.rast.survey)))
# bath.full <- cbind(values(bath.rast.survey), row.counts)
# bath <- na.omit(bath.full)
# X.full <- cbind(bath, glac.dist)
# seal.idx <- c()
# for(i in 1:nrow(X.full)){
#   if(X.full[i,2] != 0){
#     seal.idx <- c(seal.idx, rep(i, times = X.full[i,2]))
#   }
# }
# X.full <- scale(X.full[,-2])
# 
# win.idx <- which(inside.owin(full.coord[,1], full.coord[,2], footprint.win))
# X.win.full <- X.full[win.idx,]
# 
# X.obs <- X.full[seal.idx,]
# n <- length(seal.idx)

ice.idx <- which(!is.na(values(ice.rast)))
ice.full.coord <- xyFromCell(ice.rast, ice.idx)

bath.idx <- cellFromXY(bath.rast.survey, ice.full.coord)
bath <- values(bath.rast.survey)[bath.idx]
glac.dist.idx <- cellFromXY(glac.dist.rast, ice.full.coord)
glac.dist <- values(glac.dist.rast)[glac.dist.idx]

seal.full.idx <- cellFromXY(ice.rast, seal.mat)
row.counts <- table(factor(seal.full.idx, levels = 1:length(ice.rast)))
ice.full.counts <- cbind(values(ice.rast), row.counts)
ice <- na.omit(ice.full.counts)
X.full <- cbind(ice, bath, glac.dist, ice.full.coord)
seal.idx <- c()
for(i in 1:nrow(X.full)){
  if(X.full[i,2] != 0){
    seal.idx <- c(seal.idx, rep(i, times = X.full[i,2]))
  }
}
X.full <- na.omit(X.full)
full.coord <- X.full[,5:6]
X.full <- scale(X.full[,-2])

# check how many seals on 0 ice
seal.ice.idx <- cellFromXY(ice.rast, seal.mat)
num.seal.0.ice <- sum(na.omit(values(ice.rast)[seal.ice.idx] == 0))

win.idx <- which(inside.owin(full.coord[,1], full.coord[,2], footprint.win)) # 15 sec

X.win.full <- X.full[win.idx,]
X.nowin.full <- X.full[-win.idx,]

X.obs <- X.full[seal.idx,]
n <- nrow(X.obs)

# --- Fit SPP w/ Complete Likelihood -------------------------------------------
n.mcmc=100000
source(here("GlacierBay_Code", "GlacierBay_Project_Pup", "GlacierBay_Code_OLD",
            "spp_win_2D", "spp.comp.mcmc.R"))
tic()
out.comp.full=spp.comp.mcmc(seal.mat,X.obs,X.win.full,ds,win.area,n.mcmc,0.1,0.1)
toc() # 543.252 sec elapsed (~9 min)

# discard burn-in
n.burn <- 0.1*n.mcmc
beta.save.full.lik <- out.comp.full$beta.save[,-(1:n.burn)]
beta.0.save.full.lik <- out.comp.full$beta.0.save[-(1:n.burn)]

# trace plots
layout(matrix(1:2,2,1))
plot(beta.0.save.full.lik,type="l")
matplot(t(beta.save.full.lik),lty=1,type="l", col = c("black", "red"))

# posterior summary
beta.save.full <- t(rbind(beta.0.save.full.lik, beta.save.full.lik))
vioplot(data.frame(beta.save.full),
        names=expression(beta[0],beta[1],beta[2]),
        ylim = c(-10,5))
abline(h = 0, lty = 2)

apply(beta.save.full,2,mean) 
apply(beta.save.full,2,sd) 
apply(beta.save.full,2,quantile,c(0.025,.975))

# N posterior
N.comp.save <- rep(0, n.mcmc - n.burn)

for(k in 1:(n.mcmc - n.burn)){
  if(k%%10000==0){cat(k," ")}
  beta.0.tmp=beta.0.save.full.lik[k]
  beta.tmp=beta.save.full.lik[,k]
  lam.nowin.int=sum(exp(log(ds)+beta.0.tmp+X.nowin.full%*%beta.tmp)) # can parallelize
  N.comp.save[k]=n+rpois(1,lam.nowin.int)
};cat("\n")

hist(N.comp.save)

# posterior summary
mean(N.comp.save)
sd(N.comp.save)
quantile(N.comp.save, c(0.025, 0.975))


# --- Fit comp. likelihood w/ ESN ----------------------------------------------
source(here("GlacierBay_Code", "spp.comp.ESN.mcmc.R"))
q <- 10
beta.tune <- 0.1
tic()
out.comp.esn <- spp.comp.ESN.mcmc(seal.mat, X.full, win.idx, seal.idx,
                               q, ds, n.mcmc, beta.tune)
toc()


# --- Fit SPP w/ cond. likelihood (num quad stage 1) ---------------------------
source(here("GlacierBay_Code", "spp_win_2D", "spp.cond.mcmc.R"))
tic()
out.cond.full = spp.cond.mcmc(seal.mat,X.obs,X.win.full,ds,n.mcmc)
toc() # 290.419 sec (~4.7 min)

# discard burn-in
beta.save <- out.cond.full$beta.save[,-(1:n.burn)]
beta.0.save <- out.cond.full$beta.0.save[-(1:n.burn)]

# trace plots
layout(matrix(1:2,2,1))
# plot(out.cond.full$beta.0.save,type="l")
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

# --- Sample beta_0 using num quad stage 1 samples -----------------------------
beta.save <- out.cond.full$beta.save
theta.save <- rep(0,n.mcmc)

for(k in 1:n.mcmc){
  if(k%%1000==0){cat(k," ")}
  lam.int <- sum(exp(log(ds)+X.full%*%beta.save[,k]))
  theta.save[k] <- rgamma(1, 0.01 + n, 0.01 + lam.int)
};cat("\n")

beta.0.save <- log(theta.save)

plot(beta.0.save, type ="l")

# --- Fit SPP using cond. likelihood (stan glm stage 1) -----------------------
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

# prepare for second stage
out.cond.bern <- list(beta.save = t(as.matrix(out.bern.cond)[,-1]), 
                      # beta.0.save = out.cond.full$beta.0.save[1:50000],
                      n.mcmc = 50000, n = n, ds = ds, X.full = X.win.full)

# --- Fit SPP using cond. likelihood (Polya-gamma stage 1) ---------------------
X.pg <- cbind(rep(1, nrow(X.obs)), X.obs)

source(here("GlacierBay_Code_OLD", "Polya_Gamma.R"))
p <- ncol(X.pg)
mu.beta <- rep(0, p)
sigma.beta <- diag(2.25, p)
tic()
out.pg <- polya_gamma(y.binary, X.pg, mu.beta, sigma.beta, 100000)
toc() # ~ 18 min
  
# posterior summary
beta.save.pg <- out.pg$beta[,-(1:n.burn)] # discard burn-in

plot(beta.save.pg[2,], type = "l")
plot(beta.save.pg[3,], type = "l")

# compare polya-gamma and full-conditional trace plots
pdf("beta1_trace_compare.pdf")
par(mfrow = c(2,1))
plot(beta.save.full.lik[1,], type = "l")
plot(beta.save.pg[1,], type = "l")
dev.off()

pdf("beta2_trace_compare.pdf")
par(mfrow = c(2,1))
plot(beta.save.full.lik[2,], type = "l")
plot(beta.save.pg[2,], type = "l")
dev.off()

# # prepare for second stage
# beta.0.precise <- rnorm(n.mcmc, mean(out.comp.full$beta.0.save), 
#                                 sd(out.comp.full$beta.0.save))
# out.cond.pg <- list(beta.save = beta.save.pg$beta[-1,], 
#                     # beta.0.save = log(rgamma(n.mcmc, n + 0.001, 1.001)),
#                     n.mcmc = n.mcmc, n = n, ds = ds, X.full = X.win.full)

# x <- seq(0,100,length.out = 10000)
# log_gamma <- function(x){
#   return(log(dgamma(x, 1, 0.001)))
# }
# plot(x = x, y = log_gamma(x), type = "l")
# 
# plot(x = x, y = dgamma(x,1000,1), type = "l")

# --- 2nd stage - compute lambda integrals -------------------------------------
# theta.save <- rep(0,n.mcmc)
lam.int.save <- rep(0, n.mcmc - n.burn)

tic()
for(k in 1:(n.mcmc - n.burn)){
  if(k%%1000==0){cat(k," ")}
  lam.int.save[k] <- sum(exp(log(ds)+X.win.full%*%beta.save.pg[-1,k]))
  # theta.save[k] <- rgamma(1, 0.01 + n, 0.01 + lam.int)
};cat("\n")
toc()

# beta.0.save.gibbs <- log(theta.save)
# 
# plot(beta.0.save, type ="l")

# prepare for third stage
out.cond.pg <- list(beta.save = beta.save.pg[-1,], 
                    # beta.0.save = log(rgamma(n.mcmc, n + 0.001, 1.001)),
                    n.mcmc = n.mcmc - n.burn, n = n, ds = ds, X.full = X.win.full,
                    lam.int.save = lam.int.save)

# --- 3rd stage MCMC -----------------------------------------------------------
source(here("GlacierBay_Code", "spp.stg3.mcmc.R"))
tic()
out.cond.pg3 <- spp.stg3.mcmc(out.cond.pg)
toc() # ~ 1 sec

beta.save <- out.cond.pg3$beta.save
beta.0.save <- out.cond.pg3$beta.0.save

# trace plots
layout(matrix(1:2,2,1))
plot(beta.0.save,type="l")
matplot(t(beta.save),lty=1,type="l")

# posterior summary
beta.save.full <- t(rbind(beta.0.save, beta.save))
vioplot(data.frame(beta.save.full),
        names=expression(beta[0],beta[1],beta[2]),
        ylim = c(-10,5))
abline(h = 0, lty = 2)

apply(beta.save.full,2,mean) 
apply(beta.save.full,2,sd) 
apply(beta.save.full,2,quantile,c(0.025,.975))

# # compare marginal posterior
# hist(out.comp.full$beta.0.save[-(1:1000)],prob=TRUE,breaks=60,main="",xlab=bquote(beta[0]),
#      ylim = c(0,1))
# lines(density(beta.0.save.gibbs[-(1:1000)],n=1000, adjust = 3),col="red",lwd=2)
# lines(density(out.cond.2.full.pg$beta.0.save[-(1:1000)],n=1000, adjust = 3),
#       col="green",lwd=2)

# --- 1st Stage MCMC Plot Comparison -------------------------------------------
beta.save.full.stan1 <- cbind(as.matrix(out.bern.cond), rep(0, nrow(as.matrix(out.bern.cond))))[,-1]
beta.save.full.cond1 <- cbind(beta.save.full, rep(1, nrow(beta.save.full)))[,-1]
beta.save.full.pg <- cbind(t(beta.save), rep(2, nrow(beta.save)))[,-1]
beta.save.full.stage1 <- as.data.frame(rbind( beta.save.full.stan1, beta.save.full.cond1, beta.save.full.pg))

beta.save.stage1.long <- melt(beta.save.full.stage1, id.vars = "V3")

pdf("stage1_compare.pdf")
ggplot(beta.save.stage1.long, aes(x = variable, y = value, fill = as.factor(V3))) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(x = "Coefficients",
       y = "Values",
       fill = "Model") + 
  scale_x_discrete(labels = c("bathymetry", "glacier distance")) + 
  scale_fill_manual(values = c("#00BFC4", "#F8766D", "#7CAE00"),
                    labels = c("stan_glm", "num quad", "polya-gamma"))
dev.off()

# --- Fit SPP using cond. output (num quad stage 2) ----------------------------
source(here("GlacierBay_Code", "spp_win_2D", "spp.stg2.mcmc.R"))

# using num quad results
tic()
out.cond.2.full=spp.stg2.mcmc(out.cond.full)
toc() # 138.445 sec (~2.3 min)

# discard burn-in
beta.save <- out.cond.2.full$beta.save[,-(1:n.burn)]
beta.0.save <- out.cond.2.full$beta.0.save[-(1:n.burn)]

# trace plots
layout(matrix(1:2,2,1))
plot(beta.0.save,type="l")
matplot(t(beta.save),lty=1,type="l")

# posterior summary
beta.save.full <- t(rbind(beta.0.save, beta.save))
vioplot(data.frame(beta.save.full),
        names=expression(beta[0],beta[1],beta[2]),
        ylim = c(-10,5))
abline(h = 0, lty = 2)

apply(beta.save.full,2,mean) 
apply(beta.save.full,2,sd) 
apply(beta.save.full,2,quantile,c(0.025,.975))

# --- Fit SPP using cond. output (stan glm stage 2) ----------------------------
tic()
out.cond.2.bern <- spp.stg2.mcmc(out.cond.bern)
toc() 

# discard burn-in
beta.save <- out.cond.2.bern$beta.save[,-(1:n.burn)]
beta.0.save <- out.cond.2.bern$beta.0.save[-(1:n.burn)]

# trace plots
layout(matrix(1:2,2,1))
plot(beta.0.save,type="l")
matplot(t(beta.save),lty=1,type="l")

# posterior summary
beta.save.full <- t(rbind(beta.0.save, beta.save))
vioplot(data.frame(beta.save.full),
        names=expression(beta[0],beta[1],beta[2]),
        ylim = c(-10,5))
abline(h = 0, lty = 2)

apply(beta.save.full,2,mean) 
apply(beta.save.full,2,sd) 
apply(beta.save.full,2,quantile,c(0.025,.975))

# --- Fit SPP using cond. output (polya-gamma stage 2) -------------------------
tic()
out.cond.2.full.pg <- spp.stg2.mcmc(out.cond.pg)
toc() # 128.9 sec (~2 min)

# discard burn-in
beta.save <- out.cond.2.full.pg$beta.save[,-(1:n.burn)]
beta.0.save <- out.cond.2.full.pg$beta.0.save[-(1:n.burn)]

# trace plots
layout(matrix(1:2,2,1))
plot(beta.0.save,type="l")
matplot(t(beta.save),lty=1,type="l")

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
# complete likelihood samples
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


# num quad samples
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


# polya-gamma samples
N.save.pg=rep(0,n.mcmc)

tic()
for(k in 1:n.mcmc){
  if(k%%10000==0){cat(k," ")}
  beta.0.tmp=out.cond.2.full.pg$beta.0.save[k]
  beta.tmp=out.cond.2.full.pg$beta.save[,k]
  lam.nowin.int=sum(exp(log(ds)+beta.0.tmp+X.nowin.full%*%beta.tmp)) # can parallelize
  N.save.pg[k]=n+rpois(1,lam.nowin.int)
};cat("\n")
toc() # ~ 2min

par(mfrow = c(1,1))

# discard burn-in
N.save.pg <- N.save.pg[-(1:n.burn)]

plot(N.save.pg,type="l")
hist(N.save.pg,breaks=50,prob=TRUE,main="",xlab="N")

# posterior summary
mean(N.save.pg)
sd(N.save.pg)
quantile(N.save.pg, c(0.025, 0.975))


# --- Posterior Intensity Function ---------------------------------------------
beta.save <- beta.save.full.lik
beta.0.save <- beta.0.save.full.lik
beta.save.full <- cbind(beta.0.save, t(beta.save))

# posterior mean heat map
beta.post.means <- apply(beta.save.full,2,mean)
lam.full <- exp(beta.post.means[1] + X.full%*%beta.post.means[-1])
lam.full.df <- as.data.frame(cbind(full.coord, lam.full))
lam.full.rast <- rasterFromXYZ(lam.full.df)

# # get cell with highest intensity
# lam.max.s <- lam.full.df[which(lam.full.df[,3] == max(lam.full.df[,3])),][-3]

pdf("posterior_mean_heatmap.pdf")
plot(lam.full.rast, col = viridis(100))
# plot(survey.win, add = TRUE)
dev.off()


# --- Simulating seal realizations ---------------------------------------------
# # get non-windowed cells
# X.nonwin <- X.full[-win.idx,]
# lam.nonwin <- exp(beta.post.means[1] + X.nonwin%*%beta.post.means[-1])
# 
# # create unobserved window for (S_0)
# survey.poly.nonwin <- st_difference(survey.poly, footprint)
# 
# S0.windows <- list()
# for(i in 1:length(survey.poly.nonwin)){
#   S0.mat <- survey.poly.nonwin[[i]][[1]]
#   if(class(survey.poly.nonwin[[i]])[2] == "GEOMETRYCOLLECTION"){
#     S0.mat <- st_collection_extract(survey.poly.nonwin[[i]], "POLYGON")[[1]]
#   }
#   if(class(S0.mat)[1] == "list"){
#     S0.mat <- S0.mat[[1]]
#   }
#   S0.windows[[i]] <- owin(poly = data.frame(x=S0.mat[,1],
#                            y=S0.mat[,2]))
# }
# S0.win <- do.call(union.owin, S0.windows)
# 
# # simulate superpop
# M=rpois(1, max(lam.nonwin)) 
# s.superpop <- rpoint(100, win = S0.win)
# 
# # prepare X matrix
# superpop.full.idx <- cellFromXY(bath.rast.survey, seal.mat)
# row.counts <- table(factor(seal.full.idx, levels = 1:length(bath.rast.survey)))
# bath.full <- cbind(values(bath.rast.survey), row.counts)
# bath <- na.omit(bath.full)
# X.full <- cbind(bath, glac.dist)
# seal.idx <- c()
# for(i in 1:nrow(X.full)){
#   if(X.full[i,2] != 0){
#     seal.idx <- c(seal.idx, rep(i, times = X.full[i,2]))
#   }
# }
# X.full <- scale(X.full[,-2])

lam.max <- max(lam.full)
M <- rpois(1, area.owin(survey.win)*lam.max)
s.superpop.full <- rpoint(M, win = survey.win)

superpop.nonwin <- !inside.owin(s.superpop.full$x, s.superpop.full$y, footprint.win)
s.superpop.nonwin <- cbind(x = s.superpop.full$x, s.superpop.full$y)[which(superpop.nonwin == TRUE),]
M0 <- nrow(s.superpop.nonwin)

# prepare X matrix
superpop.nonwin.idx.full <- cellFromXY(ice.rast, s.superpop.nonwin)
row.counts <- table(factor(superpop.nonwin.idx.full, levels = 1:length(ice.rast)))
ice.full <- cbind(values(ice.rast), row.counts)
ice <- na.omit(ice.full)
X.superpop.full <- cbind(ice, bath, glac.dist)
X.superpop.full <- na.omit(X.superpop.full)
superpop.nonwin.idx <- rep(seq_len(nrow(X.full)), times = X.superpop.full[, 2])
X.superpop.full <- scale(cbind(X.superpop.full[,-2], full.coord))
X.superpop.nonwin <- X.superpop.full[superpop.nonwin.idx,]

gelu <- function(z){	
  z*pnorm(z)
}

A <- out.comp.esn$A
W.superpop.nonwin <- gelu(X.superpop.nonwin%*%A)
M0 <- nrow(X.superpop.nonwin)

# thin superpop
lam.superpop.nonwin=exp(beta.post.means[1] + W.superpop.nonwin%*%beta.post.means[-1])
# lam.superpop.nonwin <- values(lam.full.rast)[superpop.nonwin.idx]

superpop.nonwin.idx <- cellFromXY(lam.full.rast, s.superpop.nonwin)
lam.superpop.nonwin <- values(lam.full.rast)[superpop.nonwin.idx]
lam.superpop.nonwin.mat <- na.omit(cbind(s.superpop.nonwin, lam.superpop.nonwin))
lam.superpop.nonwin.df <- as.data.frame(lam.superpop.nonwin.mat)
colnames(lam.superpop.nonwin.df) <- c("x", "y", "fill")

obs.idx=rbinom(M0,1,lam.superpop.nonwin.mat[,3]/lam.max)==1
s.obs=lam.superpop.nonwin.mat[obs.idx,1:2] # total observed points 
lam.obs <- lam.superpop.nonwin.mat[obs.idx,3]
N0.pred=nrow(s.obs)

# check how many seals on 0 ice
seal.ice.idx <- cellFromXY(ice.rast, s.obs)
num.seal.0.ice <- sum(na.omit(values(ice.rast)[seal.ice.idx] == 0))

pdf("simulate_08132007_2.pdf")
ggplot() + 
  geom_sf(data = survey.poly) +
  # geom_tile(data = lam.full.df, aes(x = x, y = y, 
  #  fill = V3), color = NA) + 
  geom_sf(data = footprint) + 
  labs(color = "lambda") +
  geom_point(aes(x = s.obs[,1], y = s.obs[,2], color = lam.obs),
             size = 0.2) +
  geom_sf(data = seal.locs, size = 0.2, color = "red") +
  theme(axis.title = element_blank())
dev.off()
# + 
# geom_sf(data = seal.locs, size = 0.1, color = "red")

