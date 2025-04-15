setwd("/Users/rlr3795/Desktop/GlacierBay_Project")

library(sf)
library(here)
library(tidyverse)
library(terra)
library(raster)
library(tictoc)
library(spatstat)
library(rstanarm)
library(BayesLogit)
library(mvnfast)
library(pgdraw)
library(coda)
library(parallel)
library(viridis)
library(foreach)
library(doParallel)

load(here("SPP_script.RData"))
set.seed(1234)

# --- Read in NPS data ---------------------------------------------------------
year <- "2007"
date <- "20070618"

path <- here("NPS_data", paste0("HARBORSEAL_", year), "seal_locations_final",
             "pup_locs")
seal.locs <- st_read(dsn = path, layer = paste0("JHI_", date, "_pup_locs"))

path <- here("NPS_data", paste0("HARBORSEAL_", year), "footprints")
footprint <- st_read(dsn = path, layer =  paste0("JHI_", date, "_footprint"))

path <- here("NPS_data", paste0("HARBORSEAL_", year), "survey_polygons")
survey.poly <- st_read(dsn = path, layer = paste0("JHI_", date, "_surveyboundary"))

# convert CRS
survey.poly <- st_transform(survey.poly$geometry, 
                            CRS("+proj=longlat +datum=WGS84"))

seal.locs <- st_transform(seal.locs$geometry,
                          CRS("+proj=longlat +datum=WGS84"))
seal.mat <- as.matrix(st_coordinates(seal.locs))

# plot seals
ggplot() +
  geom_sf(data = survey.poly) + 
  geom_point(aes(x = seal.mat[,1], y = seal.mat[,2]), size = 0.5) + 
  labs(title = date)

footprint <- st_transform(footprint$geometry, 
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


# --- Read in covariates -------------------------------------------------------
# read in bathymetry
bath.rast <- raster(here("covariates", "bathymetry.tiff"))

ice.rast <- raster(here("covariates", "LK_ice_estimates.tiff"))
ice.rast <- raster::crop(ice.rast, extent(survey.poly.mat))
ice.rast <- raster::mask(ice.rast, as(survey.poly, 'Spatial'))
plot(ice.rast)

# ice.df <- as.data.frame(cbind(xyFromCell(ice.rast, 1:length(ice.rast)),
#                               values(ice.rast)))
# 
# ggplot() + 
#   geom_sf(data = survey.poly) + 
#   geom_sf(data = footprint) + 
#   geom_tile(data = ice.df, aes(x = x, y = y,
#                                fill = is.na(V3))) + 
#   scale_fill_manual(
#     values = c("TRUE" = "blue", "FALSE" = "red"),
#     name = "Ice NA") # + 
#   # geom_sf(data = seal.locs, size = 0.1, col = "white") + 
#   # coord_sf(xlim = c(NA, -137.06), ylim = c(NA, 58.86))

glac.dist.rast <- raster(here("covariates", "glacier_dist.tiff"))

## crop using survey boundary
bath.rast.survey <- raster::crop(bath.rast, extent(survey.poly.mat))
bath.rast.survey <- raster::mask(bath.rast.survey, as(survey.poly, 'Spatial'))
# # remove NAs
# s.bath.rast <- xyFromCell(bath.rast.survey, which(!is.na(values(bath.rast.survey))))
# bath.rast.df <- data.frame(x = s.bath.rast[,1], y = s.bath.rast[,2],
#                            z = na.omit(values(bath.rast.survey)))
# bath.rast.survey.noNA <- rasterFromXYZ(bath.rast.df)
# extent(bath.rast.survey.noNA) <- extent(survey.poly.mat)

# s.bath.rast.na <- xyFromCell(bath.rast.survey, which(is.na(values(bath.rast.survey))))
# plot(bath.rast.survey)
# plot(survey.poly, add = TRUE)
# points(x = s.bath.rast.na[,1], y = s.bath.rast.na[,2])

# calculate distance from southern boundary (glacier)
# ggplot() + 
#   geom_sf(data = survey.poly) + 
#   geom_point(aes(x = -137.1311, y = 58.84288), color = "red") # westmost point
# 
# ggplot() + 
#   geom_sf(data = survey.poly) +
#   geom_line(data = glacier.poly, aes(x = V1, y = V2), col = "blue") + 
#   geom_point(aes(x = glacier.poly[101,1], y = glacier.poly[101, 2]), color = "red")
# 
# survey.poly.df <- as.data.frame(survey.poly.mat)
# glacier.poly <- survey.poly.df %>% filter(V2 < 58.84288)
# glacier.poly <- as.matrix(glacier.poly[-(1:100),]) # found index 101 using localMinima
# 
# ggplot() + 
#   geom_sf(data = survey.poly) + 
#   geom_point(data = glacier.poly, aes(x = V1, y = V2), color = "red")

# bath.rast <- na.omit(values(bath.rast.survey))

# # calculate glacier distance
# seal.glac.dist <- dist2Line(seal.mat, glacier.poly) # in meters
# 
# bath.survey.idx <- which(!is.na(values(bath.rast.survey)))
# full.coord <- xyFromCell(bath.rast.survey, bath.survey.idx)
# 
# full.glac.dist <- dist2Line(full.coord, glacier.poly) # takes a while
# 
# glac.dist.df <- data.frame(x = s.full[,1], y = s.full[,2],
#                            z = full.glac.dist[,1])
# glac.dist.rast <- rasterFromXYZ(glac.dist.df)
# writeRaster(glac.dist.rast, filename = "glacier_dist.tiff", format = "GTiff")


# --- Calculate areas ----------------------------------------------------------
tot.area <- area.owin(survey.win)
tot.win.area <- area.owin(footprint.win)
tot.nonwin.area <- tot.area - tot.win.area
n.win <- length(footprint)
ex.win <- owin(poly = data.frame(x = footprint[[1]][[1]][,1],
                                 y = footprint[[1]][[1]][,2]))
win.area <- area.owin(ex.win) # approx. bc windows not equally sized

ds <- res(bath.rast.survey)[1]*res(bath.rast.survey)[2]
  

# --- Set X matrices -----------------------------------------------------------
# glac.dist <- full.glac.dist[,1]

# bath.survey.idx <- which(!is.na(values(bath.rast.survey)))
# full.coord <- xyFromCell(bath.rast.survey, bath.survey.idx)

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
X.full <- cbind(ice, bath, glac.dist)
seal.idx <- c()
for(i in 1:nrow(X.full)){
  if(X.full[i,2] != 0){
    seal.idx <- c(seal.idx, rep(i, times = X.full[i,2]))
  }
}
X.full <- na.omit(X.full)
X.full <- scale(X.full[,-2])

win.idx <- which(inside.owin(ice.full.coord[,1], ice.full.coord[,2], footprint.win)) # 15 sec

X.win.full <- X.full[win.idx,]

X.obs <- X.full[seal.idx,]
n <- nrow(X.obs)


# --- Fit SPP w/ Complete Likelihood -------------------------------------------
n.mcmc <- 100000
source(here("GlacierBay_Code", "spp_win_2D", "spp.comp.mcmc.R"))
theta.tune <- 0.1
beta.tune <- 0.01
tic()
out.comp.full <- spp.comp.mcmc(seal.mat, X.obs, X.win.full, ds, n.mcmc, theta.tune,
                               beta.tune)
toc() # 453.14 sec elapsed (~7.5 min)

# discard burn-in
n.burn <- 0.1*n.mcmc
beta.save.full.lik <- out.comp.full$beta.save[,-(1:n.burn)]
beta.0.save.full.lik <- out.comp.full$beta.0.save[-(1:n.burn)]

# trace plots
layout(matrix(1:2,2,1))
plot(beta.0.save.full.lik,type="l")
matplot(t(beta.save.full.lik),lty=1,type="l", col = c("black", "red"))

# # find optimal tuning parameters
# beta.0.tune <- c(0.01,0.05,0.1,0.5)
# beta.tune <- c(0.01,0.05,0.1,0.5)
# tune.params <- expand.grid(beta.0.tune, beta.tune)
# 
# out.comp.compare <- list()
# effective.size <- matrix(nrow = nrow(tune.params), ncol = 3)
# # for(i in 1:nrow(tune.params)){
# #   out.comp <- spp.comp.mcmc(seal.mat,X.obs,X.win.full,ds,win.area,
# #                                          10000,tune.params[i,1], tune.params[i,2])
# #   effective.size[i,] <- c(effectiveSize(out.comp.compare$beta.0.save[-(1:1000)]),
# #                           effectiveSize(out.comp.compare$beta.save[1,-(1:1000)]),
# #                           effectiveSize(out.comp.compare$beta.save[2,-(1:1000)]))
# #   out.comp.compare[[i]] <- out.comp
# # }
# 
# process_function <- function(i) {
#   result <- list()
#   out.comp <- spp.comp.mcmc(seal.mat, X.obs, X.win.full, ds, win.area, 
#                             10000, tune.params[i, 1], tune.params[i, 2])
#   result$out.comp <- out.comp
#   result$eff.size <- c(effectiveSize(out.comp$beta.0.save[-(1:2000)]),
#                        effectiveSize(out.comp$beta.save[1, -(1:2000)]),
#                        effectiveSize(out.comp$beta.save[2, -(1:2000)]))
#   return(result)
# }
# 
# # Run the parallel loop
# results <- mclapply(1:nrow(tune.params), process_function, mc.cores = 10)
# 
# # Combine the results
# out.comp.compare <- lapply(results, function(x) x$out.comp)
# effective.size <- do.call(rbind, lapply(results, function(x) x$eff.size))
# 
# # # posterior summary
# # beta.save.full <- t(rbind(beta.0.save, beta.save))
# # vioplot(data.frame(beta.save.full),
# #         names=expression(beta[0],beta[1],beta[2]),
# #         ylim = c(-10,5))
# # abline(h = 0, lty = 2)

beta.save.full <- t(rbind(beta.0.save.full.lik, beta.save.full.lik))
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

hist(N.comp.save, breaks = 10)

# posterior summary
mean(N.comp.save)
sd(N.comp.save)
quantile(N.comp.save, c(0.025, 0.975))


# --- Fit comp. likelihood w/ ESN ----------------------------------------------
source(here("GlacierBay_Code", "spp.comp.ESN.mcmc.R"))
X.full <- cbind(X.full, full.coord)

theta.tune <- 0.1
beta.tune <- 0.001
q <- 10
lambda <- 1/100
tic()
out.comp.esn <- spp.comp.ESN.mcmc(seal.mat, X.full,
                               win.idx, seal.idx, ds, n.mcmc, theta.tune, 
                               beta.tune, q, lambda)
toc()

matplot(t(out.comp.esn$beta.save), type = 'l')

q <- 12:20
q.lambda <- expand.grid(q,lambda)

cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)

out.comp.esn.list2 <- foreach(k = 1:length(q)) %dopar% {
  q <- q[k]
  lambda <- 1/100
  out.comp.esn <- spp.comp.ESN.mcmc(seal.mat, scale(cbind(X.full, full.coord)),
                                    win.idx, seal.idx, ds, n.mcmc, theta.tune, 
                                    beta.tune, q, lambda)
  
  return(out.comp.esn)
}

stopCluster(cl)

# --- Fit SPP w/ cond. likelihood (num quad stage 1) ---------------------------
source(here("GlacierBay_Code", "spp_win_2D", "spp.cond.mcmc.R"))
tic()
out.cond.full=spp.cond.mcmc(seal.mat,X.obs,X.win.full,ds,n.mcmc)
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

# --- Fit SPP using cond. likelihood (glm stage 1) ------------------------
# obtain background sample
n.bg <- 100000
bg.pts <- rpoint(n.bg, win = footprint.win)
  
# ggplot() + 
#   geom_sf(data = survey.poly) + 
#   geom_sf(data = footprint) + 
#   geom_point(aes(x = bg.pts$x, y = bg.pts$y), size = 0.3) + 
#   geom_sf(data = seal.locs, size = 0.3, col = "red")
   
# prepare covariates for background sample 
bg.mat <- cbind(bg.pts$x, bg.pts$y)
  
bg.full.idx <- cellFromXY(ice.rast, bg.mat)
row.counts <- table(factor(bg.full.idx, levels = 1:length(ice.rast)))
ice.full.counts <- cbind(values(ice.rast), row.counts)
ice <- na.omit(ice.full.counts)
X.temp <- na.omit(cbind(ice, bath, glac.dist))
bg.idx <- c()
for(i in 1:nrow(X.full)){
  if(X.temp[i,2] != 0){
    bg.idx <- c(bg.idx, rep(i, times = X.temp[i,2]))
  }
}

X.obs.aug <- X.full[c(seal.idx, bg.idx),]
y.obs.binary <- rep(0, n + length(bg.idx))
y.obs.binary[1:n] <- 1
logit.obs.df <- data.frame(y = y.obs.binary, ice = X.obs.aug[,1], 
                           bath = X.obs.aug[,2], glac.dist = X.obs.aug[,3])

X.full.aug <- rbind(X.full, X.full[bg.idx,])
X.win.aug <- X.full[c(win.idx, bg.idx),]

# vanilla glm
tic()
out.bern.cond <- glm(y ~ ice + bath + glac.dist, data = logit.obs.df, 
                     family=binomial(link="logit"))
toc() # 0.17
beta.glm <- coef(out.bern.cond)[-1]
se.glm <- summary(out.bern.cond)$coefficients[-1, 2]

# # vanilla Bayesian glm
# tic()
# out.bern.cond <- stan_glm(y ~ ice + bath + glac.dist, data = logit.obs.df, 
#                      family=binomial(link="logit"), iter = 100000, chains = 1)
# toc()


# # vanilla logistic bayesreg
# tic()
# out.bern.cond <- bayesreg(y ~ ice + bath + glac.dist, data = logit.obs.df, 
#                           model = "logistic", n.samples = n.mcmc, burnin = n.burn)
# toc()
# 
# # test ELM logistic glm
# n.sim <- 100
# q <- 5
# A.array <- array(0, c(p, q, n.sim))
# W.array <- array(0, c(nrow(X.obs.aug), q, n.sim))
# aic.vec <- rep(0, n.sim)
# beta.mat <- matrix(0, q+1, n.mcmc)
# 
# for(l in 1:n.sim){
#   A.array[,,l] <- matrix(rnorm(q*p), p, q)
#   W.array[,,l] <- gelu(X.obs.aug%*%A.array[,,l])
#   tmp.lm <- glm(y.obs.binary ~ W.array[,,l], family = binomial(link = "logit"))
#   aic.vec[l] <- AIC(tmp.lm)
#   beta.mat[,l] <- coef(tmp.lm)
# }
# best.idx <- (1:n.sim)[aic.vec == min(aic.vec)]
# cat("best AIC:", aic.vec[best.idx], "\n")
# 
# # ELM logistic bayesreg
# source(here("GlacierBay_Code", "spp.logit.bayesreg.ELM.R"))
# tic()
# out.bern.ELM <- spp.logit.bayesreg.ELM(X.obs.aug, X.full.aug, X.win.aug, win.idx,
#                                        seal.idx, ds, n.mcmc, q)
# toc()

# prepare for second stage
out.cond.bern <- list(beta.save = out.bern.ELM$beta.save, 
                      n.mcmc = n.mcmc, n = n, ds = ds, X.full = X.win.full)


# --- Fit SPP using cond. likelihood (Polya-gamma stage 1) ---------------------
X.pg <- cbind(rep(1, nrow(X.obs)), X.obs)

source(here("GlacierBay_Code", "Polya_Gamma.R"))
p <- ncol(X.pg)
mu.beta <- c(-5, rep(0, p-1))
sigma.beta <- diag(10, p)
# w <- 2^(1-y.binary)
# w <- rep(1, length(y.binary))
tic()
out.pg <- polya_gamma(y.binary, X.pg, mu.beta, sigma.beta, n.mcmc)
toc() # ~ 18 min
# weighted: 380 sec (~6 min)

# posterior summary
beta.save.pg <- out.pg$beta # discard burn-in

par(mfrow = c(2,1))
plot(beta.save.pg[1,], type = "l")
plot(beta.save.pg[3,], type = "l")

# compare polya-gamma and full-conditional trace plots
pdf("beta1_trace_compare.pdf")
par(mfrow = c(2,1))
plot(beta.save.full.lik[1,], type = "l")
plot(beta.save.full.pg[,1], type = "l")
dev.off()

pdf("beta2_trace_compare.pdf")
par(mfrow = c(2,1))
plot(beta.save.full.lik[2,], type = "l")
plot(beta.save.full.pg[,2], type = "l")
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
beta.save <- out.bern.ELM$beta[,-(1:n.burn)]
W.win.full <- out.bern.ELM$W.full[win.idx,]

beta.save <- mvnfast::rmvn(n.mcmc, mu = beta.glm, sigma = diag(se.glm^2))
lam.int.save <- rep(0, n.mcmc)

cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

tic()
lam.int.save <- foreach(k = 1:nrow(beta.save), .combine = c) %dopar% {
  lam.int <- sum(exp(log(ds) + X.win.full %*% beta.save[k,]))
  return(lam.int) 
}
toc()

stopCluster(cl) # 22.4 sec

# prepare for third stage
out.cond.pg2 <- list(beta.save = t(beta.save), n.mcmc = n.mcmc, n = n,
                     ds = ds, X.full = X.win.full, lam.int.save = lam.int.save)

# --- 3rd stage MCMC -----------------------------------------------------------
source(here("GlacierBay_Code", "spp.stg3.mcmc.R"))
tic()
out.cond.pg3 <- spp.stg3.mcmc(out.cond.pg2)
toc() # ~ 1 sec

# discard burn-in
beta.save <- out.cond.pg3$beta.save[,-(1:n.burn)]
beta.0.save <- out.cond.pg3$beta.0.save[-(1:n.burn)]

# trace plots
layout(matrix(1:2,2,1))
plot(beta.0.save,type="l")
plot(beta.save[1,], type = "l")
plot(beta.save[2,], type = "l")

matplot(t(beta.save),lty=1,type="l")

# posterior summary
beta.save.full <- t(rbind(beta.0.save, beta.save))
vioplot(data.frame(beta.save.full),
        names=expression(beta[0],beta[1],beta[2]),
        ylim = c(-10,5))
abline(h = 0, lty = 2)

beta.post.means <- apply(beta.save.full,2,mean) 
beta.post.sd <- apply(beta.save.full,2,sd) 
apply(beta.save.full,2,quantile,c(0.025,.975))

# --- Compare Marginal Posteriors ----------------------------------------------
layout(matrix(1:4,1,4))
# hist(out.comp.full$beta.0.save[-(1:n.burn)], prob=TRUE, breaks=60,main="", 
     xlab=bquote(beta[0]), ylim = c(0,10))
lines(density(out.comp.full$beta.0.save[-(1:n.burn)],n=1000,adj=2),col="red",lwd=2)
lines(density(out.cond.pg3$beta.0.save[-(1:n.burn)], n=1000, adj=2), col="green",
      lwd=2)

hist(out.comp.full$beta.save[1,-(1:n.burn)], prob=TRUE, breaks=60, main="", 
     xlab=bquote(beta[1]), ylim = c(0,15))
# lines(density(out.bern.cond$beta[1,-(1:n.burn)],n=1000,adj=2),col="red",lwd=2)
lines(density(out.cond.pg3$beta.save[1,-(1:n.burn)], n=1000, adj=2), col="green",
      lwd=2)

hist(out.comp.full$beta.save[2,-(1:n.burn)], prob=TRUE, breaks=60, 
     main= "" , xlab=bquote(beta[2]), ylim = c(0,10))
# lines(density(out.bern.cond$beta[2,-(1:n.burn)],n=1000,adj=2),col="red",lwd=2)
lines(density(out.cond.pg3$beta.save[2,-(1:n.burn)], n=1000, adj=2), col="green", 
      lwd=2)

hist(out.comp.full$beta.save[3,-(1:n.burn)], prob=TRUE, breaks=60, main="", 
     xlab=bquote(beta[2]), ylim = c(0,10))
# lines(density(out.bern.cond$beta[3,-(1:n.burn)],n=1000,adj=2),col="red",lwd=2)
lines(density(out.cond.pg3$beta.save[3,-(1:n.burn)], n=1000, adj=2), col="green",
      lwd=2)

# --- Posterior for N ----------------------------------------------------------
# complete likelihood samples
N.save=rep(0,n.mcmc)

X.nowin.full <- X.full[-win.idx,]
n <- nrow(seal.mat)

tic()
for(k in 1:n.mcmc){
  if(k%%10000==0){cat(k," ")}
  beta.0.tmp=out.cond.pg3$beta.0.save[k]
  beta.tmp=out.cond.pg3$beta.save[,k]
  lam.nowin.int=sum(exp(log(ds)+beta.0.tmp+X.nowin.full%*%beta.tmp)) # can parallelize
  N.save[k]=n+rpois(1,lam.nowin.int)
};cat("\n")
toc()

par(mfrow = c(1,1))

# discard burn-in
N.save <- N.save[-(1:n.burn)]

plot(N.save,type="l")
pdf("abundance_posterior_pred.pdf")
hist(N.save,breaks=50,prob=TRUE,main="",xlab="N")
abline(v = 1129, lty = 2, col = "red", lwd = 2)
dev.off()

# posterior summary
mean(N.save)
sd(N.save)
quantile(N.save, c(0.025, 0.975))


# --- Posterior Intensity Function ---------------------------------------------
out.comp.esn <- out.cond.pg3
beta.0.save <- out.comp.esn$beta.0.save
beta.save <- out.comp.esn$beta.save
beta.save.full <- cbind(beta.0.save, t(beta.save))
W.full <- out.bern.ELM$W.full[1:nrow(X.full),]

# posterior mean heat map
beta.post.means <- apply(beta.save.full,2,mean)
lam.full <- exp(beta.post.means[1] + W.full%*%beta.post.means[-1])
lam.full.df <- as.data.frame(cbind(full.coord, lam.full))
lam.full.rast <- rasterFromXYZ(lam.full.df)

# # get cell with highest intensity
# lam.max.s <- lam.full.df[which(lam.full.df[,3] == max(lam.full.df[,3])),][-3]

pdf("posterior_mean_heatmap.pdf")
plot(lam.full.rast, col = viridis(100))
# plot(survey.win, add = TRUE)
dev.off()

# --- Simulating seal realizations ---------------------------------------------
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
X.superpop.full <- scale(X.superpop.full[,-2])
X.superpop.nonwin <- X.superpop.full[superpop.nonwin.idx,]

gelu <- function(z){	
  z*pnorm(z)
}

A <- out.bern.ELM$A
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
             size = 0.5) +
  geom_sf(data = seal.locs, size = 0.5, color = "red") +
  theme(axis.title = element_blank())
dev.off()


# --- L-function p-value -------------------------------------------------------
## compute L-function for observed data
obs.ppp <- ppp(seal.mat[,1], seal.mat[,2], window = footprint.win)
obs.L <- Linhom(obs.ppp)
plot(obs.L)
# lines(x = obs.L$r, y = obs.L$theo)

## compute L-function for 
out.comp.esn <- out.comp.esn.list[[14]]
beta.0.save <- out.comp.esn$beta.0.save[-(1:n.burn)]
beta.save <- out.comp.esn$beta.save[,-(1:n.burn)]
beta.post <- cbind(beta.0.save, t(beta.save))
W.full <- out.comp.esn$W.full

# compute lambda in parallel
cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)

lam.full.list <- foreach(k = 1:100) %dopar% {
  lam.full <- exp(beta.post[k,1] + W.full%*%beta.post[k,-1])
  lam.df <- as.data.frame(cbind(full.coord, lam.full))
  
  return(lam.df)
}

stopCluster(cl)

# simulate points in parallel
sim_points <- function(lam.df, survey.win, footprint.win){
  lam <- lam.df$V3
  lam.max <- max(lam)
  M <- rpois(1, spatstat.geom::area.owin(survey.win)*lam.max)
  superpop.full <- rpoint(M, win = survey.win)
  
  superpop <- cbind(superpop.full$x, superpop.full$y)
  
  is.superpop.win <- inside.owin(superpop.full$x, superpop.full$y, footprint.win)
  superpop <- cbind(x = superpop.full$x, superpop.full$y)[which(is.superpop.win == TRUE),]
  
  lam.rast <- rasterFromXYZ(lam.df)
  superpop.idx <- cellFromXY(lam.rast, superpop)
  lam.superpop <- values(lam.rast)[superpop.idx]
  lam.superpop.mat <- na.omit(cbind(superpop, lam.superpop))
  M <- nrow(lam.superpop.mat)
  
  obs.idx <- rbinom(M,1,lam.superpop.mat[,3]/lam.max)==1 # thin
  s.obs <- lam.superpop.mat[obs.idx,1:2] 
  lam.obs <- lam.superpop.mat[obs.idx,3]
  
  return(list(s.obs, lam.obs))
}

cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)

sim.points.list <- foreach(k = 1:100) %dopar% {
  sim.point <- sim_points(lam.full.list[[k]], survey.win, footprint.win)
  return(sim.point)
}

stopCluster(cl)


