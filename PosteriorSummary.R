library(raster)
library(sf)
library(tidyverse)
library(spatstat)

# --- Calculate intensity (lambda) using beta posterior samples ----------------
n.mcmc <- nrow(beta.post)
beta.post.means <- apply(beta.post,2,mean)

calc_lambda <- function(X, beta){
  beta.0 <- beta[1]
  beta <- as.matrix(beta[-1])
  return(exp(beta.0 + X%*%beta))
}

lam.post.mean <- calc_lambda(X.full, beta.post.means)

# for plotting
s.full <- xyFromCell(bath.rast.survey, which(!is.na(values(bath.rast.survey))))
colnames(s.full) <- c("lat", "long")
lam.post.mean.df <- as.data.frame(cbind(s.full, lam.post.mean))
colnames(lam.post.mean.df) <- c("lat", "long", "lam.post.mean")
lam.post.mean.rast <- rasterFromXYZ(lam.post.mean.df)


# --- Simulate point realizations from posterior intensity ---------------------
sim_points <- function(lam, n, survey.win, footprint.win, nonwin = TRUE){
  lam.max <- max(lam)
  M <- rpois(1, area.owin(survey.win)*lam.max)
  superpop.full <- rpoint(M, win = survey.win)
  
  if(nonwin){
    is.superpop.nonwin <- !inside.owin(superpop.full$x, superpop.full$y, footprint.win)
    superpop <- cbind(x = superpop.full$x, superpop.full$y)[which(is.superpop.nonwin == TRUE),]
    
  } else{
    superpop <- cbind(superpop.full$x, superpop.full$y)
  }
  
  superpop.idx <- cellFromXY(lam.post.mean.rast, superpop)
  lam.superpop <- values(lam.post.mean.rast)[superpop.idx]
  lam.superpop.mat <- na.omit(cbind(superpop, lam.superpop))
  M <- nrow(lam.superpop.mat)

  obs.idx <- rbinom(M,1,lam.superpop.mat[,3]/lam.max)==1
  s.obs <- lam.superpop.mat[obs.idx,1:2] # total observed points 
  lam.obs <- lam.superpop.mat[obs.idx,3]
  
  if(nonwin){
    N0.pred <- nrow(s.obs)
    N.pred <- N0.pred + n
  } else{
    N.pred <- nrow(s.obs)
  }
  
  return(list(s.obs, lam.obs, N.pred))
}

post.mean.sim <- sim_points(lam.post.mean, length(seal.locs), survey.win,
                            footprint.win, nonwin = FALSE)
post.mean.points <- post.mean.sim[[1]]
post.mean.lam <- post.mean.sim[[2]]


# --- Posterior for N ----------------------------------------------------------
N.save=rep(0,n.mcmc)

ds <- res(bath.rast.survey)[1]*res(bath.rast.survey)[2]

for(k in 1:n.mcmc){
  if(k%%10000 == 0){cat(k," ")}
  beta.0.tmp <- beta.post[1,k]
  beta.tmp <- beta.post[-1,k]
  lam.nonwin.int <- sum(exp(log(ds)+beta.0.tmp+X.nonwin.full%*%beta.tmp)) # numerical quadrature to approximate integral
  N.save[k] <- n + rpois(1,lam.nonwin.int)
};cat("\n")

hist(N.save, breaks = 50, prob = TRUE, main = "", xlab = "N")
abline(v = N.pred, lty = 2,lwd = 2, col = "red")


# --- Plotting -----------------------------------------------------------------
# Plot intensity function
ggplot() + 
  geom_tile(data = lam.post.mean.df, aes(x = lat, y = long,
   fill = lam.post.mean))

# Plot seal data
ggplot() + 
  geom_sf(data = survey.poly) + 
  geom_sf(data = footprint) + 
  geom_sf(data = seal.locs, size = 0.1)

# Plot bathymetry
plot(bath.rast)

# Plot glacier distance
plot(glac.dist.rast)

# Plot point realization
ggplot() + 
  geom_sf(data = survey.poly) +
  geom_sf(data = footprint) + 
  geom_point(aes(x = post.mean.points[,1], y = post.mean.points[,2],
                 color = post.mean.lam), size = 0.1) +
  # geom_sf(data = seal.locs, size = 0.1, color = "red") +
  theme(axis.title = element_blank())


# --- Spatial Count Plot ------------------------------------------------------
survey.poly.sfc <- st_sfc(st_polygon(list(survey.poly[[1]][[1]])), crs = "WGS84")
survey.vect <- vect(survey.poly)
coarse.rast <- rast(ext(survey.vect), resolution = 0.00125)
values(coarse.rast) <- rep(0, dim(coarse.rast)[1]*dim(coarse.rast)[2], crs = "WGS84")
coarse.rast <- mask(coarse.rast, survey.vect)

# s.sim <- post.mean.sim[[1]]
# cell.counts <- table(cellFromXY(coarse.rast, xy = s.sim))
# 
s.full <- xyFromCell(coarse.rast, cell = 1:(dim(coarse.rast)[1]*dim(coarse.rast)[2]))
# count.mat <- cbind(s.full, z = values(coarse.rast))
# count.mat[as.integer(names(cell.counts)),3] <- cell.counts
# count.rast <- rasterFromXYZ(count.mat)
# 
# plot(count.rast, col = viridis(100))
# points(x = seal.mat[,1], y = seal.mat[,2], col = "red", pch = 19, cex = 0.1)

# # crop survey boundary
# count.rast <- terra::rast(count.rast)
# count.rast <- mask(count.rast, survey.vect)
# plot(count.rast, col = viridis(100), ylim = c(58.82, 58.92), xlim = c(-137.15, -137))

# compute count posterior mean and variance
s.sim.list <- list()
count.mat <- matrix(0, nrow = (n.mcmc - n.burn), ncol = length(values(coarse.rast)))
for(i in 1:(n.mcmc - n.burn)){
  lambda.full <- calc_lambda(W.full, beta.save.full[i,])
  s.sim.list[[i]] <- sim_points(lambda.full, n, survey.win,
                                footprint.win, nonwin = FALSE)[[1]]
  cell.counts <- table(cellFromXY(coarse.rast, xy = s.sim))
  count.mat[i, as.integer(names(cell.counts))] <- cell.counts
  if(i %% 100 == 0){
    print(i)
  }
}

count.means <- apply(count.mat, 2, mean)
count.mat.mean <- cbind(s.full, z = count.means)
count.rast.mean <- rasterFromXYZ(count.mat.mean)
count.rast.mean <- terra::rast(count.rast.mean)
count.rast.mean <- mask(count.rast.mean, survey.vect)
plot(count.rast.mean)

rast.df <- as.data.frame(count.rast.mean, xy = TRUE, na.rm = TRUE)

pdf("posterior_mean_count_20070621.pdf", compress = FALSE)
ggplot(data = as.data.frame(rast.df)) + 
  geom_tile(aes(x = x, y = y, fill = z, col = z)) +
  scale_color_viridis_c(guide = "none") +
  scale_fill_viridis_c(name = "count") +
  theme(
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),  
    axis.ticks = element_line(),
  ) +
  theme_bw() + 
  xlim(c(-137.14, -137)) + 
  theme(axis.title = element_blank())
dev.off()
  
count.med <- apply(count.mat, 2, median)
count.mat.med <- cbind(s.full, z = count.med)
count.rast.med <- rasterFromXYZ(count.mat.med)

count.var <- apply(count.mat, 2, sd)
count.mat.var <- cbind(s.full, z = count.var)
count.rast.var <- rasterFromXYZ(count.mat.var)
count.rast.var <- terra::rast(count.rast.var)
count.rast.var <- mask(count.rast.var, survey.vect.crop)

count.sum <- apply(count.mat, 1, sum)
pdf("N_post_pred.pdf")
hist(count.sum, main = NULL, xlab = "Total Abundance", breaks = 30, prob = TRUE,
     ylab = "Probability")
dev.off()
