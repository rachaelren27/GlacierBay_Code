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
sim_points_full <- function(lam, full.coord, survey.win){
  lam.max <- max(lam)
  M <- rpois(1, area.owin(survey.win)*lam.max)
  superpop.full <- rpoint(M, win = survey.win)
  
  superpop <- cbind(superpop.full$x, superpop.full$y)
  
  lam.df <- data.frame(x = full.coord[,1], y = full.coord[,2], z = lam)
  lam.rast <- rasterFromXYZ(lam.df)
  superpop.idx <- cellFromXY(lam.rast, superpop)
  lam.superpop <- values(lam.rast)[superpop.idx]
  lam.superpop.mat <- na.omit(cbind(superpop, lam.superpop))
  M <- nrow(lam.superpop.mat)
  
  obs.idx <- rbinom(M,1,lam.superpop.mat[,3]/lam.max)==1
  s.obs <- lam.superpop.mat[obs.idx,1:2] # total observed points 
  lam.obs <- lam.superpop.mat[obs.idx,3]
  
  return(list(s.obs, lam.obs))
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

s.full <- xyFromCell(coarse.rast, cell = 1:(dim(coarse.rast)[1]*dim(coarse.rast)[2]))

# compute count posterior mean and variance
s.sim.list <- list()
count.mat <- matrix(0, nrow = n.mcmc, ncol = length(values(coarse.rast)))
for(i in 1:n.mcmc){
  lambda.full <- calc_lambda(W.full, beta.save.full[i,])
  s.sim.list[[i]] <- sim_points_full(lambda.full, ice.full.coord, survey.win)[[1]]
  cell.counts <- table(cellFromXY(coarse.rast, xy = s.sim.list[[i]]))
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

pdf("post_count_20070621.pdf", compress = FALSE)
ggplot() + 
  geom_tile(data = as.data.frame(rast.df), aes(x = x, y = y, fill = z, col = z)) +
  scale_color_viridis_c(guide = "none") +
  scale_fill_viridis_c(name = "count") +
  # geom_point(aes(x = seal.mat[,1], y = seal.mat[,2]), size = 0.75, col = "red") + 
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

# --- L-function Model Checking ------------------------------------------------
# compute L
cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)

L.fun.list <- foreach(k = 1:n.mcmc, .packages = c('spatstat', 'raster')) %dopar% {
  sim.mat <- s.sim.list[[k]]
  sim.ppp <- ppp(sim.mat[,1], sim.mat[,2], footprint.win)
  sim.L <- Linhom(sim.ppp)
  return(cbind(sim.L$r, sim.L$iso))
}

stopCluster(cl)

# plot L functions
L.df <- lapply(1:1000, function(i) {
  df <- as.data.frame(L.fun.list[[i]])
  names(df) <- c("r", "Lr")
  df$iso <- df$Lr - df$r
  df$sim_id <- i
  df
}) %>% bind_rows()

obs.df <- data.frame(r = obs.L$r,
                     iso = obs.L$iso - obs.L$r)

baseline.df <- data.frame(r = obs.L$r,
                          iso = 0)

L.fun.plot <- ggplot() +
  geom_line(data = L.df, aes(x = r, y = iso, group = sim_id), color = "gray60", alpha = 0.6) +
  geom_line(data = baseline.df, aes(x = r, y = iso), color = "black", linetype = "dashed") +
  geom_line(data = obs.df, aes(x = r, y = iso), color = "red", linewidth = 1) +
  labs(x = "r", y = "L(r) - r") +
  theme_classic() + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))

# Bayesian p-value
sim.iso.mat <- matrix(NA, nrow = length(obs.L$r), ncol = length(L.fun.list))

for (i in 1:1000) {
  sim.iso.mat[,i] <- L.fun.list[[i]][,2]
}

bayes.p <- rowSums(sim.iso.mat < obs.L$iso)/1000

plot(bayes.p)