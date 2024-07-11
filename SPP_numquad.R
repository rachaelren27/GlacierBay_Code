library(sf)
library(here)
library(tidyverse)
library(terra)
library(tidyterra)
library(raster)

# read in NPS data
path <- here("NPS_data", "HARBORSEAL_2007", "seal_locations_final",
             "nonpup_locs")
seal.locs.20070813 <- st_read(dsn = path, layer = "JHI_20070813_nonpup_locs")

path <- here("NPS_data", "HARBORSEAL_2007", "footprints")
footprint.20070813 <- st_read(dsn = path, layer = "JHI_20070813_footprint")

path <- here("NPS_data", "HARBORSEAL_2007", "survey_polygons")
survey.poly.20070813 <- st_read(dsn = path, layer = "JHI_20070813_surveyboundary")

ggplot() + 
  geom_sf(data = survey.poly.20070813) + 
  geom_sf(data = footprint.20070813) + 
  geom_sf(data = seal.locs.20070813, size = 0.5)

# convert CRS
survey.poly <- st_transform(survey.poly.20070813$geometry, 
                            CRS("+proj=longlat +datum=WGS84"))
survey.poly.mat <- survey.poly[[1]][[1]]
survey.win <- owin(poly = data.frame(x=rev(survey.poly.mat[,1]),
                                     y=rev(survey.poly.mat[,2])))

seal.locs <- st_transform(seal.locs.20070813$geometry,
                          CRS("+proj=longlat +datum=WGS84"))
seal.coord <- st_coordinates(seal.locs)

seal.X <- seal.coord[,1]
seal.Y <- seal.coord[,2]
# all.seal.20070618 <- ppp(x = seal.X, y = seal.Y,
#                          owin = survey.win) # warning

# check if seal locs lie within survey boundary 
grid.tf <- inside.owin(x = seal.X, y = seal.Y, survey.win) 
sum(grid.tf) # all true!

# convert CRS
footprint <- st_transform(footprint.20070813$geometry, 
                          CRS("+proj=longlat +datum=WGS84"))

footprint.coord <- matrix(data = footprint[[1]][[1]], nrow = 5, ncol = 2)
for(i in 2:length(footprint)){
  footprint.coord <- rbind(footprint.coord, footprint[[i]][[1]])
}

footprint.win <- owin(poly = data.frame(x=rev(footprint.coord[,1]),
                                        y=rev(footprint.coord[,2])))

# check if seal locs lie within footprints
grid.tf <- inside.owin(x = seal.X, y = seal.Y, footprint.win)
sum(grid.tf) # all true!
# grid.f.ind <- which(grid.tf == FALSE)

# ggplot() +
#   geom_sf(data = survey.poly) +
#   geom_point(data = as.data.frame(seal.coord[grid.f.ind,]), aes(x = X, y = Y), 
#              size = 0.5)

# read in covariates
bath.rast <- raster(here("covariates", "bathymetry.tiff"))

# crop using survey polygon
bath.rast.survey <- raster::crop(bath.rast, extent(survey.poly.mat))
bath.rast.survey <- raster::mask(bath.rast.survey, as(survey.poly, 'Spatial'))

plot(bath.rast)
plot(survey.poly, add = TRUE)

pdf("survey_crop.pdf")
plot(bath.rast.survey)
plot(survey.poly, add = TRUE)
points(x = seal.X, y = seal.Y)
dev.off()

# crop using footprint
bath.rast.footprint <- raster::crop(bath.rast, extent(footprint.coord))
bath.rast.footprint <- raster::mask(bath.rast.footprint, as(footprint, 'Spatial'))

pdf("footprint_crop.pdf")
plot(bath.rast.footprint)
plot(footprint, add = TRUE)
dev.off()

pdf("footprint_crop2.pdf")
plot(bath.rast.footprint)
dev.off()

# prepare covariates and observed indices
X.grid <- scale(values(bath.rast.footprint))
p <- 1
obs <- cellFromXY(bath.rast.footprint, seal.coord)
row.counts <- table(factor(obs, levels = 1:length(X.grid)))
X.grid <- cbind(X.grid, row.counts)
X.grid <- na.omit(X.grid)
obs <- c()
for(i in 1:nrow(X.grid)){
  if(X.grid[i,p+1] != 0){
    obs <- c(obs, rep(i, times = X.grid[i,p+1]))
  }
}
X.grid <- X.grid[,-(p+1)]

length(obs) # 577 (some locs are mapped to NA in bathymetry raster)

# check how many locs mapped to NA for full bathymetry raster
obs.idx.full <- cellFromXY(bath.rast, seal.coord)
sum(is.na(values(bath.rast)[obs.idx.full])) # 0

# why length of footprint raster > length of survey raster
length(bath.rast.footprint) # 1654497
length(bath.rast.survey) # 1619618