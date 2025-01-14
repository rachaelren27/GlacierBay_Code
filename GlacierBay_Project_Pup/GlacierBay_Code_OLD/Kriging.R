library(sf)
library(sp)
library(here)
library(tidyverse)
library(stringi)
library(spatstat)
library(stelfi)
library(INLA)
library(raster)
library(tidyr)

path <- here("Data", "HARBORSEAL_2007", "seal_locations_final",
             "nonpup_locs")
np.seal.locs.20070618 <- st_read(dsn = path, layer = "JHI_20070618_nonpup_locs")

# get coordinates
geom.points <- st_coordinates(np.seal.locs.20070618$geometry) # these are not coordinates
np.seal.20070618.X <- geom.points[,1]
np.seal.20070618.Y <- geom.points[,2]

ggplot() +
  geom_sf(data = survey.poly.20070618) + 
  geom_sf(data = np.seal.locs.20070618) + 

# fit lgcp
survey.poly.mat <- survey.poly.20070618$geometry[[1]][[1]]
survey.win <- owin(poly = data.frame(x=rev(survey.poly.mat[,1]),
                                     y=rev(survey.poly.mat[,2])))
np.seal.20070618 <- ppp(x = np.seal.20070618.X, y = np.seal.20070618.Y,
                        owin = survey.win) # error

plot(survey.win)
points(x = np.seal.20070618.X, y = np.seal.20070618.Y)

# read in bathymetry data
#bath.rast <- raster(here("bathygrd", "bathyg", "w00100"))
geom.points <- as.data.frame(geom.points)

bath.rast <- raster(here("" "bathymetry.tiff"))
bath.df <- raster::extract(x = bath.rast, y = as.data.frame(geom.points), layer=2, nl=2, df = TRUE)


