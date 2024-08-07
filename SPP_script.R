setwd("/Users/rlr3795/Desktop/GlacierBay_Project")

library(sf)
library(here)
library(tidyverse)
library(terra)
library(raster)
library(hilbertSimilarity)
library(geosphere)

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

survey.poly.mat <- survey.poly[[1]][[1]]
survey.win <- owin(poly = data.frame(x=rev(survey.poly.mat[,1]),
                                     y=rev(survey.poly.mat[,2])))

# plot
ggplot() + 
  geom_sf(data = survey.poly) + 
  geom_sf(data = footprint) + 
  geom_sf(data = seal.locs, size = 0.5)

# --- Read in Covariates -------------------------------------------------------
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
