setwd("/Users/rlr3795/Desktop/GlacierBay_Project")

library(sf)
library(here)
library(tidyverse)
library(stringi)

# --- Read in NPS data ---------------------------------------------------------

years <- c(2007:2019, 2022)
n.years <- length(years)

# get nonpup survey dates
non.pup.dates <- list()
for(i in 1:(n.years-3)){
  year <- years[i]
  dirname <- here("NPS_data", paste0("HARBORSEAL_", year), "seal_locations_final",
                  "nonpup_locs")
  files <- list.files(path = dirname)
  non.pup.dates[[i]] <- unique(lapply(files, substr, start = 9, stop = 12))
}
for(i in (n.years-2):n.years){
  year <- years[i]
  dirname <- here("NPS_ata", paste0("HARBORSEAL_", year), "seal_locations_final",
                  "non_pup_locs")
  files <- list.files(path = dirname)
  non.pup.dates[[i]] <- unique(lapply(files, substr, start = 9, stop = 12))
}

# get pup survey dates
pup.dates <- list()
for(i in 1:(n.years)){
  year <- years[i]
  dirname <- here("NPS_data", paste0("HARBORSEAL_", year), "seal_locations_final",
                  "pup_locs")
  files <- list.files(path = dirname)
  pup.dates[[i]] <- unique(lapply(files, substr, start = 9, stop = 12))
}
 
# example all seal location (6/18/2007)
path <- here("NPS_data", "HARBORSEAL_2007", "seal_locations_final",
                            "all_seal_locs")
all.seal.locs.20070618 <- st_read(dsn = path, layer = "JHI_20070618_allseal_locs")

ggplot(all.seal.locs.20070618, aes(shape = seal_type)) + 
  geom_sf()

# read in non-pup data
# non.pup.locs <- list()
# for(year in years){
#   for(date in non.pup.dates){
#     path <- here("Data", paste0("HARBORSEAL_", year), "seal_locations_final",
#                  "nonpup_locs")
#     temp.df <- st_read(dsn = path, layer = paste0("JHI_", year, date, "_nonpup_locs"))
#     stri_sub(date, 3, 2)  <- "-"
#     temp.df <- temp.df %>% mutate(date = as.Date(paste0(date, '-', year), format = '%m-%d-%Y'))
#     non.pup.locs <- rbind(non.pup.locs, temp.df)
#   }
# } # error with rbind

# example footprint
path <- here("NPS_data", "HARBORSEAL_2007", "footprints")
footprint.20070618 <- st_read(dsn = path, layer = "JHI_20070618_footprint")

ggplot(footprint.20070618) +
  geom_sf()

# example image location
path <- here("NPS_data", "HARBORSEAL_2007", "image_locations")
image.loc.20070618 <- st_read(dsn = path, layer = "XYJHI_20070618_exif_NAD83")

ggplot() +
  geom_sf(data = footprint.20070618) + 
  geom_sf(data = image.loc.20070618, size = 0.3)

path <- here("NPS_data", "HARBORSEAL_2007", "image_locations")
image.loc.20070618.2 <- st_read(dsn = path, layer = "XYJHI_20070618_exif")

ggplot() +
  geom_sf(data = footprint.20070618) + 
  geom_sf(data = image.loc.20070618.2, size = 0.3)
  
# example joined image location
path <- here("NPS_data", "HARBORSEAL_2007", "joined_image_locations")
joined.image.loc.20070618 <- st_read(dsn = path, layer = "XYJHI_20070618_exif_NAD83_join")

ggplot() +
  geom_sf(data = joined.image.loc.20070618, size = 0.3) # what's the difference?

# example survey polygons
path <- here("NPS_data", "HARBORSEAL_2007", "survey_polygons")
survey.poly.20070618 <- st_read(dsn = path, layer = "JHI_20070618_surveyboundary")

ggplot() +
  geom_sf(data = survey.poly.20070618) + 
  geom_sf(data = footprint.20070618)
  
path <- here("NPS_data", "HARBORSEAL_2007", "footprints")
footprint.20070619 <- st_read(dsn = path, layer = "JHI_20070619_footprint")
path <- here("NPS_data", "HARBORSEAL_2007", "survey_polygons")
survey.poly.20070619 <- st_read(dsn = path, layer = "JHI_20070619_surveyboundary")

ggplot() +
  geom_sf(data = survey.poly.20070619) + 
  geom_sf(data = footprint.20070619)

# basemaps
path <- here("NPS_data", "HARBORSEAL_2007", "basemaps")
gbnp.cutter.NAD83 <- st_read(dsn = path, layer = "gbnp_cutter_NAD83")

ggplot() + 
  geom_sf(data = gbnp.cutter.NAD83) + 
  geom_sf(data = survey.poly.20070618) + 
  geom_sf(data = footprint.20070618)

path <- here("NPS_data", "HARBORSEAL_2007", "basemaps")
gbnp.cutter <- st_read(dsn = path, layer = "gbnp_cutter")

ggplot() + 
  geom_sf(data = gbnp.cutter) + 
  geom_sf(data = survey.poly.20070618) + 
  geom_sf(data = footprint.20070618)

# GBNP land outline
path <- here("NPS_data", "HARBORSEAL_2007", "basemaps")
gbnp.NAD83 <- st_read(dsn = path, layer = "glacier_bay_np_NAD83")

ggplot() + 
  geom_sf(data = gbnp.NAD83) + 
  geom_sf(data = survey.poly.20070618, col = "red")

# read in shapefiles
path <- here("NPS_data", "HARBORSEAL_2007", "shapefiles")
shapefile.20070620 <- st_read(dsn = path, layer = "JHI_20070620_seals")

ggplot() + 
  geom_sf(data = shapefile.20070620)
# are these just for testing?

# --- Convert to same CRS ------------------------------------------------------
survey.poly <- st_transform(survey.poly.20070618$geometry, 
                            CRS("+proj=longlat +datum=WGS84"))
survey.poly.mat <- survey.poly[[1]][[1]]
survey.win <- owin(poly = data.frame(x=rev(survey.poly.mat[,1]),
                                     y=rev(survey.poly.mat[,2])))

seal.locs <- st_transform(all.seal.locs.20070618$geometry,
                          CRS("+proj=longlat +datum=WGS84"))
seal.coord <- st_coordinates(seal.locs)

seal.X <- seal.coord[,1]
seal.Y <- seal.coord[,2]
# all.seal.20070618 <- ppp(x = seal.X, y = seal.Y,
#                          owin = survey.win) # warning

grid.tf <- inside.owin(x = seal.X, y = seal.Y, survey.win) # all true!

footprint <- st_transform(footprint.20070618$geometry, 
                          CRS("+proj=longlat +datum=WGS84"))

footprint.coord <- matrix(data = footprint[[1]][[1]], nrow = 5, ncol = 2)
for(i in 2:length(footprint)){
  footprint.coord <- rbind(footprint.coord, footprint[[i]][[1]])
}

footprint.win <- owin(poly = data.frame(x=rev(footprint.coord[,1]),
                                        y=rev(footprint.coord[,2])))
grid.tf <- inside.owin(x = seal.X, y = seal.Y, footprint.win) # all but two
grid.f.ind <- which(grid.tf == FALSE)

ggplot() +
  geom_sf(data = survey.poly) +
  geom_point(data = as.data.frame(seal.coord[grid.f.ind,]), aes(x = X, y = Y), 
             size = 0.5)

# --- Read in covariate data ---------------------------------------------------
bath.rast <- raster(here("covariates", "bathymetry.tiff"))

# convert footprint crs
footprint.20070618 <- st_transform(footprint.20070618, 
                                   CRS("+proj=longlat +datum=WGS84"))

plot(bath.rast)
plot(st_geometry(survey.poly), add = TRUE, reset = FALSE)
plot(st_geometry(footprint.20070618), add = TRUE, reset = FALSE)

