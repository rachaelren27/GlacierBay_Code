setwd("/Users/rlr3795/Desktop/GlacierBay_Project")

library(sf)
library(here)
library(tidyverse)
library(terra)
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



