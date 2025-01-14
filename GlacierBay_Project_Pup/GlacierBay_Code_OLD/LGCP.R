library(here)
library(sf)
library(lgcp)
library(tidyverse)
library(spatstat)

path <- here("Data", "HARBORSEAL_2007", "seal_locations_final",
             "nonpup_locs")
np.seal.locs.20070618 <- st_read(dsn = path, layer = "JHI_20070618_nonpup_locs")

ggplot() +
  geom_sf(data = survey.poly.20070618) + 
  geom_sf(data = np.seal.locs.20070618)

geom.points <- st_coordinates(np.seal.locs.20070618$geometry) # these are not coordinates
np.seal.20070618.X <- geom.points[,1]
np.seal.20070618.Y <- geom.points[,2]

survey.poly.mat <- survey.poly.20070618$geometry[[1]][[1]]
survey.win <- owin(poly = data.frame(x=rev(survey.poly.mat[,1]),
                                     y=rev(survey.poly.mat[,2])), check = FALSE)
pp <- ppp(x = np.seal.20070618.X, y = np.seal.20070618.Y,
                        owin = survey.win, check = FALSE)

# mcmc.pars <- mcmcpars(mala.length = 1000, burnin = 100, retain = 1,
#                       adaptivescheme = constanth(0.001))
# lg <- lgcpPredictSpatial(pp, cellwidth = 1, mcmc.control = mcmc.pars, spatial.intensity = )

# lgcp <- lgcp.estpcf(pp, startpar=c(var=1,scale=1),
#             covmodel=list(model="exponential"))


