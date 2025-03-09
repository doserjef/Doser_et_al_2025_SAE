# 4a-get-pred-data.R: script to extract the coordinates and covariates across
#                     a grid of the Southern US for prediction of biomass 
#                     across the reigon.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)
library(stars)
# devtools::install_github("mikejohnson51/AOI")
library(AOI)
# devtools::install_github("mikejohnson51/climateR")
library(climateR)
library(raster)
library(rasterVis)
library(elevatr)

# Get area of prediction --------------------------------------------------
my.crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = my.crs)
se.us <- usa %>%
  dplyr::filter(ID %in% c('north carolina', 'texas', 'oklahoma', 'virginia',
                          'kentucky', 'tennessee', 'arkansas', 'south carolina',
                          'louisiana', 'mississippi', 'alabama', 'georgia',
                          'florida'))
# Get sf object with counties of NC
usa.county <- st_as_sf(maps::map("county", fill = TRUE, plot = FALSE))
usa.county <- usa.county %>%
  st_transform(crs = my.crs)
se.us.county <- usa.county %>%
  separate_wider_delim(ID, delim = ',', names = c('state', 'county')) %>%
  dplyr::filter(state %in% c('north carolina', 'texas', 'oklahoma', 'virginia',
                             'kentucky', 'tennessee', 'arkansas', 'south carolina',
                             'louisiana', 'mississippi', 'alabama', 'georgia',
                             'florida')) %>%
  mutate(county = pull(unite(data = ., col = 'county', state, county, sep = '-'), county)) %>%
  st_as_sf()
n.counties <- nrow(se.us.county)

# Generate prediction grid and data set -----------------------------------
# The dx x dy indicates the resolution in terms of km.
# Generate a 1x1km grid. 
grid.pred <- st_as_stars(st_bbox(se.us), dx = 1, dy = 1)
# Convert to data frame
coords.pred <- as.data.frame(grid.pred, center = TRUE)
# Convert coordinates to an sf object
coords.pred.sf <- st_as_sf(coords.pred,
                           coords = c('x', 'y'),
                           crs = my.crs)
# Intersect with region of interest
coords.pred.sf <- st_intersection(coords.pred.sf, st_make_valid(se.us))
coords.0 <- as.data.frame(st_coordinates(coords.pred.sf))
coords.lat.long <- coords.pred.sf %>%
  st_transform(crs = '+proj=longlat +datum=WGS84')

# Download and process TerraClim normals ----------------------------------
# Source: TerraClimate (https://www.climatologylab.org/terraclimate.html)
# 30-year climate normals (1981-2010)
# Citation: Abatzoglou et al., (2018) https://doi.org/10.1038/sdata.2017.191

period <- "19812010"
cvars <- c("tmax","tmin","ppt","pet","aet","def","vpd")

nc.vars <- length(cvars)
# Design matrix of associated variables (+ 2 is for elevation and TCC)
# at the prediction locations
n.0 <- nrow(coords.lat.long)
X.0 <- matrix(NA, n.0, nc.vars + 2)

for (i in 1:nc.vars) {
  print(paste0("Currently on variable ", i, " out of ", nc.vars))
  cdat <- getTerraClimNormals(coords.lat.long, cvars[i], period, 1:12)[[1]]
  plt_dat <- extract(cdat, coords.lat.long)
  if (cvars[i] == "tmax"){
    val <- apply(plt_dat[, 2:ncol(plt_dat)], 1, max)
  } else if (cvars[i] == "tmin"){
    val <- apply(plt_dat[, 2:ncol(plt_dat)], 1, min)
  } else if (cvars[i] == "ppt"){
    val <- apply(plt_dat[, 2:ncol(plt_dat)], 1, sum)
  } else if (cvars[i] %in% c("pet", "aet", "def")){
    val <- apply(plt_dat[, 5:9], 1, sum)
  } else {
    val <- apply(plt_dat[, 5:9], 1, mean)
  }
  X.0[,i] <- val
}

# Download elevation data -------------------------------------------------
## Source: Amazon Web Services (AWS) Terrain Tiles (https://registry.opendata.aws/terrain-tiles/)
## Citation: Terrain Tiles was accessed on INSERT DATE from https://registry.opendata.aws/terrain-tiles.

elev <- get_elev_point(coords.lat.long, src = "aws")
X.0[, nc.vars + 1] <- elev$elevation

# Get TCC data from NLCD --------------------------------------------------
# NOTE: this will not run successfully as you will need to download the TCC data
#       from online here: https://data.fs.usda.gov/geodata/rastergateway/treecanopycover/.
# Tree Canopy Cover in 2019.
tcc.us <- read_stars("~/Dropbox/data/tcc-nlcd/science_tcc_conus_2019_v2021-4.tif")
coords.tcc <- coords.lat.long %>%
  st_transform(crs = st_crs(tcc.us))
tcc.fia.sf <- st_extract(tcc.us, at = coords.tcc)
tcc.val <- tcc.fia.sf$`science_tcc_conus_2019_v2021-4.tif`
X.0[, nc.vars + 2] <- tcc.val

# Remove sites with NA values ---------------------------------------------
# Remove sites with any NA values
bad.indx <- which(apply(X.0, 1, function(a) sum(is.na(a))) > 0)
bad.indx <- c(bad.indx, which(X.0[, nc.vars + 2] > 100))
coords.0 <- coords.0[-bad.indx, ]
coords.pred.sf <- coords.pred.sf[-bad.indx, ]
X.0 <- X.0[-bad.indx, ]

# Determine the county that each point falls into -------------------------
indx.by.county <- st_contains(se.us.county, coords.pred.sf)
# Fill in county vector with the appropriate information
county <- vector(mode = 'character', length = nrow(X.0))
for (i in 1:n.counties) {
  county[indx.by.county[[i]]] <- se.us.county$county[i]
}

# Column names
colnames(X.0) <- c(cvars, 'elev', 'tcc')

# Save results to hard drive ----------------------------------------------
save(X.0, coords.0, county, file = 'data/se-prediction-data.rda')
