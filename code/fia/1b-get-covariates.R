# 1b-get-covariates.R: script to extract the covariates used in the associated
#                      case study.
# Author: Jeffrey W. Doser, Malcolm S. Itter
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

# Load the FIA data -------------------------------------------------------
load("data/se_fia_bio_data.rda")
coords.df <- final.wide.dat %>% st_coordinates() %>% as.data.frame()
coords.sf <- st_as_sf(coords.df,
                      coords = c('X', 'Y'),
                      crs = 4326)
# Number of plots in the analysis
n <- nrow(coords.sf)
site.indx <- final.wide.dat$pltID

# Download and process TerraClim normals ----------------------------------
# Source: TerraClimate (https://www.climatologylab.org/terraclimate.html)
# 30-year climate normals (1981-2010)
# Citation: Abatzoglou et al., (2018) https://doi.org/10.1038/sdata.2017.191
# Data accessed on: October 5, 2024

period <- "19812010"
cvars <- c("tmax","tmin","ppt","pet","aet","def","vpd")

nc.vars <- length(cvars)
# Design matrix of associated variables (+ 2 is for elevation, TCC)
X <- matrix(NA, n, nc.vars + 2)
for (i in 1:nc.vars){
  print(paste0("Currently on variable ", i, " out of ", nc.vars))
  cdat <- getTerraClimNormals(coords.sf, cvars[i], period, 1:12)[[1]]
  plt_dat <- extract(cdat, coords.sf)
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
  X[,i] <- val
}

# Download elevation data -------------------------------------------------
## Source: Amazon Web Services (AWS) Terrain Tiles (https://registry.opendata.aws/terrain-tiles/)
## Citation: Terrain Tiles was accessed on INSERT DATE from https://registry.opendata.aws/terrain-tiles.
elev <- get_elev_point(coords.sf, src = "aws")
all(elev$pltID == site.indx) # quick check that data is in correct order
X[, nc.vars + 1] <- elev$elevation

# Get TCC data from NLCD --------------------------------------------------
# Tree Canopy Cover in 2019.
# NOTE: this will not run successfully as you will need to download the TCC data 
#       from online here: https://data.fs.usda.gov/geodata/rastergateway/treecanopycover/. 
tcc.us <- read_stars("~/Dropbox/data/tcc-nlcd/science_tcc_conus_2019_v2021-4.tif")
coords.tcc <- coords.sf %>%
  st_transform(crs = st_crs(tcc.us))
tcc.fia.sf <- st_extract(tcc.us, at = coords.tcc)
tcc.val <- tcc.fia.sf$`science_tcc_conus_2019_v2021-4.tif`
X[, nc.vars + 2] <- tcc.val

X <- as.data.frame(X)
colnames(X) <- c(cvars, 'elev', 'tcc')

# Save to hard drive ------------------------------------------------------
save(X, site.indx, file = 'data/covariate-data.rda')
