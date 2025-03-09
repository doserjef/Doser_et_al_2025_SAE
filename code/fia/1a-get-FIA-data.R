# 1a-get-FIA-data.R: this script extracts FIA data from rFIA for use in the
#                    analysis.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(rFIA)
library(sf)
# For preventing timeout when downloading data
options(timeout=3600)

# NOTE: if running this script, you may not produce exactly the same 
#       set of data used in the case study as a result of the timing of when 
#       you download the FIA data from the US Forest Service Data Market. 
#       For our analysis, we downloaded the data in November 2024. The 
#       object created from this script that we used in the analysis is included
#       on GitHub ('data/se_fia_bio_data.rda')

# Set directories ---------------------------------------------------------
# NOTE: change this for running on your own computer
# Directory for where to save and load FIA data
fia.data.dir <- '~/Dropbox/data/fia/'

# Read in the FIA data ----------------------------------------------------
# NOTE: uncomment this code to download the FIA data for all states used 
#       in the case study. 
# Only needs to be done once to download FIA data, then can just load 
# in the data from your local machine.
# WV <- getFIA(states = 'WV', dir = fia.data.dir)
# VA <- getFIA(states = 'VA', dir = fia.data.dir)
# NC <- getFIA(states = 'NC', dir = fia.data.dir)
# TN <- getFIA(states = 'TN', dir = fia.data.dir)
# AR <- getFIA(states = 'AR', dir = fia.data.dir)
# LA <- getFIA(states = 'LA', dir = fia.data.dir)
# MS <- getFIA(states = 'MS', dir = fia.data.dir)
# AL <- getFIA(states = 'AL', dir = fia.data.dir)
# GA <- getFIA(states = 'GA', dir = fia.data.dir)
# SC <- getFIA(states = 'SC', dir = fia.data.dir)
# FL <- getFIA(states = 'FL', dir = fia.data.dir)
# KY <- getFIA(states = 'KY', dir = fia.data.dir)
# TX <- getFIA(states = 'TX', dir = fia.data.dir)
# OK <- getFIA(states = 'OK', dir = fia.data.dir)

# Load the data -----------------------------------------------------------
fia.dat <- readFIA(dir = fia.data.dir, inMemory = FALSE,
                   states = c('NC', 'VA', 'TN', 'AR', 'LA', 'MS', 'AL',
                              'GA', 'SC', 'FL', 'KY', 'OK', 'TX'))

# Extract biomass by plot -------------------------------------------------
# Get total aboveground biomass for all live and dead stems measured greater than
# 1 in. in DBH. Live/dead includes all stems greater than 1 in. DBH which are
# live or dead (leaning less than 45 degrees).
biomass.se <- biomass(fia.dat, byPlot = TRUE, bySpecies = TRUE, returnSpatial = TRUE,
                      treeType = 'all')
# Filter data to only include the most recent record of the given plot.
biomass.se <- biomass.se %>%
  group_by(pltID) %>%
  slice_max(YEAR)
# Filter to only use data over the last 10 years (2015-2024)
biomass.se <- biomass.se %>%
  filter(YEAR >= 2015) %>%
  ungroup()
# Extract coordinates and unique plot ids ---------------------------------
my.crs <- st_crs(biomass.se)
coords.df <- biomass.se %>%
  st_coordinates()
coords.link.df <- data.frame(lon = coords.df[, 1],
                             lat = coords.df[, 2],
                             pltID = biomass.se$pltID)
coords.unique.df <- unique(coords.link.df)
# Remove lat/lons that have two different plot IDs
tmp <- coords.unique.df %>%
  mutate(site = paste(lon, lat)) %>%
  group_by(site) %>%
  summarize(vals = n()) %>%
  arrange(desc(vals))

# NOTE: this is hardcoded
bad.indx <- c(which(coords.unique.df$lon == -82.357975 &
                    coords.unique.df$lat == 29.675897),
              which(coords.unique.df$lon == -83.174799 &
                    coords.unique.df$lat == 29.989016),
              which(coords.unique.df$lon == -83.31699 &
                    coords.unique.df$lat == 36.543153), 
              which(coords.unique.df$lon == -83.174799 & 
                    coords.unique.df$lat == 29.989016))
coords.unique.df <- coords.unique.df[-bad.indx, ]
# Drop spatial geometries from objects to make the following go much faster
biomass.se <- biomass.se %>% st_drop_geometry()

# Get final data set in desired format ------------------------------------
# Add 0s for species at pltIDs that have been sampled in a given YEAR
biomass.se <- biomass.se %>%
  select(-CARB_ACRE, -PROP_FOREST, -SPCD, -SCIENTIFIC_NAME) %>%
  complete(COMMON_NAME, nesting(pltID, YEAR, PLT_CN),
           fill = list(BIO_ACRE = 0))

# Filter to only include the 20 most abundant species
bio.sum.by.sp <- biomass.se %>%
  group_by(COMMON_NAME) %>%
  summarize(total_bio = sum(BIO_ACRE)) %>%
  arrange(desc(total_bio)) %>%
  slice(1:20)
biomass.se <- biomass.se %>%
  filter(COMMON_NAME %in% bio.sum.by.sp$COMMON_NAME)

# Pivot to wider format
final.wide.dat <- biomass.se %>%
  pivot_wider(names_from = COMMON_NAME, values_from = c(BIO_ACRE))

# Join to coordinates to get spatial information --------------------------
final.wide.dat <- inner_join(final.wide.dat, coords.unique.df, by = 'pltID')

final.wide.dat <- st_as_sf(final.wide.dat,
                           coords = c('lon', 'lat'),
                           crs = my.crs)

# Save to hard drive ------------------------------------------------------
save(final.wide.dat, file = 'data/se_fia_bio_data.rda')
