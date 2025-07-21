# 6a-subset-data.R: script to subset the FIA data across the Southeast for 
#                   prediction in counties where no sample data is used. This
#                   script generates four distinct hold out and prediction sets
#                   for use in k-fold cross-validation.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)

# Load the full data set --------------------------------------------------
load('data/se_bio_stage_1_data.rda')
load('data/se_bio_stage_2_data.rda')

# Determine the county that each point falls into -------------------------
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
my.crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"
usa <- usa %>%
  st_transform(crs = my.crs)
my.state <- usa %>%
  dplyr::filter(ID %in% c('north carolina', 'virginia', 'kentucky', 'tennessee',
                          'arkansas', 'south carolina',
                          'louisiana', 'mississippi', 'alabama', 'georgia',
                          'florida', 'oklahoma', 'texas'))
coords.sf <- st_as_sf(data.frame(data.list.2$coords),
                      coords = c('X', 'Y'),
                      crs = my.crs)
usa.county <- st_as_sf(maps::map("county", fill = TRUE, plot = FALSE))
usa.county <- usa.county %>%
  st_transform(crs = my.crs)
se.county <- usa.county %>%
  separate_wider_delim(ID, delim = ',', names = c('state', 'county')) %>%
  dplyr::filter(state %in% c('north carolina', 'virginia',
                             'kentucky', 'tennessee', 'arkansas', 'south carolina',
                             'louisiana', 'mississippi', 'alabama', 'georgia',
                             'florida', 'oklahoma', 'texas')) %>%
  mutate(county = pull(unite(data = ., col = 'county', state, county, sep = '-'), county)) %>%
  st_as_sf()
indx.by.county <- st_contains(se.county, coords.sf)
# Fill in county vector with the appropriate information
n.counties <- nrow(se.county)
n.plots <- nrow(data.list.1$coords)
county <- vector(mode = 'character', length = n.plots)
for (i in 1:n.counties) {
  county[indx.by.county[[i]]] <- se.county$county[i]
}

# Number of folds and associated ------------------------------------------
# NOTE: hardcoded
n.folds <- 4
set.seed(19191)
full.indx <- sample(1:n.counties, n.counties, replace = FALSE)
indices <- split(full.indx, sort(full.indx %% n.folds))

# Get the testing and hold out data set -----------------------------------
# Loop through and create a new data set each time.
for (i in 1:n.folds) {
  n.obs.counties <- n_distinct(county)
  # Determine the counties that you're going to keep
  hold.counties.indx <- indices[[i]] 
  fit.counties.indx <- (1:n.counties)[-hold.counties.indx]
  fit.counties <- se.county$county[fit.counties.indx]
  hold.counties <- se.county$county[hold.counties.indx]
  fit.indx <- which(county %in% fit.counties)
  hold.indx <- which(county %in% hold.counties)
  # The values that aren't in fit.indx and hold.indx fall outside of a county, which is 
  # likely do to the fuzzing. Assigning those to the fit data. 
  no.county.indx <- which(county == '')
  fit.indx <- c(fit.indx, no.county.indx)
  
  # Generate the necessary data sets ----------------------------------------
  # Fit data ----------------------------
  data.fit.1 <- data.list.1
  data.fit.1$y <- data.list.1$y[, fit.indx]
  data.fit.1$coords <- data.list.1$coords[fit.indx, ]
  data.fit.1$covs <- data.list.1$covs[fit.indx, ]
  data.fit.2 <- data.list.2
  data.fit.2$y <- data.list.2$y[, fit.indx]
  data.fit.2$coords <- data.list.2$coords[fit.indx, ]
  data.fit.2$covs <- data.list.2$covs[fit.indx, ]
  data.fit.2$z <- data.list.2$z[, fit.indx]
  save(data.fit.1, file = paste0('data/cross_val_data/hold_', i, '_stage_1_data_validation_fit.rda'))
  save(data.fit.2, file = paste0('data/cross_val_data/hold_', i, '_stage_2_data_validation_fit.rda'))
  
  # Hold out data -----------------------
  data.hold.1 <- data.list.1
  data.hold.1$y <- data.list.1$y[, hold.indx]
  data.hold.1$coords <- data.list.1$coords[hold.indx, ]
  data.hold.1$covs <- data.list.1$covs[hold.indx, ]
  data.hold.2 <- data.list.2
  data.hold.2$y <- data.list.2$y[, hold.indx]
  data.hold.2$coords <- data.list.2$coords[hold.indx, ]
  data.hold.2$covs <- data.list.2$covs[hold.indx, ]
  data.hold.2$z <- data.list.2$z[, hold.indx]
  save(data.hold.1, file = paste0('data/cross_val_data/hold_', i, '_stage_1_data_validation_hold.rda'))
  save(data.hold.2, file = paste0('data/cross_val_data/hold_', i, '_stage_2_data_validation_hold.rda'))
}

