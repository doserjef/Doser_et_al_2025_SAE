# 5a-get-direct-estimates.R: script to extract the design-based estimates
#                            for county-level and species-specific
#                            biomass estimates across the South.
#                            Estimates are calculated
#                            using a standard Horvitz-Thompson estimator
#                            for each species (i.e., the basic sample mean).
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)

# Data prep ---------------------------------------------------------------
load("data/se_bio_stage_2_data.rda")
y.raw <- apply(data.list.2$y, 1, function(a) ifelse(a == 0, 0, exp(a)))
dat <- as.data.frame(y.raw)
dat$x <- data.list.2$coords[, 1]
dat$y <- data.list.2$coords[, 2]

# Get county grid for calculating SAEs ------------------------------------
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
se.county <- usa.county %>%
  separate_wider_delim(ID, delim = ',', names = c('state', 'county')) %>%
  dplyr::filter(state %in% c('north carolina', 'texas', 'oklahoma', 'virginia',
                          'kentucky', 'tennessee', 'arkansas', 'south carolina',
                          'louisiana', 'mississippi', 'alabama', 'georgia',
                          'florida')) %>%
  mutate(county = pull(unite(data = ., col = 'county', state, county, sep = '-'), county)) %>%
  st_as_sf()
n.counties <- nrow(se.county)

# Coordinates in sf
coords.sf <- st_as_sf(as.data.frame(data.list.2$coords),
                      coords = c('X', 'Y'),
                      crs = my.crs)


indx.by.county <- st_contains(se.county, coords.sf)
# Fill in county vector with the appropriate information
n.counties <- nrow(se.county)
n.plots <- nrow(dat)
county <- vector(mode = 'character', length = n.plots)
for (i in 1:n.counties) {
  county[indx.by.county[[i]]] <- se.county$county[i]
}
unique.counties <- unique(county)

dat <- dat %>%
  mutate(county = factor(county, levels = unique.counties)) %>%
  arrange(county)

# Convert to long format --------------------------------------------------
# NOTE: hardcoded with species names
dat.long <- dat %>%
  pivot_longer(cols = `longleaf pine`:`yellow-poplar`, names_to = 'species',
               values_to = 'biomass')

# Get output for comparison to model-based values -------------------------
# Mean --------------------------------
sample.mean.ests <- dat %>%
  group_by(county, .drop = FALSE) %>%
  summarize(across(`longleaf pine`:`yellow-poplar`, mean),
            n.j = n())
sample.mean.ests.long <- sample.mean.ests %>%
  pivot_longer(cols = `longleaf pine`:`yellow-poplar`,
               names_to = 'species', values_to = 'mu')
direct.ests <- sample.mean.ests.long %>%
  filter(!is.na(county))

# Standard error ----------------------
se.direct.ests <- dat %>%
  group_by(county) %>%
  summarize(across(`longleaf pine`:`yellow-poplar`, 
                   function(a) sqrt(sum((a - mean(a))^2) / (n() * (n() - 1)))))

se.ests.long <- se.direct.ests %>%
  pivot_longer(cols = `longleaf pine`:`yellow-poplar`,
               names_to = 'species', values_to = 'mu')
direct.se.ests <- se.ests.long %>%
  filter(!is.na(county))

# Save to hard drive ------------------------------------------------------
save(direct.ests, direct.se.ests, file = 'results/se_direct_estimates.rda')

