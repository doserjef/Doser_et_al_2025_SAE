# 1d-case-study-data-plots.R: some exploratory data analysis plots as well as code
#                             to generate Figure 1 in the manuscript.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)
library(patchwork)
library(viridis)
library(RColorBrewer)

# Load the data -----------------------------------------------------------
load('data/se_bio_stage_1_data.rda')
load('data/se_bio_stage_2_data.rda')

# Map of the inventory plots ---------------------------------------------- 
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
my.crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"
usa <- usa %>%
  st_transform(crs = my.crs)
my.state <- usa %>%
  dplyr::filter(ID %in% c('north carolina', 'virginia', 'kentucky', 'tennessee', 
                          'arkansas', 'south carolina',
                          'louisiana', 'mississippi', 'alabama', 'georgia',
                          'florida', 'oklahoma', 'texas'))
usa.county <- st_as_sf(maps::map("county", fill = TRUE, plot = FALSE))
usa.county <- usa.county %>%
  st_transform(crs = my.crs)
se.county <- usa.county %>%
  separate_wider_delim(ID, delim = ',', names = c('state', 'county')) %>%
  dplyr::filter(state %in% c('north carolina', 'virginia',
                             'kentucky', 'tennessee', 'arkansas', 'south carolina',
                             'louisiana', 'mississippi', 'alabama', 'georgia',
                             'florida', 'oklahoma', 'texas')) %>%
  st_as_sf()
coords.sf <- st_as_sf(data.frame(data.list.2$coords),
                      coords = c('X', 'Y'),
                      crs = my.crs)

# Map showing the fuzzed plot locations
ggplot(coords.sf) +
  geom_sf(data = se.county, fill = 'white', color = 'grey', lwd = 0.4) +
  geom_sf(size = 0.02, col = 'black') +
  theme_bw(base_size = 18) +
  labs(x = 'Longitude', y = 'Latitude') +
  theme(text = element_text(family="LM Roman 10"))
# ggsave(file = 'figures/plot-locations.png', width = 7, height = 7, units = 'in',
#        bg = 'white')

# Figure 1 (map showing number of plots by county) ------------------------
# Determine number of points in each county
indx.by.county <- st_contains(se.county, coords.sf)
se.county$n.plots <- sapply(indx.by.county, length)

se.county$n.plots <- ifelse(se.county$n.plots > 100, 100, se.county$n.plots)

ggplot(se.county) +
  geom_sf(aes(fill = n.plots)) +
  theme_bw(base_size = 18) +
  # scale_fill_stepsn(colors = rev(brewer.pal(9, 'Blues')),
  #                   breaks = c(0, 10, 30, 50, 100, 135)) +
  scale_fill_gradient(low = '#F7FBFF', high = '#08306B', 
                      labels = c('0', '25', '50', '75', '>100')) +
  # scale_fill_gradient(low = '#67000D', high = '#FFF5F0') +
  # scale_fill_brewer(palette = 'Blues') +
  labs(x = 'Longitude', y = 'Latitude', fill = 'Sample\nSize') +
  theme(text = element_text(family="LM Roman 10"),
        legend.position = 'inside',
        legend.position.inside = c(0.93, 0.4),
        legend.background = element_blank())
ggsave(file = 'figures/Figure-1.png', width = 10, height = 7, units = 'in',
       bg = 'white')
