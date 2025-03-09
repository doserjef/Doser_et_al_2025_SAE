# 3b-get-pop-means.R: script to extract the true population means
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)

# Directories -------------------------------------------------------------
# NOTE: change as needed
out.dir <- 'results/'
data.dir <- 'data/'

# Data prep ---------------------------------------------------------------
load(paste0(data.dir, "sim_pop_data.rda"))
n.sp <- nrow(biomass)
dat.long <- data.frame(biomass = c(t(biomass)), 
                       species = rep(1:n.sp, each = nrow(coords.0)),
                       x = coords.0[, 1], 
                       y = coords.0[, 2], 
                       county = rep(county, times = n.sp)) %>%
  arrange(county)

check.1 <- dat.long %>% group_by(species) %>% summarize(val = mean(biomass))
check.2 <- apply(biomass, 1, mean)
all.equal(check.1$val, check.2)

# Get true population estimates -------------------------------------------
pop.vals <- dat.long %>% 
  group_by(species, county) %>%
  summarize(bio = mean(biomass)) %>%
  ungroup()

# Save to hard drive ------------------------------------------------------
save(pop.vals, file = 'data/sim_pop_county_true.rda')
