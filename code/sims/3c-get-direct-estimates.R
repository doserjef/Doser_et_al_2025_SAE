# 3c-get-direct-estimates.R: script to extract the direct estimates for each of 
#                            the 100 replicate simulated data sets. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)

# Data prep ---------------------------------------------------------------
# Number of simulations
n.sims <- 100
direct.ests <- list()
direct.sds <- list()
for (i in 1:n.sims) {
  load(paste0('data/sim_data_replicates/sim_', i, '_data.rda'))
  n.sp <- nrow(data.list.1$y)
  biomass <- ifelse(data.list.2$y == 0, 0, exp(data.list.2$y))
  dat.long <- data.frame(biomass = c(t(biomass)), 
                         species = rep(1:n.sp, each = nrow(data.list.1$coords)),
                         county = rep(data.list.1$covs$county, times = n.sp)) %>%
    arrange(county)
  direct.ests[[i]] <- dat.long %>% 
    group_by(species, county) %>%
    summarize(bio = mean(biomass)) %>%
    ungroup()
  direct.sds[[i]] <- dat.long %>%
    group_by(species, county) %>%
    summarize(sd = sqrt(sum((biomass - mean(biomass))^2) / (n() * (n() - 1)))) %>%
    ungroup()
}


# Save to hard drive ------------------------------------------------------
save(direct.ests, direct.sds, file = 'results/sim_direct_ests.rda')

