# 1c-get-replicate-data-sets: extract the individual data sets from the overall
#                             population for running multiple simulations with 
#                             different replicates.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)

# Set seed for reproducibility --------------------------------------------
set.seed(773829)

# Read in the full data set -----------------------------------------------
load('data/sim_pop_data.rda')

# Subset the data and save as output in spOcc/spAbund format --------------
# Number of replicate data sets
n.reps <- 100
# Number of inventory plots for each replicate data set.  
n <- 1000
# Total number of plots across the population.
n.full <- nrow(coords.0)
for (i in 1:n.reps) {
  indx <- sample(1:n.full, n, replace = FALSE)
  data.list.1 <- list(y = z[, indx], 
                      covs = data.frame(X.1[indx, ], 
                                        county = county[indx]),
                      coords = coords.0[indx, ])
  data.list.2 <- list(y = y[, indx],
                      covs = X.2[indx, ],
                      coords = coords.0[indx, ],
                      z = z[, indx])
  # Save each data set individually. 
  save(data.list.1, data.list.2, file = paste0('data/sim_data_replicates/sim_', 
                                               i, '_data.rda')) 
}
