# 5b-get-kNN-estimates.R: this script implements a non-parametric approach to
#                         multivariate small area estimation using a k nearest
#                         neighbors approach.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(yaImpute)

# Read in the stage 2 (biomass) data --------------------------------------
load('data/se_bio_stage_2_data.rda')

# Run the k-means imputation ----------------------------------------------
# Back-transform the species-specific biomass variable to now hold the
# true biomass instead of the log values
y <- as.data.frame(t(ifelse(data.list.2$y != 0, exp(data.list.2$y), data.list.2$y)))
# colnames(y) <- str_replace_all(colnames(y), ' ', '.')
# colnames(y) <- str_replace_all(colnames(y), '-', '.')
covs <- data.list.2$covs[, c('tmin', 'tmax', 'tcc', 'vpd', 'ppt', 'elev')]
covs$easting <- data.list.2$coords[, 1]
covs$northing <- data.list.2$coords[, 2]

y.fit <- y
covs.fit <- covs

out.k.mal <- yai(x = covs.fit, y = y.fit, method = 'mahalanobis', k = 10)

# Predict at new locations ------------------------------------------------
load('data/se-prediction-data.rda')

covs.0 <- data.frame(tmin = X.0[, 'tmin'], 
                     tmax = X.0[, 'tmax'], 
                     tcc = X.0[, 'tcc'],
                     vpd = X.0[, 'vpd'],
                     ppt = X.0[, 'ppt'],
                     elev = X.0[, 'elev'],
                     easting = coords.0[, 1],
                     northing = coords.0[, 2])
# predict.yai() works by imputing
rownames(covs.0) <- paste0('new-', 1:nrow(covs.0))

pred.out.yai <- predict(out.k.mal, newdata = covs.0)

# Aggregate predictions at the county level -------------------------------
k.means.pred.out <- pred.out.yai %>%
  as.data.frame() %>%
  select(!ends_with('.o')) %>%
  mutate(county = county)
attr(k.means.pred.out, 'scale') <- NULL

k.means.county.ests <- k.means.pred.out %>%
  group_by(county) %>%
  summarize(across(where(is.numeric), mean))

# Save to output file -----------------------------------------------------
save(k.means.county.ests, file = 'results/k-means-sae-estimates.rda')

