# 3d-get-kNN-estimates.R: this script implements a non-parametric approach to
#                         multivariate small area estimation using a k nearest
#                         neighbors approach.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(yaImpute)

load('data/sim_pop_data.rda')
y.0 <- as.data.frame(t(ifelse(y != 0, exp(y), y)))
covs.0 <- data.frame(tmin = X.1[, 'tmin'], 
                     tmax = X.1[, 'tmax'], 
                     ppt = X.1[, 'ppt'], 
                     tcc = X.2[, 'tcc'], 
                     vpd = X.2[, 'vpd'], 
                     elev = X.2[, 'elev'], 
                     easting = coords.0[, 1],
                     northing = coords.0[, 2])
rownames(covs.0) <- paste0('new-', 1:nrow(covs.0))
n.sims <- 100
knn.ests <- list()
for (i in 1:n.sims) {
  load(paste0('data/sim_data_replicates/sim_', i, '_data.rda'))
  y.fit <- as.data.frame(t(ifelse(data.list.2$y != 0, exp(data.list.2$y), data.list.2$y)))
  covs.fit <- data.frame(tmin = data.list.1$covs[, 'tmin'], 
                         tmax = data.list.1$covs[, 'tmax'], 
                         ppt = data.list.1$covs[, 'ppt'], 
                         tcc = data.list.2$covs[, 'tcc'], 
                         vpd = data.list.2$covs[, 'vpd'], 
                         elev = data.list.2$covs[, 'elev'], 
                         easting = data.list.2$coords[, 1],
                         northing = data.list.2$coords[, 2])
  out.k.mal <- yai(x = covs.fit, y = y.fit, method = 'mahalanobis', k = 10)
  pred.out.yai <- predict(out.k.mal, newdata = covs.0)
  k.means.pred.out <- pred.out.yai %>%
    as.data.frame() %>%
    select(!ends_with('.o')) %>%
    mutate(county = county) %>% 
    arrange(county)
  attr(k.means.pred.out, 'scale') <- NULL
  k.means.county.ests <- k.means.pred.out %>%
    group_by(county) %>%
    summarize(across(where(is.numeric), mean)) %>%
    select(-county)
  ests <- t(as.matrix(k.means.county.ests))
  knn.ests[[i]] <- ests
}

# Save to output file -----------------------------------------------------
save(knn.ests, file = 'results/sim_knn_ests.rda')

