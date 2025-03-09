# 4-summary.R: summarize results from the simulation study in a couple of figures
#              that highlight the bias and uncertainty of the proposed estimates.  
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)

# Load in all data objects ------------------------------------------------
# True population data
load('data/sim_pop_county_true.rda')
# Read in direct estimates (loads direct.ests and direct.sds)
load('results/sim_direct_ests.rda')
# Read in knn estimates
load('results/sim_knn_ests.rda')
# Read in model-based estimates
n.sims <- length(direct.ests)
model.ests <- list()
model.sds <- list()
model.cvs <- list()
for (i in 1:n.sims) {
  load(paste0('results/sim_results/sim-biomass-', i, '-replicate-sae-county-means.rda'))
  model.ests[[i]] <- apply(biomass.mean.samples, c(2, 3), median)
  model.sds[[i]] <- apply(biomass.mean.samples, c(2, 3), sd)
  model.cvs[[i]] <- model.sds[[i]] / model.ests[[i]]
}

# Calculate bias and summarize --------------------------------------------
model.bias <- sapply(model.ests, function(a) c(t(a)) - pop.vals$bio)
knn.bias <- sapply(knn.ests, function(a) c(t(a)) - pop.vals$bio)
# Calculate direct estimator bias, which requires more work as some 
# counties don't always have at least one plot.
direct.bias <- array(NA, dim = dim(model.bias))
for (i in 1:n.sims) {
  tmp <- direct.ests[[i]] %>% rename(direct = bio)
  tmp.2 <- left_join(pop.vals, tmp, by = c('species', 'county'))
  direct.bias[, i] <- tmp.2$direct - tmp.2$bio
}
n.sp <- n_distinct(pop.vals$species)
sp.means <- pop.vals %>%
  group_by(species) %>%
  summarize(val = mean(bio)) %>%
  arrange(val) %>%
  pull(species)
bias.plot.df <- data.frame(species = factor(rep(rep(1:n.sp, each = 100), n.sims), 
                                            levels = sp.means, ordered = TRUE),
                           county = rep(pop.vals$county, n.sims), 
                           model.bias = c(model.bias),  
                           direct.bias = c(direct.bias), 
                           knn.bias = c(knn.bias))
bias.plot.df <- data.frame(species = factor(rep(1:n.sp, each = 100), levels = sp.means, 
                                            ordered = TRUE),
                           model.bias = apply(model.bias, 1, mean),  
                           direct.bias = apply(direct.bias, 1, mean), 
                           knn.bias = apply(knn.bias, 1, mean))
# Remove NAs situation when a species was never observed in a plot
bias.plot.df <- bias.plot.df %>%
  filter(!is.na(direct.bias))

bias.plot.df %>%
  group_by(species) %>%
  summarize(model = mean(model.bias), 
            direct = mean(direct.bias), 
            knn = mean(knn.bias))

bias.plot.df.long <- bias.plot.df %>%
  pivot_longer(cols = c('model.bias', 'direct.bias', 'knn.bias'), 
               values_to = 'bias', names_to = 'Type') %>%
  mutate(Type = ifelse(Type == 'model.bias', 'Model', 
                       ifelse(Type == 'direct.bias', 'Direct', 'kNN')))

bias.sim.plot <- ggplot(data = bias.plot.df.long, aes(fill = Type, y = bias)) +
  geom_boxplot(outliers = FALSE, alpha = 0.5) + 
  facet_wrap(vars(species), scales = 'free_y', nrow = 2, ncol = 5) + 
  scale_fill_viridis_d() + 
  theme_bw(base_size = 18) + 
  labs(y = 'Bias (Mg/ha)') +
  theme(text = element_text(family="LM Roman 10"), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        strip.background = element_blank(), 
        strip.text.x = element_blank(), 
        legend.position = 'bottom') 
# Figure S1 in Supp Information S1 ----------------------------------------        
ggsave(plot = bias.sim.plot, file = 'figures/Figure-S2.png', width = 10, 
       height = 5, units = 'in', bg = 'white')

# Calculate RMSE and summarize --------------------------------------------
model.rmse <- apply(model.bias, 1, function(a) sqrt(mean(a^2))) 
knn.rmse <- apply(knn.bias, 1, function(a) sqrt(mean(a^2))) 
direct.rmse <- apply(direct.bias, 1, function(a) sqrt(mean(a^2))) 
rmse.plot.df <- data.frame(species = factor(rep(1:n.sp, each = 100), 
                                            levels = sp.means, ordered = TRUE),
                           model.rmse = model.rmse,  
                           direct.rmse = direct.rmse, 
                           knn.rmse = knn.rmse)
# Remove NAs situation when a species was never observed in a plot
rmse.plot.df <- rmse.plot.df %>%
  filter(!is.na(direct.rmse))

rmse.plot.df %>%
  group_by(species) %>%
  summarize(model = mean(model.rmse), 
            direct = mean(direct.rmse), 
            knn = mean(knn.rmse))

rmse.plot.df.long <- rmse.plot.df %>%
  pivot_longer(cols = c('model.rmse', 'direct.rmse', 'knn.rmse'), 
               values_to = 'rmse', names_to = 'Type') %>%
  mutate(Type = ifelse(Type == 'model.rmse', 'Model', 
                       ifelse(Type == 'direct.rmse', 'Direct', 'kNN')))

rmse.sim.plot <- ggplot(data = rmse.plot.df.long, aes(fill = Type, y = rmse)) +
  geom_boxplot(outliers = FALSE, alpha = 0.5) + 
  facet_wrap(vars(species), scales = 'free_y', nrow = 2, ncol = 5) + 
  scale_fill_viridis_d() + 
  theme_bw(base_size = 18) + 
  labs(y = 'RMSE (Mg/ha)') +
  theme(text = element_text(family="LM Roman 10"), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        strip.background = element_blank(), 
        strip.text.x = element_blank(), 
        legend.position = 'bottom') 
# Figure S2 in Supp Information S1 ----------------------------------------
ggsave(plot = rmse.sim.plot, file = 'figures/Figure-S3.png', width = 10, 
       height = 5, units = 'in', bg = 'white')
