# 6-summary.R: script to summarize results from the Southeastern US case study.
# Author: Jeffrey W. Doser
rm(list = ls())
library(spAbundance)
library(spOccupancy)
library(coda)
library(tidyverse)
library(sf)
library(patchwork)
# For the simpleCap function
library(clinUtils)

# Directories -------------------------------------------------------------
out.dir <- 'results/'
data.dir <- 'data/'

# Read in and prep spatial data -------------------------------------------
load(paste0(data.dir, "se_bio_stage_1_data.rda"))
load(paste0(data.dir, "se_bio_stage_2_data.rda"))
load(paste0(data.dir, "se-prediction-data.rda"))
# Number of counties and their names in the order of the resulting estimates.
counties <- unique(county)
n.counties <- length(counties)
# Species names
sp.names <- dimnames(data.list.1$y)[[1]]

my.crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"
coords.sf <- st_as_sf(data.frame(data.list.2$coords),
                      coords = c('X', 'Y'),
                      crs = my.crs)

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
# Determine number of points in each county
indx.by.county <- st_contains(se.county, coords.sf)
se.county$n.plots <- sapply(indx.by.county, length)

# Generate histograms showing skewness in data distributions --------------
y.bio <- ifelse(data.list.2$y == 0, data.list.2$y, exp(data.list.2$y))
# Determine which species has the most biomass
# Loblolly pine, white oak, yellow poplar, sweetgum
sort(apply(y.bio, 1, sum))
# Convert biomass from short tons per acre to Mg.ha
y.bio <- y.bio * 2.2417
# Plot showing
y.curr.sp <- y.bio[which(sp.names %in% c('loblolly pine', 'white oak',
                                         'yellow-poplar', 'sweetgum')), ]
J <- ncol(y.curr.sp)
plot.df <- data.frame(bio = c(y.curr.sp),
                      species = rep(c('Loblolly pine', 'White oak',
                                      'Yellow poplar', 'Sweetgum'), times = J))
ggplot(data = plot.df, aes(x = bio, fill = species)) +
  geom_histogram() +
  facet_wrap(vars(species), scales = 'free') +
  scale_fill_viridis_d() +
  theme_bw(base_size = 18) +
  guides(fill = 'none') +
  labs(x = 'Biomass (Mg/ha)', y = 'Number of plots') +
  theme(text = element_text(family="LM Roman 10"))
ggsave(file = 'figures/Figure-S1.png', width = 10, height = 7, units = 'in',
       bg = 'white')

# Explicit comparisons between model and direct estimates -----------------
load('results/se_direct_estimates.rda')
load('results/se-biomass-sae-county-means.rda')
load('results/k-means-sae-estimates.rda')
# Biomass medians from model-based approach
sae.meds <- apply(biomass.mean.samples, c(2, 3), median)
sae.ci.low <- apply(biomass.mean.samples, c(2, 3), quantile, 0.025)
sae.ci.high <- apply(biomass.mean.samples, c(2, 3), quantile, 0.975)
sae.ci.width <- sae.ci.high - sae.ci.low
rownames(sae.meds) <- sp.names
colnames(sae.meds) <- counties
rownames(sae.ci.low) <- sp.names
colnames(sae.ci.low) <- counties
rownames(sae.ci.high) <- sp.names
colnames(sae.ci.high) <- counties

# Convert sae.meds into long format
model.ests <- as.data.frame(t(sae.meds))
model.ests <- model.ests %>%
  mutate(county = counties) %>%
  pivot_longer(data = ., cols = -county, names_to = 'species', values_to = 'model_ests')

model.low.ests <- as.data.frame(t(sae.ci.low))
model.low.ests <- model.low.ests %>%
  mutate(county = counties) %>%
  pivot_longer(data = ., cols = -county, names_to = 'species', values_to = 'model_ci_low')

model.high.ests <- as.data.frame(t(sae.ci.high))
model.high.ests <- model.high.ests %>%
  mutate(county = counties) %>%
  pivot_longer(data = ., cols = -county, names_to = 'species', values_to = 'model_ci_high')

# Join all model-based means, cis in one
model.ests <- left_join(model.ests, model.low.ests, by = c('species', 'county'))
model.ests <- left_join(model.ests, model.high.ests, by = c('species', 'county'))

direct.ests <- direct.ests %>%
  mutate(county = as.character(county)) %>%
  rename(direct_ests = mu) 
full.df <- left_join(model.ests, direct.ests, by = c('species', 'county'))
# Add in the k-means estimates as well
k.means.long.df <- k.means.county.ests %>%
  pivot_longer(data = ., cols = -county, names_to = 'species', values_to = 'k_mean_ests')
full.df <- left_join(full.df, k.means.long.df, by = c('species', 'county'))

# Convert from short tons per acre to Megagrams per hectare
full.df <- full.df %>%
  mutate(model_ests = model_ests * 2.2417, 
         model_ci_low = model_ci_low * 2.2417,
         model_ci_high = model_ci_high * 2.2417,
         direct_ests = direct_ests * 2.2417, 
         k_mean_ests = k_mean_ests * 2.2417)

# Species-level correlations
sp.cor.df <- full.df %>%
  group_by(species) %>%
  summarize(cor.direct = cor(model_ests, direct_ests, use = 'pairwise.complete.obs')) %>%
  arrange(desc(cor.direct))

full.df <- full.df %>% 
  mutate(species.plot = factor(species, levels = sp.cor.df$species, order = TRUE))


# Figure 2 ----------------------------
fig.2 <- full.df %>%
  ggplot(aes(x = model_ests, y = direct_ests)) + 
    geom_point() + 
    facet_wrap(vars(species.plot), scales = 'free', nrow = 4, ncol = 5) + 
    geom_abline(slope = 1, intercept = 0) + 
    theme_light(base_size = 18) + 
    theme(text = element_text(family="LM Roman 10"), 
          strip.text.y = element_text(color = 'black'),
          strip.text.x = element_text(color = 'black'),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank()) + 
    labs(x = 'Model estimates (Mg/ha)', y = 'Direct estimates (Mg/ha)')
ggsave(plot = fig.2, file = 'figures/Figure-2.png', width = 14,
       height = 8, units = 'in', bg = 'white')

# Uncertainty comparisons -------------------------------------------------
sae.sd <- apply(biomass.mean.samples, c(2, 3), sd)               
rownames(sae.sd) <- sp.names
colnames(sae.sd) <- counties

# Convert sae.meds into long format
model.ests <- as.data.frame(t(sae.sd))
model.ests <- model.ests %>%
  mutate(county = counties) %>%
  pivot_longer(data = ., cols = -county, names_to = 'species', values_to = 'model_ests')
direct.se.ests <- direct.se.ests %>%
  mutate(county = as.character(county)) %>%
  rename(direct_ests = mu) 
full.df.se <- left_join(model.ests, direct.se.ests, by = c('species', 'county'))

# Convert from short tons per acre to Megagrams per hectare
full.df.se <- full.df.se %>%
  mutate(model_ests = model_ests * 2.2417, 
         direct_ests = direct_ests * 2.2417)

# Look at relative values 
full.df.se <- full.df.se %>%
  mutate(relative_model_sd = model_ests/full.df$model_ests, 
         relative_direct_se = direct_ests/full.df$direct_ests)

# Combine SE values with means all in one
complete.df <- data.frame(county = full.df$county, 
                          species = full.df$species, 
                          model_ests = full.df$model_ests, 
                          model_ci_low = full.df$model_ci_low,
                          model_ci_high = full.df$model_ci_high,
                          direct_ests = full.df$direct_ests,
                          k_mean_ests = full.df$k_mean_ests, 
                          model_se = full.df.se$model_ests,
                          direct_se = full.df.se$direct_ests,
                          model_cv = full.df.se$relative_model_sd, 
                          direct_cv = full.df.se$relative_direct_se)

# Overall percentage of mean improvement
complete.df %>%
  filter(direct_ests != 0) %>%
  mutate(indicator = model_cv < direct_cv) %>%
  summarize(mean.improve = mean(indicator, na.rm = TRUE) * 100) %>%
  arrange(mean.improve)

plot.cv.df <- complete.df %>%
  mutate(indicator = model_cv < direct_cv) %>%
  group_by(species) %>%
  summarize(mean.improve = mean(indicator, na.rm = TRUE) * 100) %>%
  arrange(mean.improve) %>%
  mutate(species = factor(species, order = TRUE, levels = species))

# Figure 4 ----------------------------
figure.4 <- plot.cv.df %>%
  ggplot(aes(x = mean.improve, y = species)) + 
  geom_point(size = 3, pch = 19) + 
  geom_vline(xintercept = 50, lty = 2) +
  theme_bw(base_size = 15) +                        
  labs(x = 'Percentage of counties with lower model-based CV', 
       y = 'Species') + 
  theme(text = element_text(family="LM Roman 10"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 12))
ggsave(plot = figure.4, file = 'figures/Figure-4.png', width = 8,
       height = 6, units = 'in', bg = 'white')


# Correlations between all three estimators -------------------------------
# By species 
complete.df %>% 
  group_by(species) %>%
  summarize(design.cor = cor(model_ests, direct_ests, use = 'pairwise.complete.obs'), 
            k.cor = cor(model_ests, k_mean_ests, use = 'pairwise.complete.obs')) %>%
  arrange(desc(design.cor))

# Overall
complete.df %>% 
  group_by(species) %>%
  summarize(design.cor = cor(model_ests, direct_ests, use = 'pairwise.complete.obs'), 
            k.cor = cor(model_ests, k_mean_ests, use = 'pairwise.complete.obs')) %>%
  arrange(desc(design.cor)) %>% 
  summarize(design = mean(design.cor), 
            kNN = mean(k.cor))

# Overall averages
complete.df %>%
  group_by(species) %>%
  summarize(design.mean = mean(direct_ests, na.rm = TRUE), 
            model.mean = mean(model_ests, na.rm = TRUE),
            k.mean = mean(k_mean_ests, na.rm = TRUE))

# Add relative efficiency for later plots
complete.df <- complete.df %>%
  mutate(rel_eff = direct_cv / model_cv)

# Figure 3 (example map of species-specific estimates ---------------------
# Most common species
indx <- which(sp.names == 'loblolly pine')
# Note the conversion from short tons/acre to megagrams per ha
plot.df <- data.frame(mu_mg_ha = sae.meds[indx, ] * 2.2417,
                      ci_mg_ha = sae.ci.width[indx, ] * 2.2417,
                      county = dimnames(sae.meds)[[2]])
plot.sf <- inner_join(se.county, plot.df, by = 'county')
lob.mean.plot <- ggplot() +
  geom_sf(data = plot.sf, aes(fill = mu_mg_ha)) +
  scale_fill_viridis_c(na.value = 'grey') +
  labs(x = 'Longitude', y = 'Latitude', fill = 'Aboveground\nBiomass (Mg/ha)',
       title = '(a) Loblolly pine median biomass') +
  theme_bw(base_size = 10) +
  theme(text = element_text(family="LM Roman 10"),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        # axis.text.y = element_text(size = 18),
        legend.key.size = unit(0.5, 'cm'),
        legend.direction="horizontal",
        legend.position = 'inside',
        legend.position.inside=c(0.25,0.85),
        legend.background = element_rect(fill = NA))
lob.ci.plot <- ggplot() +
  geom_sf(data = plot.sf, aes(fill = ci_mg_ha)) +
  scale_fill_viridis_c(na.value = 'grey') +
  labs(x = 'Longitude', y = 'Latitude', fill = 'Aboveground\nBiomass (Mg/ha)',
       title = '(b) Loblolly pine 95% CI width') +
  theme_bw(base_size = 10) +
  theme(text = element_text(family="LM Roman 10"),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        # axis.text.y = element_text(size = 18),
        legend.key.size = unit(0.5, 'cm'),
        legend.direction="horizontal",
        legend.position = 'inside',
        legend.position.inside=c(0.25,0.85),
        legend.background = element_rect(fill = NA))
# Least abundant species
indx <- which(sp.names == 'American beech')
# Note the conversion from short tons/acre to megagrams per ha
plot.df <- data.frame(mu_mg_ha = sae.meds[indx, ] * 2.2417,
                      ci_mg_ha = sae.ci.width[indx, ] * 2.2417,
                      county = dimnames(sae.meds)[[2]])
plot.sf <- inner_join(se.county, plot.df, by = 'county')
beech.mean.plot <- ggplot() +
  geom_sf(data = plot.sf, aes(fill = mu_mg_ha)) +
  scale_fill_viridis_c(na.value = 'grey') +
  labs(x = 'Longitude', y = 'Latitude', fill = 'Aboveground\nBiomass (Mg/ha)',
       title = '(c) American beech median biomass') +
  theme_bw(base_size = 10) +
  theme(text = element_text(family="LM Roman 10"),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        # axis.text.y = element_text(size = 18),
        legend.key.size = unit(0.5, 'cm'),
        legend.direction="horizontal",
        legend.position = 'inside',
        legend.position.inside=c(0.25,0.85),
        legend.background = element_rect(fill = NA))
beech.ci.plot <- ggplot() +
  geom_sf(data = plot.sf, aes(fill = ci_mg_ha)) +
  scale_fill_viridis_c(na.value = 'grey') +
  labs(x = 'Longitude', y = 'Latitude', fill = 'Aboveground\nBiomass (Mg/ha)',
       title = '(d) American beech 95% CI width') +
  theme_bw(base_size = 10) +
  theme(text = element_text(family="LM Roman 10"),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        # axis.text.y = element_text(size = 18),
        legend.key.size = unit(0.5, 'cm'),
        legend.direction="horizontal",
        legend.position = 'inside',
        legend.position.inside=c(0.25,0.85),
        legend.background = element_rect(fill = NA))
full.plot <- lob.mean.plot + lob.ci.plot + beech.mean.plot + beech.ci.plot
ggsave(plot = full.plot, file = 'figures/Figure-3.png', 
       width = 10, height = 7, units = 'in', bg = 'white')

# Figure 5 (example relative efficiency plots) ----------------------------
# Two species with the worst gains in efficiency
longleaf.pine.plot.sf <- complete.df %>%
  filter(species == 'longleaf pine')
longleaf.pine.plot.sf <- inner_join(se.county, longleaf.pine.plot.sf, by = 'county')

longleaf.pine.rel.eff.plot <- ggplot() +
  geom_sf(data = longleaf.pine.plot.sf, aes(fill = rel_eff), col = 'lightgray') +
  scale_fill_gradient2(midpoint = 1, high = '#2166AC', mid = 'white', low = '#B2182B',
                       na.value = 'grey') +
  labs(x = 'Longitude', y = 'Latitude', fill = 'Relative\nefficiency',
       title = '(a) Longleaf pine') +
  theme_bw(base_size = 10) +
  theme(text = element_text(family="LM Roman 10"),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        # axis.text.y = element_text(size = 18),
        legend.key.size = unit(0.5, 'cm'),
        legend.direction="horizontal",
        legend.position = 'inside',
        legend.position.inside=c(0.25,0.85),
        legend.background = element_rect(fill = NA))

loblolly.pine.plot.sf <- complete.df %>%
  filter(species == 'loblolly pine')
loblolly.pine.plot.sf <- inner_join(se.county, loblolly.pine.plot.sf, by = 'county')

loblolly.pine.rel.eff.plot <- ggplot() +
  geom_sf(data = loblolly.pine.plot.sf, aes(fill = rel_eff), col = 'lightgray') +
  scale_fill_gradient2(midpoint = 1, high = '#2166AC', mid = 'white', low = '#B2182B',
                       na.value = 'grey') +
  labs(x = 'Longitude', y = 'Latitude', fill = 'Relative\nefficiency',
       title = '(b) Loblolly pine') +
  theme_bw(base_size = 10) +
  theme(text = element_text(family="LM Roman 10"),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        # axis.text.y = element_text(size = 18),
        legend.key.size = unit(0.5, 'cm'),
        legend.direction="horizontal",
        legend.position = 'inside',
        legend.position.inside=c(0.25,0.85),
        legend.background = element_rect(fill = NA))

# Two species with the best gains in efficiency
mock.plot.sf <- complete.df %>%
  filter(species == 'mockernut hickory')
mock.plot.sf <- inner_join(se.county, mock.plot.sf, by = 'county')

mock.rel.eff.plot <- ggplot() +
  geom_sf(data = mock.plot.sf, aes(fill = rel_eff), col = 'lightgray') +
  scale_fill_gradient2(midpoint = 1, high = '#2166AC', mid = 'white', low = '#B2182B',
                       na.value = 'grey') +
  labs(x = 'Longitude', y = 'Latitude', fill = 'Relative\nefficiency',
       title = '(c) Mockernut hickory') +
  theme_bw(base_size = 10) +
  theme(text = element_text(family="LM Roman 10"),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        # axis.text.y = element_text(size = 18),
        legend.key.size = unit(0.5, 'cm'),
        legend.direction="horizontal",
        legend.position = 'inside',
        legend.position.inside=c(0.25,0.85),
        legend.background = element_rect(fill = NA))

ypop.plot.sf <- complete.df %>%
  filter(species == 'yellow-poplar')
ypop.plot.sf <- inner_join(se.county, ypop.plot.sf, by = 'county')

ypop.rel.eff.plot <- ggplot() +
  geom_sf(data = ypop.plot.sf, aes(fill = rel_eff), col = 'lightgray') +
  scale_fill_gradient2(midpoint = 1, high = '#2166AC', mid = 'white', low = '#B2182B',
                       na.value = 'grey') +
  labs(x = 'Longitude', y = 'Latitude', fill = 'Relative\nefficiency',
       title = '(d) Yellow-poplar') +
  theme_bw(base_size = 10) +
  theme(text = element_text(family="LM Roman 10"),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        # axis.text.y = element_text(size = 18),
        legend.key.size = unit(0.5, 'cm'),
        legend.direction="horizontal",
        legend.position = 'inside',
        legend.position.inside=c(0.25,0.85),
        legend.background = element_rect(fill = NA))

full.plot <- longleaf.pine.rel.eff.plot + loblolly.pine.rel.eff.plot + 
             mock.rel.eff.plot + ypop.rel.eff.plot
ggsave(plot = full.plot, file = 'figures/Figure-5.png', 
       width = 10, height = 7, units = 'in', bg = 'white')

# Relative efficiency summarized across all species -----------------------
rel.plot.all.df <- complete.df %>%
  mutate(improved = ifelse(rel_eff > 1, 1, 0)) %>%
  group_by(county) %>%
  summarize(prop.improved = mean(improved, na.rm = TRUE))
plot.sf <- inner_join(se.county, rel.plot.all.df, by = 'county')
fig.6 <- ggplot() +
    geom_sf(data = plot.sf, aes(fill = prop.improved), col = NA) +
    scale_fill_gradient2(midpoint = 0.5, high = '#2166AC', mid = 'white', low = '#B2182B',
                         na.value = 'grey') +
    labs(x = 'Longitude', y = 'Latitude', fill = 'Proportion of species\nwith improved\nprecision') +
    theme_bw(base_size = 10) +
    theme(text = element_text(family="LM Roman 10"),
          axis.title.x=element_blank(), axis.title.y=element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          # axis.text.y = element_text(size = 18),
          legend.key.size = unit(0.5, 'cm'),
          legend.direction="horizontal",
          legend.position = 'inside',
          legend.position.inside=c(0.25,0.85),
          legend.background = element_rect(fill = NA))
ggsave(plot = fig.6, file = 'figures/Figure-6.png', 
       width = 6, height = 4, units = 'in', bg = 'white')


# Spatial maps for each species -------------------------------------------
# These are the plots shown in Supplemental Information S2
# Relative efficiency plots -----------
for (i in 1:length(sp.names)) {
  sp.plot.sf <- complete.df %>%
    filter(species == sp.names[i])
  sp.plot.sf <- inner_join(se.county, sp.plot.sf, by = 'county')
  
  sp.rel.eff.plot <- ggplot() +
    geom_sf(data = sp.plot.sf, aes(fill = rel_eff), col = 'lightgray') +
    scale_fill_gradient2(midpoint = 1, high = '#2166AC', mid = 'white', low = '#B2182B',
                         na.value = 'grey') +
    labs(x = 'Longitude', y = 'Latitude', fill = 'Relative\nefficiency',
           title = simpleCap(str_replace_all(sp.names[i], '-', ' '))) +
    theme_bw(base_size = 10) +
    theme(text = element_text(family="LM Roman 10"),
          axis.title.x=element_blank(), axis.title.y=element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          # axis.text.y = element_text(size = 18),
          legend.key.size = unit(0.5, 'cm'),
          legend.direction="horizontal",
          legend.position = 'inside',
          legend.position.inside=c(0.25,0.85),
          legend.background = element_rect(fill = NA))
    ggsave(plot = sp.rel.eff.plot, file = paste('figures/species-maps/', str_replace_all(sp.names[i], " ", "-"),
	      '-rel-efficiency.png', sep = ''), height = 6, width = 6, units = 'in'  )
}

# Mean and 95% credible interval ------
for (i in 1:length(sp.names)) {
  sp.curr <- sp.names[i]
  indx <- which(sp.names == sp.curr)
  # Note the conversion from short tons/acre to megagrams per ha
  plot.df <- data.frame(mu_mg_ha = sae.meds[indx, ] * 2.2417,
                        ci_mg_ha = sae.ci.width[indx, ] * 2.2417,
                        county = dimnames(sae.meds)[[2]])
  plot.sf <- inner_join(se.county, plot.df, by = 'county')
  mean.plot <- ggplot() +
    geom_sf(data = plot.sf, aes(fill = mu_mg_ha)) +
    scale_fill_viridis_c(na.value = 'grey') +
    labs(x = 'Longitude', y = 'Latitude', fill = 'Aboveground\nBiomass (Mg/ha)',
         title = paste0(simpleCap(str_replace_all(sp.names[i], '-', ' ')), ' mean')) +
    theme_bw(base_size = 10) +
    theme(text = element_text(family="LM Roman 10"),
          axis.title.x=element_blank(), axis.title.y=element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          # axis.text.y = element_text(size = 18),
          legend.key.size = unit(0.5, 'cm'),
          legend.direction="horizontal",
          legend.position = 'inside',
          legend.position.inside=c(0.25,0.85),
          legend.background = element_rect(fill = NA))
  ci.plot <- ggplot() +
    geom_sf(data = plot.sf, aes(fill = ci_mg_ha)) +
    scale_fill_viridis_c(na.value = 'grey') +
    labs(x = 'Longitude', y = 'Latitude', fill = 'Aboveground\nBiomass (Mg/ha)',
         title = paste0(simpleCap(str_replace_all(sp.names[i], '-', ' ')), ' 95% CI Width')) +
    theme_bw(base_size = 10) +
    theme(text = element_text(family="LM Roman 10"),
          axis.title.x=element_blank(), axis.title.y=element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          # axis.text.y = element_text(size = 18),
          legend.key.size = unit(0.5, 'cm'),
          legend.direction="horizontal",
          legend.position = 'inside',
          legend.position.inside=c(0.25,0.85),
          legend.background = element_rect(fill = NA))
    full.plot <- mean.plot + ci.plot
    ggsave(plot = full.plot, file = paste('figures/species-maps/', str_replace_all(sp.names[i], " ", "-"),
	      '-county-ests.png', sep = ''), height = 6, width = 12, units = 'in'  )
}

# Direct estimate maps for each species -----------------------------------
for (i in 1:length(sp.names)) {
  sp.plot.sf <- complete.df %>%
    filter(species == sp.names[i])
  sp.plot.sf <- inner_join(se.county, sp.plot.sf, by = 'county')
  mean.plot <- ggplot() +
    geom_sf(data = sp.plot.sf, aes(fill = direct_ests)) +
    scale_fill_viridis_c(na.value = 'grey') +
    labs(x = 'Longitude', y = 'Latitude', fill = 'Aboveground\nBiomass (Mg/ha)',
         title = paste0(simpleCap(str_replace_all(sp.names[i], '-', ' ')), ' mean')) +
    theme_bw(base_size = 10) +
    theme(text = element_text(family="LM Roman 10"),
          axis.title.x=element_blank(), axis.title.y=element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          # axis.text.y = element_text(size = 18),
          legend.key.size = unit(0.5, 'cm'),
          legend.direction="horizontal",
          legend.position = 'inside',
          legend.position.inside=c(0.25,0.85),
          legend.background = element_rect(fill = NA))
    ggsave(plot = mean.plot, file = paste('figures/species-maps/direct/', 
           str_replace_all(sp.names[i], " ", "-"),
           '-direct-ests.png', sep = ''), height = 6, width = 6, units = 'in'  )
}
