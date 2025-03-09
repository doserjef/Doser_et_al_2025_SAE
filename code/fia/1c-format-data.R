# 1c-format-data.R: script to format the FIA data into the necessary format
#                   log-normal hurdle model in spOccupancy and spAbundance.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)

# Load data ---------------------------------------------------------------
# FIA data and coordinates
load("data/se_fia_bio_data.rda")
# Covariates
load("data/covariate-data.rda")

# Get coordinates in AEA --------------------------------------------------
my.crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"
coords.sf.aea <- final.wide.dat %>%
  select(geometry) %>%
  st_transform(crs = my.crs)
coords <- coords.sf.aea %>% st_coordinates()

# Extract data for formatting with the models -----------------------------
# Number of sites
n <- nrow(coords)
# Get data in y matrix format
y.df <- final.wide.dat %>%
  select(-pltID, -YEAR, -PLT_CN) %>%
  st_drop_geometry()
sp.names <- colnames(y.df)
# Biomass matrix
y <- y.df %>% as.matrix() %>% t()

# Presence-absence matrix
z <- ifelse(y > 0, 1, y)

# Remove sites with missing covariate values ------------------------------
bad.indx <- c(which(apply(X, 1, function(a) sum(is.na(a)) > 0)),
              which(X$tcc > 100))
coords <- coords[-bad.indx, ]
X <- X[-bad.indx, ]
y <- y[, -bad.indx]
z <- z[, -bad.indx]

# Form the data lists -----------------------------------------------------
# Data list for Stage 1
data.list.1 <- list(y = z,
                    coords = coords,
                    covs = X)

# Data list for Stage 2. Note the biomass data are supplied in the log format.
y.logged <- ifelse(y == 0, y, log(y))
data.list.2 <- list(y = y.logged,
                    coords = coords,
                    covs = X,
                    z = z)

str(data.list.1)
str(data.list.2)

# Reformat species ordering to add in mixing and convergence --------------
# Reorder species to aid in mixing of MCMC chains.
# Using 5 factors
sp.names <- dimnames(data.list.1$y)[[1]]
# Longleaf, northern red oak, loblolly, laurel oak, chestnut oak
start.sp <- c('longleaf pine', 'northern red oak', 'loblolly pine', 'laurel oak',
              'chestnut oak')
# Other species codes
indices <- rep(NA, length(start.sp))
for (i in 1:length(start.sp)) {
  indices[i] <- which(sp.names == start.sp[i])
}
indices.other <- 1:nrow(data.list.1$y)
indices.other <- indices.other[-indices]
# Create the ordered y data frame
data.list.1$y <- data.list.1$y[c(indices, indices.other), ]
data.list.2$y <- data.list.2$y[c(indices, indices.other), ]
data.list.2$z <- data.list.2$z[c(indices, indices.other), ]

# Save data lists to hard drive -------------------------------------------
save(data.list.1, file = 'data/se_bio_stage_1_data.rda')
save(data.list.2, file = 'data/se_bio_stage_2_data.rda')

