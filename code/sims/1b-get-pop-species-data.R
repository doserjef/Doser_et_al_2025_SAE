# 1b-get-pop-species-data.R: script to simulate the actual species-level
#                            biomass values for the entire simulated population.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)

# Functions ---------------------------------------------------------------
rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p)))){stop("Dimension problem!")}
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

# Set seed for exact results ----------------------------------------------
set.seed(38392)

# Load the covariate data -------------------------------------------------
load('data/sim_pop_covariate_data.rda')
# Number of species to simulate
n.sp <- 10

# Presence/absence simulation ---------------------------------------------
n <- nrow(coords.0)
# Design matrix for stage 1
X.1 <- as.matrix(data.frame(intercept = 1, 
                            tmin = scale(X.0[, 'tmin']),
                            tmin.2 = scale(X.0[, 'tmin'])^2, 
                            tmax = scale(X.0[, 'tmax']),
                            tmax.2 = scale(X.0[, 'tmax'])^2, 
                            ppt = scale(X.0[, 'ppt']),
                            ppt.2 = scale(X.0[, 'ppt'])^2))
# Average occupancy coefficients                            
beta.mean <- c(-2, 0, -0.5, -0.75, -0.7, 0.6, -0.4)
# Variances in occupancy coefficients across the species
tau.sq.beta <- c(4, 5, 0.3, 1.5, 0.5, 0.5, 0.1)
# Simulate species-level parameters from the means and variances
p.abund <- length(beta.mean)
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
for (i in 1:n.sp) {
  beta[i, ] <- rnorm(p.abund, beta.mean, sqrt(tau.sq.beta))
}
# Simulate spatial factors ------------
# A range of spatial decay parameters (all pretty broad-scale)
phi <- c(3 / 700, 3 / 400, 3 / 300, 3 / 200)
n.factors <- 4
# Matrix of spatial locations
w.star <- matrix(0, n.sp, n)
w <- matrix(NA, n.factors, n)
# Factor loadings matrix
lambda <- matrix(rnorm(n.sp * n.factors, 0, 1), n.sp, n.factors)
# Set diagonals to 1
diag(lambda) <- 1
# Set upper tri to 0
lambda[upper.tri(lambda)] <- 0
# Generate the spatial factors from a Gaussian Process. 
for (ll in 1:n.factors) {
  print(ll)
  Sigma <- spBayes::mkSpCov(as.matrix(coords.0), as.matrix(1), as.matrix(0),
                            phi[ll], 'exponential')
  w[ll, ] <- rmvn(1, rep(0, n), Sigma)
  rm(Sigma)
  gc()
}
for (i in 1:n) {
  w.star[, i] <- lambda %*% w[, i]
}
# Generate the latent process ---------
psi <- logit.inv(beta %*% t(X.1) + w.star)
z <- matrix(NA, nrow = n.sp, ncol = n)
for (j in 1:n.sp) {
  z[j, ] <- rbinom(n = n, size = 1, prob = psi[j, ])
}

# Generate biomass on the log scale ---------------------------------------
# Design matrix for stage 2
set.seed(111)
X.2 <- as.matrix(data.frame(intercept = 1, 
                            tcc = scale(X.0[, 'tcc']),
                            vpd = scale(X.0[, 'vpd']),
                            ppt = scale(X.0[, 'ppt']), 
                            elev = scale(X.0[, 'elev']),
                            elev.2 = scale(X.0[, 'elev'])^2))
# Notice that here we simulate species-specific coefficients from a hierarchical
# distribution, but this is not what is done when fitting the model. The individual
# species estimates for stage 2 are estimated independently. 
# Average community coefficients                            
beta.mean <- c(3, 0.2, -0.4, 0, 0.4, -0.4)
# Variances across the species
tau.sq.beta <- c(0.75, 0.4, 0.4, 0.3, 0.2, 0.2)
# Simulate species-level parameters from the means and variances
p.abund <- length(beta.mean)
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
for (i in 1:n.sp) {
  beta[i, ] <- rnorm(p.abund, beta.mean, sqrt(tau.sq.beta))
}
# Simulate spatial factors ------------
# A range of spatial decay parameters
phi <- c(3 / 500, 3 / 200, 3 / 100, 3 / 60)
n.factors <- 4
# Matrix of spatial locations
w.star <- matrix(0, n.sp, n)
w <- matrix(NA, n.factors, n)
lambda <- matrix(rnorm(n.sp * n.factors, 0, 0.2), n.sp, n.factors)
# Set diagonals to 1
diag(lambda) <- 1
# Set upper tri to 0
lambda[upper.tri(lambda)] <- 0
for (ll in 1:n.factors) {
  print(ll)
  Sigma <- spBayes::mkSpCov(as.matrix(coords.0), as.matrix(1), as.matrix(0),
                            phi[ll], 'exponential')
  w[ll, ] <- rmvn(1, rep(0, n), Sigma)
  rm(Sigma)
  gc()
}
for (i in 1:n) {
  w.star[, i] <- lambda %*% w[, i]
}
# Generate the latent process ---------
mu <- beta %*% t(X.2) + w.star
# Simulate a substantial amount of unexplained variation. 
tau.sq <- runif(n.sp, 1, 3)
# Biomass on th elog scale. 
y <- matrix(NA, nrow = n.sp, ncol = n)
for (j in 1:n.sp) {
  y[j, ] <- rnorm(n, mu[j, ], sqrt(tau.sq[j]))
}

# y is on the log scale. Actual biomass is below
biomass <- exp(y) 

# Convert values in y and biomass to 0 when z == 0
y <- ifelse(z == 0, 0, y)
biomass <- ifelse(z == 0, 0, biomass)

# Quickly plot the biomass of one of the species just for fun -------------
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
my.crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"
usa.county <- st_as_sf(maps::map("county", fill = TRUE, plot = FALSE))
usa.county <- usa.county %>%
  st_transform(crs = my.crs)
nc.county <- usa.county %>%
  separate_wider_delim(ID, delim = ',', names = c('state', 'county')) %>%
  dplyr::filter(state %in% c('north carolina')) %>%
  st_as_sf()

coords.sf <- st_as_sf(coords.0,
                      coords = c('X', 'Y'),
                      crs = my.crs)

# Temporary plotting code
i <- 5
coords.sf$y <- y[i, ]
ggplot() + 
  geom_sf(data = coords.sf, aes(col = y)) + 
  scale_color_viridis_c() +
  theme_bw()


# Determine county that each point falls into -----------------------------
indx.by.county <- st_contains(nc.county, coords.sf)
# Fill in county vector with the appropriate information
n.counties <- nrow(nc.county)
n.plots <- nrow(coords.0)
county <- vector(mode = 'character', length = n.plots)
for (i in 1:n.counties) {
  county[indx.by.county[[i]]] <- nc.county$county[i]
}
unique.counties <- unique(county)

# Save data to hard drive -------------------------------------------------
save(y, biomass, mu, z, psi, coords.0, X.1, X.2, unique.counties, 
     county, file = 'data/sim_pop_data.rda')


