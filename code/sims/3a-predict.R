# 3a-predict.R: script to predict county-level biomass for each of the simulated
#               datasets.
rm(list = ls())
library(spOccupancy)
library(spAbundance)

# Determine current data set ----------------------------------------------
# NOTE: this in combination with the submit script is used to run multiple
#       of the simulated data sets one at a time. Alternatively, can comment
#       out the line specifying "sim.indx" below and this will allow you to
#       simply run 1 data set at a time on your local computer.
args <- commandArgs(trailingOnly = TRUE)
# Current iteration
sim.indx <- as.numeric(args[1])
# Alternatively, if not running the script from the command line, you can
# specify the individual data sets manually below
# NOTE: uncomment this for running locally.
# sim.indx <- 10
if(length(args) == 0) base::stop('Need to tell spOccupancy the current simulation')
print(sim.indx)

# Directories -------------------------------------------------------------
# NOTE: this is used to determine the directories for reading in and writing out
#       data. You should set out.dir and data.dir depending on where you are
#       running this code.
# Determine what machine you're on and change directories as needed.
machine.name <- Sys.info()['nodename']
if (machine.name == 'pop-os') {
  out.dir <- 'results/sim_results/'
  data.dir <- 'data/sim_data_replicates/'
  full.data.dir <- 'data/'
} else { # Running on NCSU HPC
  out.dir <- '/share/doserlab/jwdoser/DIDF24/results/sim_results/'
  data.dir <- '/share/doserlab/jwdoser/DIDF24/data/sim_data_replicates/'
  full.data.dir <- '/share/doserlab/jwdoser/DIDF24/data/'
}

# Load spOccupancy formatted data -----------------------------------------
load(paste0(data.dir, "sim_", sim.indx, '_data.rda'))

# Functions and CRS -------------------------------------------------------
quants <- function(x){
  quantile(x, prob=c(0.5, 0.025, 0.25, 0.75, 0.975), na.rm = TRUE)
}
my.crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"

# Load model objects ------------------------------------------------------
# Stage 1
load(paste0(out.dir, 'sim-stage-1-50000-samples-', sim.indx, '-replicate-2025-02-23.rda'))
# Make a smaller object with only the necessary components to reduce memory load
out.small.1 <- list()
out.small.1$n.post <- out$n.post
out.small.1$n.chains <- out$n.chains
out.small.1$X <- out$X
out.small.1$y <- out$y
out.small.1$coords <- out$coords
out.small.1$q <- out$q
out.small.1$theta.samples <- out$theta.samples
out.small.1$beta.samples <- out$beta.samples
out.small.1$lambda.samples <- out$lambda.samples
out.small.1$w.samples <- out$w.samples
out.small.1$n.neighbors <- out$n.neighbors
out.small.1$cov.model.indx <- out$cov.model.indx
out.small.1$std.by.sp <- out$std.by.sp
out.small.1$species.means <- out$species.means
out.small.1$species.sds <- out$species.sds
out.small.1$type <- 'NNGP'
out.small.1$psiRE <- out$psiRE
out.small.1$re.level.names <- out$re.level.names
out.small.1$shared.spatial <- out$shared.spatial
class(out.small.1) <- 'sfJSDM'
# Remove the full, much larger object
rm(out)
gc()
# Stage 2
load(paste0(out.dir, 'sim-stage-2-50000-samples-', sim.indx, '-replicate-2025-02-23.rda'))
# Create a smaller object containing only the components needed for prediction.
out.small.2 <- list()
out.small.2$X <- out$X
out.small.2$dist <- out$dist
out.small.2$n.post <- out$n.post
out.small.2$n.chains <- out$n.chains
out.small.2$y <- out$y
out.small.2$coords <- out$coords
out.small.2$q <- out$q
out.small.2$theta.samples <- out$theta.samples
out.small.2$beta.samples <- out$beta.samples
out.small.2$lambda.samples <- out$lambda.samples
out.small.2$n.neighbors <- out$n.neighbors
out.small.2$cov.model.indx <- out$cov.model.indx
out.small.2$re.cols <- out$re.cols
out.small.2$type <- out$type
out.small.2$re.level.names <- out$re.level.names
out.small.2$muRE <- out$muRE
out.small.2$beta.star.samples <- out$beta.star.samples
out.small.2$X.re <- out$X.re
out.small.2$sigma.sq.mu.samples <- out$sigma.sq.mu.samples
out.small.2$tau.sq.samples <- out$tau.sq.samples
out.small.2$w.samples <- out$w.samples
class(out.small.2) <- 'sfMsAbund'
# Remove the full, much larger object
rm(out)
gc()

# Read in prediction object -----------------------------------------------
load(paste0(full.data.dir, "sim_pop_data.rda"))
# Predict one county at a time --------------------------------------------
n.samples <- min(out.small.1$n.post, out.small.2$n.post)
n.counties <- length(unique(county))
unique.counties <- sort(unique(county))
N <- nrow(y)
# Save posterior samples for each species in each county. 
biomass.mean.samples <- array(NA, dim = c(n.samples, N, n.counties))
for (j in 1:n.counties) {
  # Stage 1 ---------------------------
  print(paste0("Currently on county ", j, " out of ", n.counties)) 
  curr.indx <- which(county == unique.counties[j])
  n <- length(curr.indx)
  biomass.samples <- array(NA, dim = c(n.samples, N, n))
  out.pred.1 <- predict(out.small.1, X.1[curr.indx, ], coords.0[curr.indx, ], 
                        n.omp.threads = 1, verbose = FALSE)
  # Stage 2 ---------------------------
  out.pred.2 <- predict(out.small.2, X.2[curr.indx, ], coords.0[curr.indx, ], 
                        n.omp.threads = 1, verbose = FALSE, 
                        z.0.samples = out.pred.1$z.0.samples[1:n.samples, , ])
  # Get biomass samples ---------------
  for (i in 1:N) {
    biomass.samples[, i, ] <- ifelse(out.pred.1$z.0.samples[1:n.samples, i, ] == 1, 
                                     exp(out.pred.2$y.0.samples[, i, ]), 
                                     out.pred.2$y.0.samples[, i, ])
  }
  biomass.mean.samples[, , j] <- apply(biomass.samples, c(1, 2), mean)
  rm(out.pred.1, out.pred.2, biomass.samples)
  gc()
}
# Save results to file
save(biomass.mean.samples, file = paste0(out.dir, 'sim-biomass-', sim.indx, 
                                         '-replicate-sae-county-means.rda'))
