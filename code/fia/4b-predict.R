# 4b-predict.R: script to predict county-level biomass for each of the 20 species
#               across the southeastern US. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)
library(spAbundance)

# NOTE: this script will not run using the files available on GitHub, because
# the results files are too large to include on GitHub. Instead, if you run
# the previous main model running scripts, you can read in the output file into this file
# to use in subsequent predictions.

# Specify directory -------------------------------------------------------
# NOTE: this is used to determine the directories for reading in and writing out
#       data. You should set out.dir and data.dir depending on where you are 
#       running this code. 
machine.name <- Sys.info()['nodename']
if (machine.name == 'pop-os') {
  out.dir <- 'results/'
  data.dir <- 'data/'
} else { # Running on NCSU HPC
  out.dir <- '/share/doserlab/jwdoser/DIDF24/results/'
  data.dir <- '/share/doserlab/jwdoser/DIDF24/data/'
}

# Functions and CRS -------------------------------------------------------
quants <- function(x){
  quantile(x, prob=c(0.5, 0.025, 0.25, 0.75, 0.975), na.rm = TRUE)
}
my.crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"

# Load Stage 1 data -------------------------------------------------------
load(paste0(data.dir, "se_bio_stage_1_data.rda"))
tmin.mean <- mean(data.list.1$covs$tmin)
tmin.sd <- sd(data.list.1$covs$tmin)
tmax.mean <- mean(data.list.1$covs$tmax)
tmax.sd <- sd(data.list.1$covs$tmax)
ppt.mean <- mean(data.list.1$covs$ppt)
ppt.sd <- sd(data.list.1$covs$ppt)

# Load Stage 2 data -------------------------------------------------------
load(paste0(data.dir, "se_bio_stage_2_data.rda"))

tcc.mean <- mean(data.list.2$covs$tcc)
tcc.sd <- sd(data.list.2$covs$tcc)
vpd.mean <- mean(data.list.2$covs$vpd)
vpd.sd <- sd(data.list.2$covs$vpd)
ppt.mean <- mean(data.list.2$covs$ppt)
ppt.sd <- sd(data.list.2$covs$ppt)
elev.mean <- mean(data.list.2$covs$elev)
elev.sd <- sd(data.list.2$covs$elev)

# Load resulting model objects --------------------------------------------
# Stage 1
load(paste0(out.dir, 'se-stage-1-2e+05-samples-chain-1-5-factors-2025-02-23.rda'))
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
load(paste0(out.dir, 'se-stage-2-1e+05-samples-5-factors-2025-02-19.rda'))
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
load(paste0(data.dir, "se-prediction-data.rda"))
tmin.pred <- (X.0[, 'tmin'] - tmin.mean) / tmin.sd
tmax.pred <- (X.0[, 'tmax'] - tmax.mean) / tmax.sd
ppt.pred <- (X.0[, 'ppt'] - ppt.mean) / ppt.sd
elev.pred <- (X.0[, 'elev'] - elev.mean) / elev.sd
tcc.pred <- (X.0[, 'tcc'] - tcc.mean) / tcc.sd
vpd.pred <- (X.0[, 'vpd'] - vpd.mean) / vpd.sd
X.0.stage.1 <- cbind(1, tmin.pred, tmin.pred^2, tmax.pred, tmax.pred^2, 
                     ppt.pred, ppt.pred^2)
dimnames(X.0.stage.1)[[2]] <- dimnames(out.small.1$X)[[2]]
X.0.stage.2 <- cbind(1, tcc.pred, vpd.pred, ppt.pred, elev.pred, 
                     elev.pred^2)
dimnames(X.0.stage.2)[[2]] <- c(dimnames(out.small.2$X)[[2]])
# Predict one county at a time --------------------------------------------
n.samples <- min(out.small.1$n.post, out.small.2$n.post)
n.counties <- length(unique(county))
unique.counties <- unique(county)
N <- nrow(data.list.1$y)
# Save posterior samples for each species in each county. 
biomass.mean.samples <- array(NA, dim = c(n.samples, N, n.counties))
for (j in 1:n.counties) {
  # Stage 1 ---------------------------
  print(paste0("Currently on county ", j, " out of ", n.counties)) 
  curr.indx <- which(county == unique.counties[j])
  J.0 <- length(curr.indx)
  biomass.samples <- array(NA, dim = c(n.samples, N, J.0))
  out.pred.1 <- predict(out.small.1, X.0.stage.1[curr.indx, ], coords.0[curr.indx, ], 
                        n.omp.threads = 5, verbose = FALSE)
  # Stage 2 ---------------------------
  out.pred.2 <- predict(out.small.2, X.0.stage.2[curr.indx, ], coords.0[curr.indx, ], 
                        n.omp.threads = 5, verbose = FALSE, 
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
# Save results to hard drive ----------------------------------------------
save(biomass.mean.samples, file = paste0(out.dir, 'se-biomass-sae-county-means.rda'))

