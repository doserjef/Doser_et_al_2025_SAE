# 6c-main-val-stage-2.R: script to fit the lognormal portion of the
#                        log-normal hurdle model for the four-fold 
#                        cross-validation assessment. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(spAbundance)
library(spOccupancy)

# Get info from command line ----------------------------------------------
# This code is to extract the current k-fold data set from the command line
# to easily run the script for different folds 
args <- commandArgs(trailingOnly = TRUE)
# Current species
indx <- as.numeric(args[1])
# If not running the script from the command line, set the fold number manually 
# NOTE: uncomment this line if running locally
# TODO:
# indx <- 1
if(length(args) == 0) base::stop('Need to tell spOccupancy the data set')
print(indx)

# Directories -------------------------------------------------------------
# Determine what machine you're on and change directories as needed.
machine.name <- Sys.info()['nodename']
if (machine.name == 'pop-os') {
  out.dir <- 'results/'
  data.dir <- 'data/'
} else { # Running on NCSU HPC
  out.dir <- '/share/doserlab/jwdoser/DIDF24/results/'
  data.dir <- '/share/doserlab/jwdoser/DIDF24/data/'
}

# Load spOccupancy formatted data -----------------------------------------
load(paste0(data.dir, 'cross_val_data/hold_', indx, '_stage_2_data_validation_fit.rda'))

# Prep the model ----------------------------------------------------------
# Priors
high.dist <- 2000
low.dist <- 50
# Might need further restrictions on phi later on.
prior.list <- list(beta.comm.normal = list(mean = 0, var = 1000),
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                   tau.sq.ig = list(a = 2, b = 1),
                   independent.betas = TRUE,
                   phi.unif = list(a = 3 / high.dist, b = 3 / low.dist))

# MCMC specifications
n.batch <- 4000
n.burn <- 40000
n.thin <- 30
batch.length <- 25

# Specify spatial parameters
n.neighbors <- 5
cov.model <- "exponential"
n.samples <- n.batch * batch.length
n.chains <- 1

# Tuning and initial values
# Load initial values and tuning values from previous model fit.
tuning.list <- list(phi = 0.5)

load(paste0(out.dir, "inits-stage-2.rda"))
inits.list$w <- NULL

# Number of spatial factors
n.factors <- 5

# Run the model -----------------------------------------------------------
out <- sfMsAbund(formula = ~ scale(tcc) + scale(vpd) + scale(ppt) + 
                             scale(elev) + I(scale(elev)^2),
                 data = data.fit.2, priors = prior.list, inits = inits.list,
                 tuning = tuning.list, n.neighbors = n.neighbors, 
                 family = 'zi-Gaussian', cov.model = cov.model, NNGP = TRUE,
                 n.factors = n.factors, n.batch = n.batch,
                 batch.length = batch.length, n.burn = n.burn,
                 accept.rate = 0.43, n.thin = n.thin, n.chains = n.chains,
                 n.report = 1, n.omp.threads = 1)

# Save the results to a hard drive ----------------------------------------
save(out, file = paste0(out.dir, 'se-validation-stage-2-', n.samples, '-samples-',
                        n.factors, '-factors-', indx, '-fold-', Sys.Date(), '.rda'))

# Save a small subset -----------------------------------------------------
beta.samples <- out$beta.samples
theta.samples <- out$theta.samples
lambda.samples <- out$lambda.samples
beta.comm.samples <- out$beta.comm.samples
tau.sq.samples <- out$tau.sq.samples
tau.sq.beta.samples <- out$tau.sq.beta.samples
mu.means <- apply(out$mu.samples, c(2, 3), mean)
mu.sds <- apply(out$mu.samples, c(2, 3), sd)
w.means <- apply(out$w.samples, c(2, 3), mean)
w.sds <- apply(out$w.samples, c(2, 3), sd)
y.rep.quants <- apply(out$y.rep.samples, c(2, 3), quantile, c(0.025, 0.5, 0.975))
save(beta.samples, theta.samples, lambda.samples, beta.comm.samples,
     tau.sq.beta.samples, mu.means, w.means, w.sds, mu.sds, y.rep.quants, tau.sq.samples,
     file = paste0(out.dir, 'se-small-validation-stage-2-', n.samples, '-samples-',
                   n.factors, '-factors-', indx, '-fold-', Sys.Date(), '.rda'))

