# 2b-main-stage-2.R: script to run the second stage models for all simulation
#                    replicates
# Author: Jeffrey W. Doser
rm(list = ls())
library(spAbundance)
library(spOccupancy)

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
# sim.indx <- 1
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
} else { # Running on NCSU HPC
  out.dir <- '/share/doserlab/jwdoser/DIDF24/results/sim_results/'
  data.dir <- '/share/doserlab/jwdoser/DIDF24/data/sim_data_replicates/'
}

# Load spOccupancy formatted data -----------------------------------------
load(paste0(data.dir, "sim_", sim.indx, '_data.rda'))

# Prep the model ----------------------------------------------------------
# Priors
high.dist <- 500
low.dist <- 50
# Might need further restrictions on phi later on.
prior.list <- list(beta.comm.normal = list(mean = 0, var = 1000),
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                   tau.sq.ig = list(a = 2, b = 1),
                   independent.betas = TRUE,
                   phi.unif = list(a = 3 / high.dist, b = 3 / low.dist))

# MCMC specifications
n.batch <- 2000
batch.length <- 25
n.burn <- 20000
n.thin <- 30 

# Specify spatial
n.neighbors <- 5
cov.model <- "exponential"
n.samples <- n.batch * batch.length
n.chains <- 1

# Tuning and initial values
# Load initial values and tuning values from previous model fit.
tuning.list <- list(phi = 0.5)

# Fix initial values of community-level parameters, which will serve as the
# prior values.
inits.list <- list(beta.comm = 0, tau.sq.beta = 10)

# Number of spatial factors
n.factors <- 4

# Run the model -----------------------------------------------------------
out <- sfMsAbund(formula = ~ tcc + vpd + ppt + elev + elev.2,
                 data = data.list.2, priors = prior.list, inits = inits.list,
                 tuning = tuning.list, n.neighbors = n.neighbors, 
                 family = 'zi-Gaussian', cov.model = cov.model, NNGP = TRUE,
                 n.factors = n.factors, n.batch = n.batch,
                 batch.length = batch.length, n.burn = n.burn,
                 accept.rate = 0.43, n.thin = n.thin, n.chains = n.chains,
                 n.report = 1, n.omp.threads = 1)

# Save the results to a hard drive ----------------------------------------
save(out, file = paste0(out.dir, 'sim-stage-2-', n.samples, '-samples-',
                        sim.indx, '-replicate-', Sys.Date(), '.rda'))
