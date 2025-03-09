# 2a-main-stage-1.R: script to run the first stage of the simulation 
#                    individually for each of the data sets. 
# Author: Jeffrey W. Doser
rm(list = ls())
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
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                   independent.betas = FALSE,
                   phi.unif = list(a = 3 / high.dist, b = 3 / low.dist))

# MCMC specifications
n.batch <- 2000
batch.length <- 25
n.burn <- 10000
n.thin <- 40

# Specify spatial parameters
n.neighbors <- 5
cov.model <- "exponential"
n.samples <- n.batch * batch.length
n.chains <- 1

# Tuning and initial values
# Using default initial values
tuning.list <- list(phi = 0.1)

# Number of spatial factors to use
n.factors <- 4

# Run the model -----------------------------------------------------------
out <- sfJSDM(formula = ~ tmin + tmin.2 +
                          tmax + tmax.2 +
                          ppt + ppt.2,
              data = data.list.1, priors = prior.list,
              tuning = tuning.list,
              n.neighbors = n.neighbors, cov.model = cov.model, NNGP = TRUE,
              n.factors = n.factors, n.batch = n.batch, batch.length = batch.length,
              n.burn = n.burn, accept.rate = 0.43, n.thin = n.thin,
              n.chains = n.chains, n.report = 1, n.omp.threads = 1)
summary(out)

# Save the results to hard drive ------------------------------------------
save(out, file = paste0(out.dir, 'sim-stage-1-', n.samples, '-samples-',
                        sim.indx, '-replicate-', Sys.Date(), '.rda'))

