# 3a-extract-inits.R: extract initial values for a subsequent run of the model 
#                     to aid in convergence. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)

# NOTE: this script will not run using the files available on GitHub, because
# the results files are too large to include on GitHub. Instead, if you run 
# the previous main scripts, you can read in the output file into this file
# to extract the initial values (note that you'll need to change the file 
# name). 

# Load the data -----------------------------------------------------------
load('data/se_bio_stage_1_data.rda')
load('data/se_bio_stage_2_data.rda')

# Stage 1 -----------------------------------------------------------------
# Load model results
load('results/se-small-stage-1-1e+05-samples-5-factors-2025-02-13.rda')

N <- nrow(data.list.1$y)
theta.inits <- apply(theta.samples, 2, median)
beta.comm.inits <- apply(beta.comm.samples, 2, median)
lambda.inits <- matrix(apply(lambda.samples, 2, median), nrow = N)
tau.sq.beta.inits <- apply(tau.sq.beta.samples, 2, median)
beta.inits <- matrix(apply(beta.samples, 2, median), nrow = N)
w.inits <- w.means

inits.list <- list(phi = theta.inits, beta.comm = beta.comm.inits, 
                   lambda = lambda.inits, tau.sq.beta = tau.sq.beta.inits,
                   beta = beta.inits, w = w.inits)
# Save initial values for Stage 1.
save(inits.list, file = 'results/inits-stage-1.rda')

# Stage 2 -----------------------------------------------------------------
load('results/se-small-stage-2-50000-samples-5-factors-2025-02-12.rda')
N <- nrow(data.list.2$y)
theta.inits <- apply(theta.samples, 2, median)
beta.comm.inits <- apply(beta.comm.samples, 2, median)
lambda.inits <- matrix(apply(lambda.samples, 2, median), nrow = N)
tau.sq.beta.inits <- apply(tau.sq.beta.samples, 2, median)
beta.inits <- matrix(apply(beta.samples, 2, median), nrow = N)
w.inits <- w.means

inits.list <- list(phi = theta.inits, beta.comm = beta.comm.inits, 
                   lambda = lambda.inits, tau.sq.beta = tau.sq.beta.inits,
                   beta = beta.inits, w = w.inits)
# Save initial values for Stage 2.                    
save(inits.list, file = 'results/inits-stage-2.rda')

