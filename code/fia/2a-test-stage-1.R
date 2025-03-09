# 5a-main-stage-1.R: script to fit the presence-absence portion of the
#                    log-normal hurdle model using the top-performing model
#                    based on the WAIC comparison.
rm(list = ls())
library(spOccupancy)

# Directories -------------------------------------------------------------
# Determine what machine you're on and change directories as needed.
machine.name <- Sys.info()['nodename']
if (machine.name == 'pop-os') {
  out.dir <- 'results/'
  data.dir <- 'data/'
} else if (machine.name == 'birch') {
  out.dir <- '/mnt/disk4/jeff/DDF24/results/'
  data.dir <- 'data/'
} else { # Running on NCSU HPC
  out.dir <- '/share/doserlab/jwdoser/DIDF24/results/'
  data.dir <- '/share/doserlab/jwdoser/DIDF24/data/'
}

# Load spOccupancy formatted data -----------------------------------------
load(paste0(data.dir, "se_bio_stage_1_data.rda"))

# Prep the model ----------------------------------------------------------
# Priors
# high.dist <- 2000
# low.dist <- 50
high.dist <- 1100
low.dist <- 1000
# Might need further restrictions on phi later on.
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                   independent.betas = FALSE,
                   phi.unif = list(a = 3 / high.dist, b = 3 / low.dist))

# MCMC specifications
# n.batch <- 8000
# n.burn <- 140000
# n.thin <- 30
# batch.length <- 25

# n.batch <- 4000
# n.burn <- 40000
# n.thin <- 30
# batch.length <- 25

n.batch <- 2000
batch.length <- 25
n.burn <- 30000
n.thin <- 10

# n.batch <- 1000
# batch.length <- 25
# n.burn <- 15000
# n.thin <- 5

# n.batch <- 400
# batch.length <- 25
# n.burn <- 5000
# n.thin <- 5

# n.batch <- 10
# batch.length <- 25
# n.burn <- 0
# n.thin <- 1

# Specify spatial parameters
n.neighbors <- 5
cov.model <- "exponential"
n.samples <- n.batch * batch.length
n.chains <- 1

# Tuning and initial values
# Using default initial values
tuning.list <- list(phi = 0.1)

# Number of spatial factors to use
n.factors <- 5

# Load initial values
load(paste0(out.dir, "inits-stage-1.rda"))
inits.list$phi <- 3 / mean(c(low.dist, high.dist))

# Run the model -----------------------------------------------------------
out <- sfJSDM(formula = ~ scale(tmin) + I(scale(tmin)^2) +
                          scale(tmax) + I(scale(tmax)^2) +
                          scale(ppt) + I(scale(ppt)^2),
              data = data.list.1, priors = prior.list,
              tuning = tuning.list, inits = inits.list,
              n.neighbors = n.neighbors, cov.model = cov.model, NNGP = TRUE,
              n.factors = n.factors, n.batch = n.batch, batch.length = batch.length,
              n.burn = n.burn, accept.rate = 0.43, n.thin = n.thin,
              n.chains = n.chains, n.report = 1, n.omp.threads = 1)
summary(out)

# Save the results to hard drive
save(out, file = paste0(out.dir, 'se-stage-1-', n.samples, '-samples-',
                        n.factors, '-factors-', Sys.Date(), '.rda'))

# Save a small subset -----------------------------------------------------
beta.samples <- out$beta.samples
theta.samples <- out$theta.samples
lambda.samples <- out$lambda.samples
beta.comm.samples <- out$beta.comm.samples
tau.sq.beta.samples <- out$tau.sq.beta.samples
psi.means <- apply(out$psi.samples, c(2, 3), mean)
psi.sds <- apply(out$psi.samples, c(2, 3), sd)
w.means <- apply(out$w.samples, c(2, 3), mean)
w.sds <- apply(out$w.samples, c(2, 3), sd)
save(beta.samples, theta.samples, lambda.samples, beta.comm.samples,
     tau.sq.beta.samples, psi.means, w.means, w.sds, psi.sds,
     file = paste0(out.dir, 'se-small-high-range-stage-1-', n.samples, '-samples-',
                   n.factors, '-factors-', Sys.Date(), '.rda'))
