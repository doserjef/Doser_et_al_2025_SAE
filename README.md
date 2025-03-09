# Multivariate spatial models for small area estimation of species-specific forest inventory parameters 

### [Jeffrey W. Doser](https://www.doserlab.com/), [Malcolm S. Itter](https://www.mitter-forestecology.com/), [Grant M. Domke](https://scholar.google.com/citations?user=6ke2YZ4AAAAJ&hl=en), [Andrew O. Finley](https://www.finley-lab.com/) 

### In review 

### Please contact the first author for questions: Jeffrey W. Doser (jwdoser@ncsu.edu)

---------------------------------

## Abstract

National Forest Inventories (NFIs) provide statistically reliable information on forest resources at national and other large spatial scales. As forest management and conservation needs become increasingly complex, NFIs are being called upon to provide forest parameter estimates at spatial scales smaller than current design-based estimation procedures can provide. This is particularly true when estimates are desired by species or species groups, which is often required to inform wildlife habitat management, sustainable forestry certifications, or timber product assessments. Here we propose a multivariate spatial model for small area estimation of species-specific forest inventory parameters. The hierarchical Bayesian modeling framework accounts for key complexities in species-specific forest inventory data, such as zero-inflation, correlations among species, and residual spatial autocorrelation. Importantly, by fitting the model directly to the individual plot-level data, the framework enables estimates of species-level forest parameters, with associated uncertainty, across any user-defined small area of interest. A simulation study revealed minimal bias and higher accuracy of the proposed model-based approach compared to the design-based estimator and a non-parametric k-nearest neighbor (kNN) estimator. We applied the model to estimate species-specific county-level aboveground biomass for the 20 most abundant tree species in the southern United States using Forest Inventory and Analysis (FIA) data. Biomass estimates from the proposed model had high correlations with design-based estimates and kNN estimates. Importantly, the proposed model provided large gains in precision across all 20 species. On average across species, 91.5\% of county-level biomass estimates had higher precision compared to the design-based estimates. The proposed framework improves the ability of NFI data users to generate species-level forest parameter estimates with reasonable precision at management-relevant spatial scales.  

## Repository Directory

Certain scripts have special instructions for running, or may not be able to run without running previous scripts. These points in the scripts are indicated with NOTEs.  

### [code/fia](./code/fia/)

Contains all code to format and extract data, fit models, and summarize results from the analysis of species-specific biomass across the Southern US. Note that scripts should be run in the order indicated by the numbers in the file names. Fitting the models takes a large amount of RAM (~150GM), so the model-fitting scripts and prediction scripts would likely need to be run on a large server or in a high performance computing environment.  

+ `1a-get-FIA-data.R`: extracts and FIA data from `rFIA` for use in the analysis. 
+ `1b-get-covariates.R`: extracts data from Terraclimate, NLCD, and Terrain Tiles for use as auxiliary variables in the analysis.
+ `1c-format-data.R`: formats the FIA data into the necessary format for fitting the multivariate log-normal hurdle model in `spOccupancy` and `spAbundance`.  
+ `1d-case-study-data-plots.R`: produces some exploratory data anlysis plots as well as code to generate Figure 1 in the manuscript. 
+ `2a-main-stage-1.R`: fits the presence-absence portion of the multivariate log-normal hurdle model. 
+ `2a-submit.sh`: LSF batch script for running script `2a-main-stage-1.R` on the NC State HPC. 
+ `2b-main-stage-2.R`: fits Stage 2 of the multivariate log-normal hurdle model. 
+ `2b-submit.sh`: LSF batch script for running script `2b-main-stage-2.R` on the NC State HPC.
+ `3-extract-inits.R`: extract initial values from a previous run of the model to aid in convergence of subsequent model runs. 
+ `4a-get-pred-data.R`: extract coordinates and covariates across a grid of the Southern US for prediction of biomass across the region.
+ `4b-predict.R`: predicts county-level biomass for each of the 20 species across the Southern US. 
+ `4b-submit.sh`: LSF batch script for running script `4b-predict.R` on the NC State HPC. 
+ `5a-get-direct-estimates.R`: extracts the design-based estimates for county-level and species-specific biomass estimates across the South. 
+ `5b-get-kNN-estimates.R`: extracts the non-parametric k nearest neighbors estimates of species-specific county-level biomass across the Southern US.  
+ `6a-summary.R`: summarize results from the southern US case study and generate all figures included in the manuscript.  

### [code/sims](./code/sims/)

Contains all code to implement the simulation study. 

+ `1a-get-pop-cov-data.R`: extracts the plot locations and values of the covariates at those locations of the entire simulated population. This is the population from which we need to simulate each of the individual replicate data sets.  
+ `1b-get-pop-species-data.R`: simulates the actual species-level biomass values from the entire simulated population. 
+ `1c-get-replicate-data-sets.R`: extracts the individual data sets from the overall population for running multiple simulations with different replicates. 
+ `2a-main-stage-1.R`: runs Stage 1 of the simulation individually for each of the 100 replicate data sets. 
+ `2a-submit.sh`: bash script for submitting Stage 1 models for multiple replicate data sets simultaneously on the NC State HPC. 
+ `2b-main-stage-2.R`: runs Stage 2 of the simulation individually for each of the 100 replicate data sets. 
+ `2b-submit.sh`: bash script for submitting Stage 2 models for multiple replicate data sets simultaneously on the NC State HPC. 
+ `3a-predict.R`: predicts county-level biomass for each of the simulated data sets. 
+ `3a-submit.sh`: bash script for submitting the prediction script for each of the 100 replicate data sets. 
+ `3b-get-pop-means.R`: extracts the true population county-level and species-specific means for the simulated population.
+ `3c-get-get-direct-estimates.R`: calculates the direct estimates for each of the 100 simulated data sets. 
+ `3d-get-kNN-estimates.R`: calculates the kNN-based estimates for each of the 100 simulated data sets for each species and each county.
+ `4-summary.R`: summarizes results from the simulation study and produces two figures that highlight the bias and uncertainty of the proposed estimates. 

### [data](./data/)

Contains all data used in the Southern US FIA analysis as well as the simulation study. Note that the directory does not contain raw FIA data files.  

+ `covariate-data.rda`: covariate data at US FIA plot locations across the Southern US. 
+ `se_bio_stage_1_data.rda`: complete data object needed for fitting the Stage 1 model using `spOccupancy`. 
+ `se_bio_stage_2_data.rda`: complete data object needed for fitting the Stage 2 model using `spAbundance`. 
+ `se_fia_bio_data.rda`: temporary data object containing the US FIA biomass data for the Southern US. 
+ `se-prediction-data.rda`: coordinates and covariate values for the 1x1 km prediction grid across the Southern US used for generating county-level species-specific biomass estimates. 
+ `sim_data_replicates/`: directory containing the completely formatted data for the 100 replicate data sets used in the simulation study. Each individual data object consists of the `spOccupancy`-formatted data for Stage 1 and the `spAbundance`-formatted data for Stage 2. 
+ `sim_pop_county_true.rda`: the true county-level species-specific biomass values for the simulated population. 
+ `sim_pop_covariate_data.rda`: temporary object containing the covariate data for the entire simulated population. 
+ `sim_pop_data.rda`: complete data for the simulated population across the state of North Carolina.

### [results](./results/)

Note that model-based results files are too large to store on GitHub, but these files can be generated with the scripts in the `code/` directory.  

+ `inits-stage-1.rda`: initial values for Stage 1 obtained from a smaller model run used to aid in convergence in the final model run. 
+ `inits-stage-2.rda`: initial values for Stage 2 obtained from a smaller model run used to aid in convergence in the final model run. 
+ `k-means-sae-estimates.rda`: kNN-based estimates for the FIA analysis.
+ `se_direct_estimates.rda`: direct estimates for the FIA analysis. 
+ `sim_direct_ests.rda`: direct estimates for the simulation study. 
+ `sim_knn_ests.rda`: kNN-based estimates for the simulation study. 
+ `sim_pop_county_true.rda`: true population values for the simulation study. 
+ `sim_results/`: contains species-specific biomass small area estimates for each of the 100 replicate data sets analyzed in the simulation study. 

### [figures](./figures/)

Contains all figures included in the manuscript and supplemental information. The `species-maps/` sub-directory contains all the species-specific maps included in Supplemental Information S2. Additionally, for visual comparison with the direct estimates, the sub-directory `species-maps/direct/` contains county-level direct estimates for each species (where gray locations correspond to counties without any FIA plots). These figures can all be reproduced from the scripts included in this repository. 
