#!/bin/bash
#BSUB -n 1
#BSUB -W 7200
#BSUB -J se-1
#BSUB -o stdout.%J
#BSUB -e stderr.%J
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=170]"
#BSUB -q cnr
module load openmpi-gcc
module load R
Rscript 2a-main-stage-1.R
