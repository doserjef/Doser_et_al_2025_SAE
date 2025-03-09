#!/bin/bash
#BSUB -n 1
#BSUB -W 7200
#BSUB -J se-2
#BSUB -o stdout.%J
#BSUB -e stderr.%J
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=150]"
#BSUB -q cnr
module load openmpi-gcc
module load R
Rscript 2b-main-stage-2.R
