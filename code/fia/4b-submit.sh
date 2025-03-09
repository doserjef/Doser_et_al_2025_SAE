#!/bin/bash
#BSUB -n 5
#BSUB -W 7200
#BSUB -J se-1
#BSUB -o stdout.%J
#BSUB -e stderr.%J
# #BSUB -x 
#BSUB -R span[hosts=1]
#BSUB -R "rusage[mem=70]"
#BSUB -q cnr
module load openmpi-gcc
module load R
Rscript 4b-predict.R
