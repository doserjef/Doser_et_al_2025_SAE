#!/bin/bash
#
# Script:  6b-submit.sh
# Usage: For submitting multiple batch jobs to the 
#        NCSU HPC.
# Author: Jeffrey W. Doser (adapted from NCSU HPC)
#
## To run, type:
#     ./6b-submit.sh [# folds] 
#  Script must have execute permissions, i.e.,
#     chmod u+x 3a-submit.sh

module load openmpi-gcc
module load R

if [ $# -ne 1 ]; then
        echo "Usage: You need to feed one argument to this program which is"
        echo "the number of folds in the cross-validation. For example,"
        echo "./6b-submit.sh 4"
        exit 1
fi

# Specify number of jobs to submit
numSpecies=$1

# Initialize loop counter
species=1

while [ $species -le $numSpecies ]
do

  echo "Submit job fold = $species"

  bsub -n 1 -W 7200 -R span[hosts=1] -R "rusage[mem=130]" -q cnr -oo out.6b.fold.$species -eo err.6b.fold.$species "Rscript 6b-main-val-stage-1.R $species"

  ((species++))

done

