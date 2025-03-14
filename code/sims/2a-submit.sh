#!/bin/bash
#
# Script:  2a-submit.sh
# Usage: For submitting multiple batch jobs to the NCSU HPC.
# Author: Jeffrey W. Doser (adapted from NCSU HPC)
#
# To run, type:
#     ./2a-submit.sh [# replicates] 
#  Script must have execute permissions, i.e.,
#     chmod u+x 2a-submit.sh

module load openmpi-gcc
module load R

if [ $# -ne 1 ]; then
        echo "Usage: You need to feed one argument to this program which is"
        echo "the number of replicates. For example,"
        echo "./2a-submit.sh 10"
        exit 1
fi

# Specify number of jobs to submit
numReps=$1

# Initialize replicate loop counter
rep=1

while [ $rep -le $numReps ]
do

  echo "Submit job rep = $rep"

  bsub -n 1 -W 7200 -R span[hosts=1] -R "rusage[mem=10]" -q cnr -oo out.2a.rep.$rep -eo err.2a.rep.$rep "Rscript 2a-main-stage-1.R $rep"

  ((rep++))

done
