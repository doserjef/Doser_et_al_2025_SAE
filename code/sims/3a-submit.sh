#!/bin/bash
#
# Script:  3a-submit.sh
# Usage: For submitting multiple batch jobs to the NCSU HPC.
# Author: Jeffrey W. Doser (adapted from NCSU HPC)
#
## To run, type:
#     ./3a-submit.sh [# replicates] 
#  Script must have execute permissions, i.e.,
#     chmod u+x 3a-submit.sh

module load openmpi-gcc
module load R

if [ $# -ne 1 ]; then
        echo "Usage: You need to feed one argument to this program which is"
        echo "the number of replicates. For example,"
        echo "./3a-submit.sh 10"
        exit 1
fi

# Specify number of jobs to submit
numReps=$1

# Initialize loop counter
rep=1

while [ $rep -le $numReps ]
do

  echo "Submit job rep = $rep"

  bsub -n 1 -W 7200 -R span[hosts=1] -R "rusage[mem=10]" -q cnr -oo out.3a.rep.$rep -eo err.3a.rep.$rep "Rscript 3a-predict.R $rep"

  ((rep++))

done
