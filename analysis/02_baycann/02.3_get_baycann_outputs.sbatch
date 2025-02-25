#!/usr/bin/bash
#SBATCH --job-name=get_baycann_outputs
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --ntasks-per-node=8
#SBATCH --output=output/baycann_outputs_sherlock.log
#SBATCH --mail-user=sjpi@stanford.edu
#SBATCH --mail-type=ALL

#* To send the job to Sherlock use:
# - sbatch analysis/02_baycann/02.3_get_baycann_outputs.sbatch

# load R module
ml R/4.3.2

# run R code
Rscript analysis/02_baycann/02.3_get_baycann_outputs.R
