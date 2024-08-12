#!/usr/bin/bash
#SBATCH --job-name=generate_BayCANN_sample
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --mem=20G
#SBATCH --ntasks-per-node=8
#SBATCH --output=output/calibration/BayCANN/generate_BayCANN_sample.log
#SBATCH --mail-user=sjpi@stanford.edu
#SBATCH --mail-type=ALL

#* To send the job to Sherlock use:
# - sbatch analysis/calibration/BayCANN/generate_BayCANN_sample.sbatch

# load R module
ml R/4.3.2

# run R code
Rscript analysis/calibration/BayCANN/generate_BayCANN_sample.R