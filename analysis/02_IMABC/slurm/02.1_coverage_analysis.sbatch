#!/usr/bin/bash
#SBATCH --job-name=coverage_analysis
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mem=24G
#SBATCH --ntasks-per-node=8
#SBATCH --output=output/logs/02.1_coverage_analysis.log
#SBATCH --mail-user=sjpi@stanford.edu
#SBATCH --mail-type=ALL

#* To send the job to Sherlock use:
# - sbatch analysis/02_IMABC/slurm/02.1_coverage_analysis.sbatch

#* To get email updates about the job, add the following lines to the header and replace ___ with your email address:
#SBATCH --mail-user=___
#SBATCH --mail-type=ALL

# load R module
ml R/4.4.2

# run R code
Rscript analysis/02_IMABC/02.1_coverage_analysis.R