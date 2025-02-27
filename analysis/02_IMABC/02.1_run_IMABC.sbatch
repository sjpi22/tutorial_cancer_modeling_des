#!/usr/bin/bash
#SBATCH --job-name=run_IMABC
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mem=20G
#SBATCH --ntasks-per-node=8
#SBATCH --output=output/imabc_sherlock.log
#SBATCH --mail-user=sjpi@stanford.edu
#SBATCH --mail-type=ALL

#* To send the job to Sherlock use:
# - sbatch analysis/02_IMABC/02.1_run_IMABC.sbatch

# load R module
ml R/4.3.2

# run R code
Rscript analysis/02_IMABC/02.1_run_IMABC.R