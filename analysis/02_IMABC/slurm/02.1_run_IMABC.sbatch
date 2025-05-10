#!/usr/bin/bash
#SBATCH --job-name=run_IMABC
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --mem=8G
#SBATCH --ntasks-per-node=12
#SBATCH --output=output/imabc_sherlock.log

#* To send the job to Sherlock use:
# - sbatch analysis/02_IMABC/slurm/02.1_run_IMABC.sbatch

# load R module
ml R/4.3.2

# run R code
Rscript analysis/02_IMABC/02.1_run_IMABC.R