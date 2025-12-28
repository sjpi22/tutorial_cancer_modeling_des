#!/usr/bin/bash
#SBATCH --job-name=run_IMABC
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=80G
#SBATCH --ntasks-per-node=24
#SBATCH --output=output/logs/02.2_run_IMABC_simulated.log

#* To send the job to Sherlock use:
# - sbatch analysis/02_IMABC/slurm/02.2_run_IMABC.sbatch

# load R module
ml R/4.4.2

# run R code
Rscript analysis/02_IMABC/02.2_run_IMABC.R simulated