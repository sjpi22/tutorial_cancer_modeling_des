#!/usr/bin/bash
#SBATCH --job-name=run_IMABC
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=16
#SBATCH --output=output/logs/02.2_run_IMABC.log
#SBATCH --mail-user=sjpi
#SBATCH --mail-type=ALL

#* To send the job to Sherlock use:
# - sbatch analysis/02_IMABC/slurm/02.2_run_IMABC.sbatch

#* To get email updates about the job, add the following lines to the header and replace ___ with your email address:
#SBATCH --mail-user=___
#SBATCH --mail-type=ALL

# load R module
ml R/4.4.2

# run R code
Rscript analysis/02_IMABC/02.2_run_IMABC.R