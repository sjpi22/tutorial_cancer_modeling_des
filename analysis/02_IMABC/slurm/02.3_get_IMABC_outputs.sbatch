#!/usr/bin/bash
#SBATCH --job-name=get_IMABC_outputs
#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --mem=10G
#SBATCH --ntasks-per-node=8
#SBATCH --output=output/imabc_outputs_sherlock.log

#* To send the job to Sherlock use:
# - sbatch analysis/02_IMABC/slurm/02.3_get_IMABC_outputs.sbatch

# load R module
ml R/4.3.2

# run R code
Rscript analysis/02_IMABC/02.3_get_IMABC_outputs.R
