#!/usr/bin/bash
#SBATCH --job-name=get_IMABC_outputs
#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --ntasks-per-node=8
#SBATCH --output=output/logs/02.4_get_IMABC_outputs_bladder.log
#SBATCH --mail-user=sjpi@stanford.edu
#SBATCH --mail-type=ALL

#* To send the job to Sherlock use:
# - sbatch analysis/02_IMABC/slurm/02.4_get_IMABC_outputs.sbatch

#* To get email updates about the job, add the following lines to the header and replace ___ with your email address:
#SBATCH --mail-user=___
#SBATCH --mail-type=ALL

# load R module
ml R/4.4.2

# run R code
Rscript analysis/02_IMABC/02.4_get_IMABC_outputs.R bladder
