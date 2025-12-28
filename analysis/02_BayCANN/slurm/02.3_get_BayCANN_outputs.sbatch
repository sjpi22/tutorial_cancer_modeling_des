#!/usr/bin/bash
#SBATCH --job-name=get_baycann_outputs
#SBATCH --nodes=1
#SBATCH --time=16:00:00
#SBATCH --mem=30G
#SBATCH --ntasks-per-node=12
#SBATCH --output=output/logs/02.3_get_BayCANN_outputs_simulated.log
#SBATCH --mail-user=sjpi@stanford.edu
#SBATCH --mail-type=ALL

#* To send the job to Sherlock use:
# - sbatch analysis/02_BayCANN/slurm/02.3_get_BayCANN_outputs.sbatch

#* To get email updates about the job, add the following lines to the header and replace ___ with your email address:
#SBATCH --mail-user=___
#SBATCH --mail-type=ALL

# load R module
ml R/4.4.2

# run R code
Rscript analysis/02_BayCANN/02.3_get_BayCANN_outputs.R simulated
