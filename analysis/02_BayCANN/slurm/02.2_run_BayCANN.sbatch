#!/usr/bin/bash
#SBATCH --job-name=run_baycann
#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --mem=20G
#SBATCH --ntasks-per-node=4
#SBATCH --output=output/logs/02.2_run_BayCANN.log
#SBATCH --mail-user=sjpi@stanford.edu
#SBATCH --mail-type=ALL

#* To send the job to Sherlock use:
# - sbatch analysis/02_BayCANN/slurm/02.2_run_BayCANN.sbatch

#* To get email updates about the job, add the following lines to the header and replace ___ with your email address:
#SBATCH --mail-user=___
#SBATCH --mail-type=ALL

# load modules for BayCANN and activate Python virtual environment
ml python/3 gcc R/4.4.2
source $HOME/tensorflow/bin/activate

# run R code
Rscript analysis/02_BayCANN/02.2_run_BayCANN.R
