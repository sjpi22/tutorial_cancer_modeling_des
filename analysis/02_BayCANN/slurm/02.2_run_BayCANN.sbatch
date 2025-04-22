#!/usr/bin/bash
#SBATCH --job-name=run_baycann
#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --mem=4G
#SBATCH --ntasks-per-node=1
#SBATCH --output=output/baycann_run.log
#SBATCH --mail-user=sjpi@stanford.edu
#SBATCH --mail-type=ALL

#* To send the job to Sherlock use:
# - sbatch analysis/02_BayCANN/slurm/02.2_run_BayCANN.sbatch

# load modules for BayCANN and activate Python virtual environment
ml python/3 gcc R/4.3.2
source $HOME/tensorflow/bin/activate

# run R code
Rscript analysis/02_BayCANN/02.2_run_BayCANN.R
