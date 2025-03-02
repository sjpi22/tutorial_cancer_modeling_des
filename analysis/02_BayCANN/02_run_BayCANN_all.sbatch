#!/usr/bin/bash
#SBATCH --job-name=run_baycann_all
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem=20G
#SBATCH --ntasks-per-node=4
#SBATCH --output=output/baycann_log_sherlock.log
#SBATCH --mail-user=sjpi@stanford.edu
#SBATCH --mail-type=ALL

#* To send the job to Sherlock use:
# - sbatch analysis/02_BayCANN/02_run_BayCANN_all.sbatch

# load R module
ml R/4.3.2

# run R code
Rscript analysis/02_BayCANN/02.1_generate_BayCANN_sample.R

# load modules for BayCANN and activate Python virtual environment
ml python/3 gcc R/4.3.2
source $HOME/tensorflow/bin/activate

# run R code
Rscript analysis/02_BayCANN/02.2_run_BayCANN.R
Rscript analysis/02_BayCANN/02.3_get_BayCANN_outputs.R
