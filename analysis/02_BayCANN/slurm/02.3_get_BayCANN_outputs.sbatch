#!/usr/bin/bash
#SBATCH --job-name=get_baycann_outputs
#SBATCH --nodes=1
#SBATCH --time=16:00:00
#SBATCH --mem=30G
#SBATCH --ntasks-per-node=12
#SBATCH --output=output/baycann_outputs_sherlock.log

#* To send the job to Sherlock use:
# - sbatch analysis/02_BayCANN/slurm/02.3_get_BayCANN_outputs.sbatch

# load R module
ml R/4.3.2

# run R code
Rscript analysis/02_BayCANN/02.3_get_BayCANN_outputs.R
