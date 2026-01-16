#!/usr/bin/bash
#SBATCH --job-name=get_baycann_outputs
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=8
#SBATCH --output=output/logs/02.3_get_BayCANN_outputs_colorectal.log

#* To send the job to Sherlock use:
# - sbatch analysis/02_BayCANN/slurm/02.3_get_BayCANN_outputs.sbatch

# load R module
ml R/4.4.2

# run R code
Rscript analysis/02_BayCANN/02.3_get_BayCANN_outputs.R colorectal
