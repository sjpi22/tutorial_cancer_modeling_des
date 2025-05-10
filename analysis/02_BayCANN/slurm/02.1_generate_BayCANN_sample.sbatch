#!/usr/bin/bash
#SBATCH --job-name=generate_baycann_sample
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mem=3G
#SBATCH --ntasks-per-node=4
#SBATCH --output=output/baycann_sample.log

#* To send the job to Sherlock use:
# - sbatch analysis/02_BayCANN/slurm/02.1_generate_BayCANN_sample.sbatch

# load R module
ml R/4.3.2

# run R code
Rscript analysis/02_BayCANN/02.1_generate_BayCANN_sample.R
