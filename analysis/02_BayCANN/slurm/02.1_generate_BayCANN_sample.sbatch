#!/usr/bin/bash
#SBATCH --job-name=generate_baycann_sample
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --ntasks-per-node=8
#SBATCH --output=output/logs/02.1_generate_BayCANN_sample_simulated.log

#* To send the job to Sherlock use:
# - sbatch analysis/02_BayCANN/slurm/02.1_generate_BayCANN_sample.sbatch

# load R module
ml R/4.4.2

# run R code
Rscript analysis/02_BayCANN/02.1_generate_BayCANN_sample.R simulated
