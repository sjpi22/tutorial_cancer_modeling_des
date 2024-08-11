#!/usr/bin/bash
#SBATCH --job-name=generate_BayCANN_sample
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH --ntasks-per-node=8
#SBATCH --output=generate_BayCANN_sample.%j.log
#SBATCH --mail-user=sjpi@stanford.edu
#SBATCH --mail-type=ALL  

#* To send the job to Sherlock use:
# - sbatch --export=ALL sbatch/03_job_array.sh

# load R module
ml R/4.3.2

# run R code
Rscript analysis/generate_BayCANN_sample.R