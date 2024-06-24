#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH --ntasks-per-node=8
#SBATCH --output=generate_BayCANN_sample.%j.log
#SBATCH --error=generate_BayCANN_sample.%j.err

# load R module
ml R/4.3.2

# run R code
Rscript analysis/generate_BayCANN_sample.R