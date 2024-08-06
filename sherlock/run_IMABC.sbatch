#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=20G
#SBATCH --ntasks-per-node=24
#SBATCH --output=run_IMABC.%j.log
#SBATCH --error=run_IMABC.%j.err

# load R module
ml R/4.3.2

# run R code
Rscript analysis/run_IMABC.R