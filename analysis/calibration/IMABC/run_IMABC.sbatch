#!/usr/bin/bash
#SBATCH --job-name=run_IMABC
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=20G
#SBATCH --ntasks-per-node=24
#SBATCH --output=output/calibration/IMABC/run_IMABC.log
#SBATCH --mail-user=sjpi@stanford.edu
#SBATCH --mail-type=ALL  

# load R module
ml R/4.3.2

# run R code
Rscript analysis/calibration/IMABC/run_IMABC.R