#!/usr/bin/bash
#SBATCH --job-name=get_true_outcomes
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mem=3G
#SBATCH --ntasks-per-node=4
#SBATCH --output=output/true_outputs.log
#SBATCH --mail-user=sjpi@stanford.edu
#SBATCH --mail-type=ALL

#* To send the job to Sherlock use:
# - sbatch _ground_truth/slurm/02_get_true_outcomes.sbatch

# load R module
ml R/4.3.2

# run R code
Rscript _ground_truth/02_get_true_outcomes.R
