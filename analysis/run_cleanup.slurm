#!/bin/bash
#SBATCH --job-name=cleanup
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/lustre/users/psych/dmm56/IllusionGameReliability/cleanup_%j.out
#SBATCH --error=/mnt/lustre/users/psych/dmm56/IllusionGameReliability/cleanup_%j.err

# Set working directory
cd /mnt/lustre/users/psych/dmm56/IllusionGameReliability || exit 1

# Load R module
module load R/4.4.1-gfbf-2023b

# Combine and save 
Rscript -e "library(brms); library(cmdstanr); library(rstan); m <- brms::combine_models(mlist = lapply(list.files(pattern = 'modeltoy.*rds$'), readRDS)); saveRDS(m, 'modeltoy.rds')"

# Clean up temporary files, output, and error logs (all modeltoy_* files)
rm -f /mnt/lustre/users/psych/dmm56/IllusionGameReliability/modeltoy_*.out
rm -f /mnt/lustre/users/psych/dmm56/IllusionGameReliability/modeltoy_*.err
# rm -f /mnt/lustre/users/psych/dmm56/IllusionGameReliability/modeltoy_task_*.rds

echo "COMPLETED."