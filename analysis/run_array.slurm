#!/bin/bash
#SBATCH --job-name=toy_brms
#SBATCH --array=1-6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/lustre/users/psych/dmm56/IllusionGameReliability/modeltoy_%A_%a.out
#SBATCH --error=/mnt/lustre/users/psych/dmm56/IllusionGameReliability/modeltoy_%A_%a.err

# Set working directory
cd /mnt/lustre/users/psych/dmm56/IllusionGameReliability || exit 1

echo "Running on $(hostname)"
echo "Current directory: $(/bin/pwd)"
echo "Files present: $(ls -lh)"
module load R/4.4.1-gfbf-2023b
Rscript models.R $SLURM_ARRAY_TASK_ID

echo "COMPLETED."