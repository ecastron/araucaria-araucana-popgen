#!/bin/bash
#SBATCH --job-name=migrlong
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=96:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/migrate_long_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/migrate_long_%j.err
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate /home/ecastron/envs/migrate
cd /home/ecastron/araucaria/results/Q1/migrate
echo "start $(date)"
echo "Theta prior widened to Uniform(0, 0.5); 400,000 steps x 3 replicates; Gelman-Rubin on."
migrate-n parmfile.a200long
echo "end $(date)"
