#!/bin/bash
#SBATCH --job-name=arau-Q5
#SBATCH --partition=basic   # ASV tables are small; leaves g001 free for Q1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=48:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/Q5_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/Q5_%j.err
set -uo pipefail
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate dada2
Rscript /home/ecastron/araucaria/scripts/q5_community.R
echo "Q5 COMPLETE"
