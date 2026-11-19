#!/bin/bash
#SBATCH --job-name=arau-fam2
#SBATCH --partition=basic
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/family_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/family_%j.err
set -euo pipefail
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate dada2
Rscript /home/ecastron/araucaria/scripts/q6_family2.R
