#!/bin/bash
#SBATCH --job-name=arau-its-r1
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/its_r1_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/its_r1_%j.err
set -euo pipefail
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate dada2
Rscript /home/ecastron/araucaria/scripts/dada2_its_r1.R
