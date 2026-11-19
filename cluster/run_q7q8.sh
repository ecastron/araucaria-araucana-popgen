#!/bin/bash
#SBATCH --job-name=arau-Q7Q8
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/Q7Q8_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/Q7Q8_%j.err
set -uo pipefail
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate dada2
Rscript /home/ecastron/araucaria/scripts/q7_q8.R
echo "Q7Q8 COMPLETE"
