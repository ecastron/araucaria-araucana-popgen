#!/bin/bash
#SBATCH --job-name=q6soil
#SBATCH --partition=basic
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=6:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/q6lit_%j.out
set -euo pipefail
module load R
Rscript /home/ecastron/araucaria/scripts/q6_paired_cordillera_lit.R
