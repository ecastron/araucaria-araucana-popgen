#!/bin/bash
#SBATCH --job-name=dosectl
#SBATCH --partition=basic
#SBATCH --cpus-per-task=2
#SBATCH --mem=60G
#SBATCH --time=4:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/dosectl_%j.out
set -euo pipefail
module load Python/3.12.2
python3 /home/ecastron/araucaria/scripts/ref/dosage_control.py
