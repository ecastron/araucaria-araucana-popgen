#!/bin/bash
#SBATCH --job-name=dosetest
#SBATCH --partition=basic
#SBATCH --cpus-per-task=2
#SBATCH --mem=60G
#SBATCH --time=8:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/dosetest_%j.out
set -euo pipefail
module load Python/3.12.2
python3 /home/ecastron/araucaria/scripts/ref/dosage_test.py
