#!/bin/bash
#SBATCH --job-name=fbfst
#SBATCH --partition=basic
#SBATCH --cpus-per-task=2
#SBATCH --mem=30G
#SBATCH --time=8:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/fbfst_filt_%j.out
set -euo pipefail
module load Python/3.12.2
python3 /home/ecastron/araucaria/scripts/ref/fb_fst_pairs.py
