#!/bin/bash
#SBATCH --job-name=arau-filter
#SBATCH --partition=basic
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/filtergrid_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/filtergrid_%j.err
set -uo pipefail
python3 /home/ecastron/araucaria/scripts/snp_filter_grid.py
