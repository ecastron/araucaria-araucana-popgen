#!/bin/bash
#SBATCH -J f4n7
#SBATCH --partition=basic
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH -o /home/ecastron/araucaria/logs/f4n7_%j.out
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate /home/ecastron/envs/feems
python -u /home/ecastron/araucaria/scripts/abc/f4_numbers7.py
