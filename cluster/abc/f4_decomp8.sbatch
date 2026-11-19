#!/bin/bash
#SBATCH --job-name=dec8
#SBATCH --partition=basic
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/dec8_%j.out
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate /home/ecastron/envs/feems
cd /home/ecastron/araucaria
python -u scripts/abc/f4_decomp8.py
