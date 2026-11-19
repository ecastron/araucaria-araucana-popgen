#!/bin/bash
#SBATCH --job-name=fig9
#SBATCH --partition=basic
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/fig9_%j.out
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate /home/ecastron/envs/feems
cd /home/ecastron/araucaria
python -u scripts/abc/f4_fig9.py
