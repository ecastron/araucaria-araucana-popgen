#!/bin/bash
#SBATCH --job-name=feas8
#SBATCH --partition=gpu,basic
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/feas8_%j.out
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate /home/ecastron/envs/feems
cd /home/ecastron/araucaria
python -u scripts/abc/feas8.py
