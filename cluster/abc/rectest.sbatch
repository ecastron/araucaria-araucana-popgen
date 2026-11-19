#!/bin/bash
#SBATCH --job-name=rectest
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/rectest_%j.out
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate /home/ecastron/envs/feems
cd /home/ecastron/araucaria
python scripts/abc/rectest.py
