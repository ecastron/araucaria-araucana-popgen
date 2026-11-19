#!/bin/bash
#SBATCH --job-name=abcdiag4
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/abcdiag4_%j.out
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate /home/ecastron/envs/feems
cd /home/ecastron/araucaria
python -u scripts/abc/diag4.py
