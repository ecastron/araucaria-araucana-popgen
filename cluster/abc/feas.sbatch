#!/bin/bash
#SBATCH --job-name=abcfeas
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/abcfeas_%j.out
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate /home/ecastron/envs/feems
cd /home/ecastron/araucaria
python scripts/abc/feas.py
