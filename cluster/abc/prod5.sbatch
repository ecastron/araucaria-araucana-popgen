#!/bin/bash
#SBATCH --job-name=abc5
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=20:00:00
#SBATCH --array=1-40%20
#SBATCH --output=/home/ecastron/araucaria/logs/abc5_%A_%a.out
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate /home/ecastron/envs/feems
cd /home/ecastron/araucaria
python -u scripts/abc/run_sims5.py ${SLURM_ARRAY_TASK_ID} 35 \
   /home/ecastron/araucaria/results/Q1/abc/s5_${SLURM_ARRAY_TASK_ID}.npz H1,H1b,H2,H3
