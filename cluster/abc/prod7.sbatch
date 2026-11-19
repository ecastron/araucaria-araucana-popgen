#!/bin/bash
#SBATCH --job-name=abc7
#SBATCH --partition=gpu,basic
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=20:00:00
#SBATCH --array=1-40%20
#SBATCH --output=/home/ecastron/araucaria/logs/abc7_%A_%a.out
# Identical to prod5.sbatch except: sim7 (T1 prior floor 10 generations) and a new seed base.
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate /home/ecastron/envs/feems
cd /home/ecastron/araucaria
python -u scripts/abc/run_sims7.py ${SLURM_ARRAY_TASK_ID} 35    /home/ecastron/araucaria/results/Q1/abc_floor10/s7_${SLURM_ARRAY_TASK_ID}.npz H1,H1b,H2,H3
