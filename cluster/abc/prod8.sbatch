#!/bin/bash
#SBATCH --job-name=abc8
#SBATCH --partition=gpu,basic
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --array=1-40%20
#SBATCH --output=/home/ecastron/araucaria/logs/abc8_%A_%a.out
# Clock ABC (sim8): literature priors on mu per year and generation time; 35 sims per model per task.
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate /home/ecastron/envs/feems
cd /home/ecastron/araucaria
mkdir -p results/Q1/abc_clock
python -u scripts/abc/run_sims8.py ${SLURM_ARRAY_TASK_ID} 35 results/Q1/abc_clock/s8_${SLURM_ARRAY_TASK_ID}.npz H1,H1b,H2,H3
