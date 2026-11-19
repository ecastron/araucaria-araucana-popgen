#!/bin/bash
#SBATCH --job-name=abc9
#SBATCH --partition=gpu,basic
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --array=1-40%20
#SBATCH --output=/home/ecastron/araucaria/logs/abc9_%A_%a.out
# Clock ABC (sim9): mu log-uniform 4.1e-10 to 1.4e-9 per year (De La Torre et al. 2017) and generation time; 35 sims per model per task.
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate /home/ecastron/envs/feems
cd /home/ecastron/araucaria
mkdir -p results/Q1/abc_clock_lit
python -u scripts/abc/run_sims9.py ${SLURM_ARRAY_TASK_ID} 35 results/Q1/abc_clock_lit/s9_${SLURM_ARRAY_TASK_ID}.npz H1,H1b,H2,H3
