#!/bin/bash
#SBATCH --job-name=arau-dada2
#SBATCH --partition=basic
#SBATCH --array=1-2
#SBATCH --cpus-per-task=12
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/dada2_%A_%a.out
#SBATCH --error=/home/ecastron/araucaria/logs/dada2_%A_%a.err
set -uo pipefail
AMPS=(16S ITS)
AMP=${AMPS[$((SLURM_ARRAY_TASK_ID-1))]}
E=/opt/miniconda/miniconda-py312_24.3.0-0/envs/dada2
export PATH=$E/bin:$PATH
echo "=== DADA2 $AMP ==="; date
$E/bin/Rscript /home/ecastron/araucaria/scripts/dada2_step1.R "$AMP"
echo "exit=$?"; date
