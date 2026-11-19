#!/bin/bash
#SBATCH --job-name=dosecov
#SBATCH --partition=basic
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=12:00:00
#SBATCH --array=1-8
#SBATCH --output=/home/ecastron/araucaria/logs/dosecov_%A_%a.out
set -euo pipefail
module load SAMTOOLS/samtools.1.20
R=/home/ecastron/araucaria/results/sv/dosage
i=$SLURM_ARRAY_TASK_ID
# -Q 20 matches the mapping-quality filter used everywhere else in this project
samtools bedcov -Q 20 "$R/sc_loci.part$i.bed" \
  $(cat /home/ecastron/araucaria/results/sv/bamlist177.txt | tr '\n' ' ') \
  > "$R/bedcov.part$i.tsv"
echo "part $i rows: $(wc -l < "$R/bedcov.part$i.tsv")"
