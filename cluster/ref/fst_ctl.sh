#!/bin/bash
#SBATCH --job-name=fstctl
#SBATCH --partition=basic
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=8:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/fstctl_%j.out
set -euo pipefail
module load Python/3.12.2
P=/home/ecastron/araucaria
python3 $P/scripts/ref/fst_estimator_control.py \
  $P/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf \
  $P/metadata/popmap_177.tsv \
  $P/results/sv/fst_estimator_control.tsv 6
