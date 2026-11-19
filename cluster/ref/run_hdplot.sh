#!/bin/bash
#SBATCH --job-name=hdplot
#SBATCH --partition=basic
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/hdplot_%j.out
set -euo pipefail
python3 /home/ecastron/araucaria/scripts/ref/hdplot.py \
  --vcf /home/ecastron/araucaria/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf \
  --minind 30 --out /home/ecastron/araucaria/results/refmap/hdplot.tsv
