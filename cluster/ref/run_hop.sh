#!/bin/bash
#SBATCH --job-name=hoptest
#SBATCH --partition=basic
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/hoptest_%j.out
set -euo pipefail
python3 /home/ecastron/araucaria/scripts/ref/lane_hopping_test.py \
  --vcf /home/ecastron/araucaria/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf \
  --lanes /home/ecastron/araucaria/results/refmap/lane_assignment.tsv \
  --popmap /home/ecastron/araucaria/metadata/popmap_177.tsv \
  --mindp 20 --maxsites 40000
