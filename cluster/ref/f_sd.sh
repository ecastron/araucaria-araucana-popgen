#!/bin/bash
#SBATCH --job-name=fsd
#SBATCH --partition=basic
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=6:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/fsd_%j.out
set -euo pipefail
module load Python/3.12.2
P=/home/ecastron/araucaria
python3 $P/scripts/ref/F_by_sd.py \
  $P/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf \
  $P/metadata/popmap_177.tsv \
  $P/results/sv/loci_copies.tsv \
  $P/results/sv/loci_sd.tsv
