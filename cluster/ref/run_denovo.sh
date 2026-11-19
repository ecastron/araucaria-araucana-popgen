#!/bin/bash
#SBATCH --job-name=denovobal
#SBATCH --partition=basic
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/denovobal_%j.out
set -euo pipefail
python3 /home/ecastron/araucaria/scripts/ref/balance_denovo.py \
  /home/ecastron/araucaria/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf \
  --libs /home/ecastron/araucaria/scripts/ref/selected_libs.txt \
  > /home/ecastron/araucaria/results/refmap/denovo_balance_by_depth.txt
echo "done"; head -30 /home/ecastron/araucaria/results/refmap/denovo_balance_by_depth.txt
