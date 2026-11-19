#!/bin/bash
#SBATCH --job-name=extloci
#SBATCH --partition=basic
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/extloci_%j.out
set -euo pipefail
python3 /home/ecastron/araucaria/scripts/ref/loci_to_fasta.py \
  /home/ecastron/araucaria/rad/ipyrad_full/full_ct091_outfiles/full_ct091.loci \
  /home/ecastron/araucaria/results/sv/rad_loci.fa
grep -c "^>" /home/ecastron/araucaria/results/sv/rad_loci.fa
