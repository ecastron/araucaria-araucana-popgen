#!/bin/bash
#SBATCH --job-name=indivdist
#SBATCH --partition=basic
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/indivdist_%j.out
set -eo pipefail
python3 /home/ecastron/araucaria/scripts/ref/indiv_distance.py \
  /home/ecastron/araucaria/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf \
  /home/ecastron/araucaria/metadata/popmap_177.tsv \
  /home/ecastron/araucaria/results/sv/indiv_dist.tsv
head -2 /home/ecastron/araucaria/results/sv/indiv_dist.tsv | cut -c1-120
