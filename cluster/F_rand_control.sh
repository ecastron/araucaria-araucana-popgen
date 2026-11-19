#!/bin/bash
#SBATCH -J Frand
#SBATCH -p gpu
#SBATCH --nodelist=g001
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 3:00:00
#SBATCH -o /home/ecastron/araucaria/logs/Frand_%j.out
#SBATCH -e /home/ecastron/araucaria/logs/Frand_%j.err
set -uo pipefail
P=/home/ecastron/araucaria
R=$P/results/Q1/F_by_depth
export PATH=/opt/miniconda/miniconda-py312_24.3.0-0/envs/popgen/bin:$PATH
V=$P/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf
DBG=$R/F_rand_debug.log; : > "$DBG"
# retention fractions matched to the depth arms: dp20 = 47.6%, dp30 = 24.6%
for spec in "rand476 0.476" "rand246 0.246"; do
  set -- $spec; NAME=$1; KEEP=$2
  B=$R/gl_${NAME}.beagle
  python3 $P/scripts/catg_to_beagle_rand.py "$V" "$B" 0.01 0 "$KEEP" 1 >> "$DBG" 2>&1 || continue
  gzip -f "$B"
  NS=$(zcat "$B.gz" | head -1 | awk '{print (NF-3)/3}')
  NL=$(( $(zcat "$B.gz" | wc -l) - 1 ))
  zcat "$B.gz" | head -1 | cut -f4- | tr '\t' '\n' | awk 'NR%3==1' > "$R/samples_${NAME}.txt"
  ngsRelate -G "$B.gz" -n "$NS" -L "$NL" -z "$R/samples_${NAME}.txt" \
    -O "$R/ngsrelate_${NAME}.tsv" >> "$DBG" 2>&1
done
