#!/bin/bash
#SBATCH --job-name=locimap
#SBATCH --partition=basic
#SBATCH --cpus-per-task=12
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/locimap_%j.out
set -euo pipefail
module load Minimap2/minimap2
D=/home/ecastron/araucaria/refs/araucana
S=/home/ecastron/araucaria/results/sv
cd "$S"
echo "node $(hostname) start $(date)"
# -N 50 keeps up to 50 alignments per locus; -p 0.6 keeps secondaries scoring
# >=60% of the primary, i.e. genuine duplicate copies rather than weak partial hits.
minimap2 -x sr -t 12 -N 50 -p 0.6 --secondary=yes -c \
  "$D/v1.sr.mmi" rad_loci.fa > loci_vs_genome.paf 2> loci_vs_genome.log
echo "alignments: $(wc -l < loci_vs_genome.paf)"
# count DISTINCT genomic placements per locus, at >=80% query coverage and >=90% identity
awk -F'\t' '($4-$3)/$2>=0.8 && $10/$11>=0.90 {print $1"\t"$6"\t"int($8/1000)}' loci_vs_genome.paf \
  | sort -u | cut -f1 | sort | uniq -c | awk '{print $2"\t"$1}' | sort > loci_copies.tsv
echo "loci with >=1 confident placement: $(wc -l < loci_copies.tsv)"
echo
echo "=== distribution of genomic copies per RAD locus ==="
awk '{c[$2>=6?"6+":$2]++} END{for(k in c) printf "  %-4s copies : %7d loci\n",k,c[k]}' loci_copies.tsv | sort
echo "done $(date)"
