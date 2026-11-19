#!/bin/bash
#SBATCH --job-name=nupt
#SBATCH --partition=basic
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/nupt_%j.out
set -euo pipefail
module load Minimap2/minimap2
D=/home/ecastron/araucaria/refs/araucana
O=/home/ecastron/araucaria/results/refmap
cd "$D"
# Align the 146 kb chloroplast against the 19.3 Gb nuclear assembly.
# Two very different things will hit:
#   (a) genuine plastid contigs that assembled alongside the nuclear genome
#       -- near-full-length, ~100% identity
#   (b) NUPTs, nuclear insertions of plastid DNA -- partial, degraded identity
# Only (b) creates mixed templates at a nuclear RAD locus, so they must be split.
minimap2 -t 12 -x asm10 -c --secondary=yes -N 50 v1.fasta v1_chloroplast.fasta \
  > "$O/cp_vs_nuclear.paf" 2> "$O/cp_vs_nuclear.log"
echo "total alignments: $(wc -l < "$O/cp_vs_nuclear.paf")"
awk -F'\t' '{
  ident=$10/$11; alen=$11; tgt=$6; tlen=$7;
  if (ident>=0.99 && alen>=100000) cat="plastid_contig";
  else if (alen>=1000) cat="NUPT_long";
  else cat="NUPT_short";
  n[cat]++; bp[cat]+=alen; if(!(tgt in seen[cat])){seen[cat][tgt]=1; c[cat]++}
} END {
  printf "%-16s %8s %14s %10s\n","class","n_aln","aligned_bp","n_contigs";
  for(k in n) printf "%-16s %8d %14d %10d\n",k,n[k],bp[k],c[k];
}' "$O/cp_vs_nuclear.paf"
# BED of candidate NUPT regions (excluding whole plastid contigs)
awk -F'\t' '($10/$11)<0.99 || $11<100000 {if($11>=300) print $6"\t"$8"\t"$9"\t"sprintf("%.3f",$10/$11)}' \
  "$O/cp_vs_nuclear.paf" | sort -k1,1 -k2,2n > "$O/nupt_regions.bed"
echo "NUPT candidate regions (>=300 bp): $(wc -l < "$O/nupt_regions.bed")"
awk '{s+=$3-$2} END{printf "NUPT bp total: %d (%.4f%% of 19.32 Gb)\n",s,100*s/19318822303}' "$O/nupt_regions.bed"
echo "identity distribution of NUPT hits:"
awk '{printf "%.1f\n",$4*10}' "$O/nupt_regions.bed" | sort -n | uniq -c | awk '{printf "  %.0f-%.0f%%  %d\n",$2*10,($2+1)*10,$1}'
