#!/bin/bash
#SBATCH -J qc_rev
#SBATCH --partition=basic
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH -o /home/ecastron/araucaria/logs/qc__rev%j.out
# Revision (25 Nov 2026): per-tree sequencing and mapping QC for the 163 analyzed trees.
# Columns: tree, stand, reads_raw, reads_passed_filter, clusters_hidepth, reads_consens,
#   loci_in_assembly (ipyrad full_ct091_stats.txt, "Final Sample stats summary"),
#   call_rate (share of the 78,026 SNPs called; results/revision/callrate_by_tree.tsv),
#   bt2_overall_alignment (bowtie2 --very-sensitive to AraucariaDB v1, results/refmap177/<run>.bt2.stats),
#   properly_paired_pct and mapq20_pct (samtools flagstat / view -q 20 on results/refmap177/<run>.bam,
#   share of mapped reads).
# PCR duplicates are not reported: single-digest RAD reads all start at the PstI site, so
# position-based duplicate marking cannot separate PCR copies from independent reads.
# Output: results/revision/qc_per_tree.tsv
set -euo pipefail
P=/home/ecastron/araucaria; OUT=$P/results/revision
ST=$P/rad/ipyrad_full/full_ct091_outfiles/full_ct091_stats.txt
module load SAMTOOLS/samtools.1.20
command -v samtools >/dev/null || { echo "FATAL: samtools unavailable" >&2; exit 1; }
awk '/^## Final Sample stats summary/{f=1; getline; next} /^## /{f=0} f && NF==10 {sub(/\.[12]$/,"",$1); print $1"\t"$3"\t"$4"\t"$6"\t"$9"\t"$10}' $ST > $OUT/ipyrad_s7.tmp
printf "tree\tstand\treads_raw\treads_passed_filter\tclusters_hidepth\treads_consens\tloci_in_assembly\tcall_rate\tbt2_overall_alignment_pct\tproperly_paired_pct\tmapq20_pct\n" > $OUT/qc_per_tree.tsv
tail -n +2 $OUT/callrate_by_tree.tsv | while IFS=$'\t' read -r t s r cr; do
  ip=$(awk -v t="$t" '$1==t{print $2"\t"$3"\t"$4"\t"$5"\t"$6}' $OUT/ipyrad_s7.tmp)
  aln=$(grep -h "overall alignment rate" $P/results/refmap177/$t.bt2.stats | awk '{sub("%","",$1); print $1}')
  B=$P/results/refmap177/$t.bam
  fs=$(samtools flagstat -@ 8 $B)
  mapped=$(echo "$fs" | awk '/ mapped \(/{print $1; exit}')
  pp=$(echo "$fs" | awk '/properly paired/{gsub(/[(%]/,"",$6); print $6; exit}')
  q20=$(samtools view -@ 8 -c -F 4 -q 20 $B)
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%.4f\n" "$t" "$s" "$ip" "$cr" "$aln" "$pp" "$(awk -v a="$q20" -v b="$mapped" 'BEGIN{print 100*a/b}')" >> $OUT/qc_per_tree.tsv
done
rm -f $OUT/ipyrad_s7.tmp
echo "rows: $(($(wc -l < $OUT/qc_per_tree.tsv)-1))"
