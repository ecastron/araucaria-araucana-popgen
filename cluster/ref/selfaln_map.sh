#!/bin/bash
#SBATCH --job-name=selfaln
#SBATCH --partition=basic
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --array=1-40%4
#SBATCH --output=/home/ecastron/araucaria/logs/selfaln_%A_%a.out
set -euo pipefail
module load Minimap2/minimap2
module load SAMTOOLS/samtools.1.20
D=/home/ecastron/araucaria/refs/araucana
R=/home/ecastron/araucaria/results/sv/selfaln
i=$(printf "%03d" "$SLURM_ARRAY_TASK_ID")
Q="$R/chunk_$i.fa"
echo "node $(hostname) task $i start $(date)"
# extract this chunk's contigs
samtools faidx "$D/v1.fasta" -r "$R/chunk_$i.txt" > "$Q"
# No -c: base-level alignment is unaffordable at 19.3 Gb and PAF cols 10/11
# are populated regardless, which is all the identity filter needs.
# Self-hits (same contig, overlapping interval) are dropped here, not by -X,
# so the behaviour is explicit rather than relying on ava-mode semantics.
minimap2 -x asm20 -t 8 -N 20 -p 0.2 --secondary=yes "$D/v1.asm20.mmi" "$Q" 2> "$R/mm_$i.log" \
 | awk -F'\t' 'BEGIN{OFS="\t"}
     {
       if($11<1000) next;
       if($1==$6 && $3<$9 && $8<$4) next;   # self-overlap on the same contig
       # WITHOUT -c, PAF col10 is a chain-derived estimate and col10/col11 badly
       # understates identity -- filtering on it discarded nearly everything.
       # minimap2 emits dv:f (approximate per-base divergence) for approximate
       # mappings; use that. Keep records where the tag is absent.
       dv=-1;
       for(i=13;i<=NF;i++) if($i ~ /^dv:f:/){ split($i,a,":"); dv=a[3]+0; break }
       if(dv>=0 && dv>0.10) next;
       print $1,$3,$4,$5,$6,$8,$9,$10,$11,dv }' \
 | gzip -1 > "$R/self_$i.paf.gz"
rm -f "$Q"
echo "task $i alignments kept: $(zcat "$R/self_$i.paf.gz" | wc -l)"
echo "done $(date)"
