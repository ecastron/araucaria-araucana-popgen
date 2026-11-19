#!/bin/bash
#SBATCH --job-name=fbcall
#SBATCH --partition=basic
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=72:00:00
#SBATCH --array=0-199%24
#SBATCH --output=/home/ecastron/araucaria/logs/fbcall_%A_%a.out
set -eo pipefail
FB=/home/ecastron/envs/freebayes/bin/freebayes
D=/home/ecastron/araucaria/refs/araucana
S=/home/ecastron/araucaria/results/sv
O=$S/fb_out; mkdir -p "$O"
C=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
BED=$S/fb_chunks/chunk_$C
OUT=$O/fb_$C.vcf
[ -s "$OUT" ] && { echo "chunk $C already done"; exit 0; }
ls /home/ecastron/araucaria/results/refmap177/*.bam > $S/bamlist177.txt
n=$(wc -l < $S/bamlist177.txt)
echo "node $(hostname) chunk $C | $n BAMs | $(wc -l < $BED) regions | start $(date)"
ulimit -n 4096 || true
# --pooled-continuous: treat each sample as a pool, report continuous allele
#   frequencies rather than forcing a diploid genotype -- appropriate given the
#   documented unreliability of diploid calls in this dataset.
# --use-best-n-alleles 4: cap multiallelic sites at 4 alleles.
$FB -f "$D/v1.fasta" -L "$S/bamlist177.txt" -t "$BED" \
    --use-best-n-alleles 4 --pooled-continuous \
    --min-alternate-count 2 --min-mapping-quality 20 --min-base-quality 20 \
  > "$OUT"
echo "chunk $C variants: $(grep -vc '^#' "$OUT")  done $(date)"
