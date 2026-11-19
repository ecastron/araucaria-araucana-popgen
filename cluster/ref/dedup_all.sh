#!/bin/bash
#SBATCH --job-name=dedupall
#SBATCH --partition=basic
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-8%3
#SBATCH --output=/home/ecastron/araucaria/logs/dedupall_%A_%a.out
set -euo pipefail
module load SAMTOOLS/samtools.1.20
O=/home/ecastron/araucaria/results/refmap
D=/home/ecastron/araucaria/refs/araucana
L=/home/ecastron/araucaria/scripts/ref/selected_libs.txt
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$L" | awk '{print $1}')
ST=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$L" | awk '{print $2}')
[ -f "$O/$SRR.bam" ] || { echo "no BAM for $SRR yet"; exit 1; }
cd "$O"
[ -f "$SRR.nodup.bam" ] || {
  samtools collate -@ 8 -O -u "$SRR.bam" tmp_col_$SRR \
    | samtools fixmate -@ 8 -m -u - - \
    | samtools sort -@ 8 -m 2G -u - \
    | samtools markdup -@ 8 -f "$SRR.markdup.stats" - "$SRR.dedup.bam"
  samtools view -@ 4 -b -F 1024 "$SRR.dedup.bam" > "$SRR.nodup.bam"
  samtools index "$SRR.nodup.bam"
}
echo "### $SRR ($ST) duplicate stats"
grep -E "DUPLICATE TOTAL|ESTIMATED_LIBRARY_SIZE|READ:|WRITTEN" "$SRR.markdup.stats" || cat "$SRR.markdup.stats"
samtools mpileup -Q 20 -q 30 -d 2000 -f "$D/v1.fasta" "$SRR.nodup.bam" \
  | python3 /home/ecastron/araucaria/scripts/ref/balance_from_pileup.py \
      --label "$SRR $ST dedup MAPQ>=30" > "$SRR.q30.nodup.balance.txt"
sed -n '1,14p' "$SRR.q30.nodup.balance.txt"
