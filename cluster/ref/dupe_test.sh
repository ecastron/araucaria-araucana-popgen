#!/bin/bash
#SBATCH --job-name=duptest
#SBATCH --partition=basic
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/duptest_%j.out
set -euo pipefail
module load SAMTOOLS/samtools.1.20
O=/home/ecastron/araucaria/results/refmap
D=/home/ecastron/araucaria/refs/araucana
SRR=SRR11841615
cd "$O"
# R2 is unanchored (random sonication shear points), so a PCR duplicate is
# identifiable as a pair sharing BOTH ends. De novo assembly could never see this.
samtools collate -@ 8 -O -u "$SRR.bam" tmp_collate_$SRR \
  | samtools fixmate -@ 8 -m -u - - \
  | samtools sort -@ 8 -m 2G -u - \
  | samtools markdup -@ 8 -f "$SRR.markdup.stats" - "$SRR.dedup.bam"
samtools index "$SRR.dedup.bam"
echo "############ DUPLICATE STATS ############"
cat "$SRR.markdup.stats"
echo
echo "############ BALANCE AFTER REMOVING DUPLICATES ############"
# -F 1024 drops reads flagged as duplicates
samtools view -@ 4 -b -F 1024 "$SRR.dedup.bam" > "$SRR.nodup.bam"
samtools index "$SRR.nodup.bam"
samtools mpileup -Q 20 -q 30 -d 2000 -f "$D/v1.fasta" "$SRR.nodup.bam" \
  | python3 /home/ecastron/araucaria/scripts/ref/balance_from_pileup.py \
      --label "$SRR VA MAPQ>=30 DUPLICATES REMOVED" \
  > "$SRR.q30.nodup.balance.txt"
sed -n '1,14p' "$SRR.q30.nodup.balance.txt"
