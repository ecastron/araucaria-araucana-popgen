#!/bin/bash
#SBATCH --job-name=refmap
#SBATCH --partition=basic
#SBATCH --cpus-per-task=12
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --array=1-8%2
#SBATCH --output=/home/ecastron/araucaria/logs/refmap_%A_%a.out
set -euo pipefail
module load Minimap2/minimap2
module load SAMTOOLS/samtools.1.20
D=/home/ecastron/araucaria/refs/araucana
R=/home/ecastron/araucaria/raw/rad
O=/home/ecastron/araucaria/results/refmap
LIST=/home/ecastron/araucaria/scripts/ref/selected_libs.txt
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | awk '{print $1}')
ST=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | awk '{print $2}')
echo "node $(hostname)  lib $SRR  stand $ST  start $(date)"

minimap2 -ax sr -t 12 "$D/v1.sr.mmi" "$R/${SRR}_1.fastq.gz" "$R/${SRR}_2.fastq.gz" \
  | samtools sort -@ 4 -m 2G -o "$O/${SRR}.bam" -
samtools index "$O/${SRR}.bam"
samtools flagstat "$O/${SRR}.bam" > "$O/${SRR}.flagstat"
echo "=== flagstat ==="; cat "$O/${SRR}.flagstat"

# q30 = confidently unique placements. If de novo paralog collapse was the cause of
# the flat allele balance, THIS arm is where a real 0.5 heterozygote mode must appear.
samtools mpileup -Q 20 -q 30 -d 2000 -f "$D/v1.fasta" "$O/${SRR}.bam" \
  | python3 /home/ecastron/araucaria/scripts/ref/balance_from_pileup.py \
      --label "${SRR} ${ST} MAPQ>=30" --depths "$O/${SRR}.q30.depths" \
  > "$O/${SRR}.q30.balance.txt"

# q0 = every placement including multi-mappers; the control arm. Paralogs that the
# reference DOES contain will pile up here and should look like the de novo result.
samtools mpileup -Q 20 -q 0 -d 2000 -f "$D/v1.fasta" "$O/${SRR}.bam" \
  | python3 /home/ecastron/araucaria/scripts/ref/balance_from_pileup.py \
      --label "${SRR} ${ST} MAPQ>=0" \
  > "$O/${SRR}.q0.balance.txt"

echo "=== q30 ==="; cat "$O/${SRR}.q30.balance.txt"
echo "done $(date)"
