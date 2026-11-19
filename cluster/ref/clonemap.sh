#!/bin/bash
#SBATCH --job-name=clonemap
#SBATCH --partition=basic
#SBATCH --cpus-per-task=12
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --array=1-4%2
#SBATCH --output=/home/ecastron/araucaria/logs/clonemap_%A_%a.out
set -euo pipefail
module load Minimap2/minimap2
module load SAMTOOLS/samtools.1.20
D=/home/ecastron/araucaria/refs/araucana
O=/home/ecastron/araucaria/results/refmap
C=/home/ecastron/araucaria/results/clonefilter
L=/home/ecastron/araucaria/scripts/ref/selected_libs.txt
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$L" | awk '{print $1}')
ST=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$L" | awk '{print $2}')
R1=$(ls $C/$SRR/out/*_1.1.fq.gz); R2=$(ls $C/$SRR/out/*_2.2.fq.gz)
echo "$SRR ($ST) mapping SEQUENCE-deduplicated reads"
minimap2 -ax sr -t 12 "$D/v1.sr.mmi" "$R1" "$R2" \
  | samtools sort -@ 4 -m 2G -o "$O/$SRR.clonefilt.bam" -
samtools index "$O/$SRR.clonefilt.bam"
samtools mpileup -Q 20 -q 30 -d 2000 -f "$D/v1.fasta" "$O/$SRR.clonefilt.bam" \
  | python3 /home/ecastron/araucaria/scripts/ref/balance_from_pileup.py \
      --label "$SRR $ST clone_filter(sequence) MAPQ>=30" \
  > "$O/$SRR.q30.clonefilt.balance.txt"
sed -n '1,14p' "$O/$SRR.q30.clonefilt.balance.txt"
