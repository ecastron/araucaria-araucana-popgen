#!/bin/bash
#SBATCH --job-name=bt2map
#SBATCH --partition=basic
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --array=1-8%2
#SBATCH --output=/home/ecastron/araucaria/logs/bt2map_%A_%a.out
export JAVA_HOME="${JAVA_HOME:-}"
set -eo pipefail
module load Bowtie2/bowtie2
module load SAMTOOLS/samtools.1.20
D=/home/ecastron/araucaria/refs/araucana
O=/home/ecastron/araucaria/results/refmap
R=/home/ecastron/araucaria/raw/rad
L=/home/ecastron/araucaria/scripts/ref/selected_libs.txt
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$L" | awk '{print $1}')
ST=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$L" | awk '{print $2}')
echo "node $(hostname)  lib $SRR ($ST)  start $(date)"
# --no-unal keeps the BAM small; XS:i is emitted whenever a second valid alignment
# exists, which is the read-level multi-mapping readout bowtie2 adds over minimap2.
bowtie2 -p 16 -x "$D/v1_bt2" --very-sensitive --no-unal \
        -1 "$R/${SRR}_1.fastq.gz" -2 "$R/${SRR}_2.fastq.gz" \
        2> "$O/${SRR}.bt2.stats" \
  | samtools sort -@ 4 -m 2G -o "$O/${SRR}.bt2.bam" -
samtools index "$O/${SRR}.bt2.bam"
echo "=== bowtie2 alignment summary ==="; cat "$O/${SRR}.bt2.stats"
samtools flagstat "$O/${SRR}.bt2.bam" > "$O/${SRR}.bt2.flagstat"
echo "=== multi-mapping reads (XS present = a second valid alignment exists) ==="
tot=$(samtools view -c -F 0x900 "$O/${SRR}.bt2.bam")
xs=$(samtools view -F 0x900 "$O/${SRR}.bt2.bam" | grep -c "XS:i:" || true)
awk -v t="$tot" -v x="$xs" 'BEGIN{printf "  primary %d | with XS %d (%.1f%%)\n",t,x,100*x/t}'
echo "=== allele balance, same pipeline as minimap2/bwa-mem2 ==="
samtools mpileup -Q 20 -q 30 -d 2000 -f "$D/v1.fasta" "$O/${SRR}.bt2.bam" \
  | python3 /home/ecastron/araucaria/scripts/ref/balance_from_pileup.py \
      --label "$SRR $ST bowtie2 --very-sensitive MAPQ>=30" > "$O/${SRR}.q30.bt2.balance.txt"
sed -n '1,14p' "$O/${SRR}.q30.bt2.balance.txt"
echo "done $(date)"
