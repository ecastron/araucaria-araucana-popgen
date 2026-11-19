#!/bin/bash
#SBATCH --job-name=bwacheck
#SBATCH --partition=basic
#SBATCH --cpus-per-task=24
#SBATCH --mem=900G          # >257G pins to n001 (2 TB); bwa-mem2 index ~540 GB for 19.3 Gb
#SBATCH --time=72:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/bwacheck_%j.out
set -euo pipefail
module load Bwa-mem2/bwa-mem2
module load SAMTOOLS/samtools.1.20
D=/home/ecastron/araucaria/refs/araucana
O=/home/ecastron/araucaria/results/refmap
R=/home/ecastron/araucaria/raw/rad
SRR=SRR11841615
cd "$D"
echo "node $(hostname) start $(date)"
if [ ! -f v1.fasta.bwt.2bit.64 ]; then
  echo ">>> building bwa-mem2 index (slow, one-off)"
  bwa-mem2 index v1.fasta
fi
echo ">>> index done $(date)"; ls -lh v1.fasta.bwt.2bit.64 2>/dev/null || true
bwa-mem2 mem -t 24 v1.fasta "$R/${SRR}_1.fastq.gz" "$R/${SRR}_2.fastq.gz" \
  | samtools sort -@ 6 -m 4G -o "$O/${SRR}.bwa.bam" -
samtools index "$O/${SRR}.bwa.bam"
samtools flagstat "$O/${SRR}.bwa.bam" > "$O/${SRR}.bwa.flagstat"
echo "############ BWA-MEM2 vs MINIMAP2 mapping rate ############"
echo "--- bwa-mem2 ---"; grep -E "in total|[^y] mapped \(|properly paired" "$O/${SRR}.bwa.flagstat"
echo "--- minimap2 ---"; grep -E "in total|[^y] mapped \(|properly paired" "$O/${SRR}.flagstat"
echo
echo "############ does the balance artifact survive bwa-mem2? ############"
samtools mpileup -Q 20 -q 30 -d 2000 -f v1.fasta "$O/${SRR}.bwa.bam" \
  | python3 /home/ecastron/araucaria/scripts/ref/balance_from_pileup.py \
      --label "$SRR VA bwa-mem2 MAPQ>=30" > "$O/${SRR}.bwa.q30.balance.txt"
sed -n '1,14p' "$O/${SRR}.bwa.q30.balance.txt"
echo "done $(date)"
