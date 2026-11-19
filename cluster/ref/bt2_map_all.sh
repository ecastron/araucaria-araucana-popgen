#!/bin/bash
#SBATCH --job-name=bt2all
#SBATCH --partition=basic
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=72:00:00
#SBATCH --array=1-177%6
#SBATCH --output=/home/ecastron/araucaria/logs/bt2all_%A_%a.out
export JAVA_HOME="${JAVA_HOME:-}"
set -eo pipefail
module load Bowtie2/bowtie2
module load SAMTOOLS/samtools.1.20
D=/home/ecastron/araucaria/refs/araucana
O=/home/ecastron/araucaria/results/refmap177
R=/home/ecastron/araucaria/raw/rad
L=/home/ecastron/araucaria/scripts/ref/all177.txt
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$L" | cut -f1)
ST=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$L" | cut -f2)
[ -s "$O/${SRR}.bam" ] && { echo "$SRR already done"; exit 0; }
echo "node $(hostname) lib $SRR ($ST) start $(date)"
# read groups so FreeBayes can keep the 177 samples separate from one merged call
bowtie2 -p 16 -x "$D/v1_bt2" --very-sensitive --no-unal \
        --rg-id "$SRR" --rg "SM:${SRR}" --rg "PL:ILLUMINA" --rg "LB:${SRR}" \
        -1 "$R/${SRR}_1.fastq.gz" -2 "$R/${SRR}_2.fastq.gz" 2> "$O/${SRR}.bt2.stats" \
  | samtools sort -@ 4 -m 2G -o "$O/${SRR}.bam" -
samtools index "$O/${SRR}.bam"
grep "overall alignment rate" "$O/${SRR}.bt2.stats" | sed "s/^/  $SRR /"
echo "done $(date)"
