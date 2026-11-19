#!/bin/bash
#SBATCH --job-name=clonefilt
#SBATCH --partition=basic
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=1-8%4
#SBATCH --output=/home/ecastron/araucaria/logs/clonefilt_%A_%a.out
set -euo pipefail
module load Stacks/stacks-2.68
R=/home/ecastron/araucaria/raw/rad
O=/home/ecastron/araucaria/results/clonefilter
L=/home/ecastron/araucaria/scripts/ref/selected_libs.txt
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$L" | awk '{print $1}')
ST=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$L" | awk '{print $2}')
W=$O/$SRR; mkdir -p "$W/in" "$W/out"
# fasterq-dump deflines carry a 'length=' field that DIFFERS between mates, which
# breaks Stacks' pair matching. Truncate to the first token (spec, RAD section).
for m in 1 2; do
  zcat "$R/${SRR}_${m}.fastq.gz" \
    | awk 'NR%4==1{split($0,a," "); print a[1]; next} {print}' \
    | gzip -1 > "$W/in/${SRR}_${m}.fq.gz"
done
echo "### $SRR ($ST) clone_filter -- SEQUENCE-based, no reference, no aligner"
clone_filter -1 "$W/in/${SRR}_1.fq.gz" -2 "$W/in/${SRR}_2.fq.gz" \
  -i gzfastq -o "$W/out" 2>&1 | tail -20
rm -rf "$W/in"
