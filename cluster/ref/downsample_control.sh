#!/bin/bash
#SBATCH --job-name=dsctrl
#SBATCH --partition=basic
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/dsctrl_%j.out
set -euo pipefail
module load SAMTOOLS/samtools.1.20
D=/home/ecastron/araucaria/refs/araucana
O=/home/ecastron/araucaria/results/refmap
cd "$O"
SRR=SRR11841753    # BP, the library with all three treatments already computed
# Random removal ignores duplicate structure. If it reproduces markdup's apparent
# rescue, the rescue is a depth artifact (smaller n -> larger expected SD -> ratio
# pulled toward 1), NOT evidence that duplicates caused the spread.
for FRAC in 0.111 0.365; do
  TAG=$(echo $FRAC | tr -d '.')
  echo "############ random downsample to ${FRAC} of reads ############"
  samtools view -@ 8 -b -s 42${FRAC#0} "$SRR.bam" > "$SRR.rand$TAG.bam"
  samtools index "$SRR.rand$TAG.bam"
  echo "reads retained: $(samtools view -c -F 0x900 -F 4 $SRR.rand$TAG.bam)"
  samtools mpileup -Q 20 -q 30 -d 2000 -f "$D/v1.fasta" "$SRR.rand$TAG.bam" \
    | python3 /home/ecastron/araucaria/scripts/ref/balance_from_pileup.py \
        --label "$SRR BP RANDOM ${FRAC} of reads MAPQ>=30" \
    > "$SRR.q30.rand$TAG.balance.txt"
  sed -n '1,13p' "$SRR.q30.rand$TAG.balance.txt"
  echo
done
echo "### reference points for the same library ###"
echo "markdup (11% kept):"; awk '/^[0-9]/&&NF>=9{printf "  %-8s %s\n",$1,$8}' "$SRR.q30.nodup.balance.txt"
echo "clone_filter (36% kept):"; awk '/^[0-9]/&&NF>=9{printf "  %-8s %s\n",$1,$8}' "$SRR.q30.clonefilt.balance.txt"
echo "all reads:"; awk '/^[0-9]/&&NF>=9{printf "  %-8s %s\n",$1,$8}' "$SRR.q30.balance.txt"
