#!/bin/bash
#SBATCH --job-name=arau-sweep
#SBATCH --partition=basic
#SBATCH --array=1-8%4
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/sweep_%A_%a.out
#SBATCH --error=/home/ecastron/araucaria/logs/sweep_%A_%a.err

# Step 2, stage 2: de novo assembly parameter sweep.
#
# One array task per value of M (mismatches allowed between stacks within an individual),
# with n = M as in Paris et al. 2017, and m = 3 (minimum stack depth) held fixed.
#
# The sweep runs on the 30 replicate libraries only: 15 trees, each sequenced as two
# independent libraries. That subset does double duty. It gives the r80 curve (loci shared
# by >=80% of samples) AND, because each tree appears twice, a direct genotyping-error rate
# per parameter setting, which is the criterion the analysis spec requires (Mastretta-Yanes
# et al. 2015). Running the grid on all 173 trees would cost far more and yield no error
# estimate.
#
# --paired makes gstacks assemble contigs from the randomly sheared R2 reads, which is the
# reason Stacks 2 was chosen for this library type.

set -uo pipefail

PROJ=/home/ecastron/araucaria
SAMPLES=$PROJ/rad/cleaned
POPMAP=$PROJ/metadata/popmap_replicates.tsv

module load Stacks/stacks-2.68

M=$SLURM_ARRAY_TASK_ID          # 1..8
n=$M                            # n = M, per Paris et al. 2017
m=3                             # minimum depth to create a stack

OUT=$PROJ/rad/sweep/M${M}
mkdir -p "$OUT"

if [[ -s "$OUT/populations.sumstats.tsv" ]]; then
  echo "SKIP  M=$M (already complete)"
  exit 0
fi

echo "START M=$M n=$n m=$m  $(date)"

# -X passes options through to the named stage. r80 is applied here so the curve can be
# read straight off populations.log for each M.
if ! denovo_map.pl \
      --samples "$SAMPLES" \
      --popmap "$POPMAP" \
      -o "$OUT" \
      -M "$M" -n "$n" -m "$m" \
      --paired \
      -T "$SLURM_CPUS_PER_TASK" \
      -X "populations: -r 0.80 --vcf --hwe" \
      > "$OUT/denovo_map.stdout" 2>&1; then
  echo "FAIL  M=$M denovo_map.pl failed" >&2
  tail -25 "$OUT/denovo_map.stdout" >&2
  exit 1
fi

if [[ ! -s "$OUT/populations.snps.vcf" ]]; then
  echo "FAIL  M=$M no VCF produced" >&2
  exit 1
fi

loci=$(grep -c "^[^#]" "$OUT/populations.snps.vcf" 2>/dev/null || echo 0)
echo "OK    M=$M n=$n m=$m r80_snps=$loci  $(date)"
