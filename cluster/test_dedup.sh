#!/bin/bash
#SBATCH --job-name=arau-dedup
#SBATCH --partition=basic
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/dedup_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/dedup_%j.err
set -uo pipefail
PROJ=/home/ecastron/araucaria
module load Stacks/stacks-2.68
SRC=$PROJ/rad/sweep/M5
DST=$PROJ/rad/dedup_test
rm -rf "$DST"; mkdir -p "$DST"
cp "$SRC"/*.bam "$DST"/ 2>/dev/null
cp "$SRC"/catalog.fa.gz "$DST"/ 2>/dev/null
echo "BAMs: $(ls $DST/*.bam | wc -l)"
gstacks -P "$DST" -M $PROJ/metadata/popmap_replicates_onepop.tsv \
        --rm-pcr-duplicates --threads "$SLURM_CPUS_PER_TASK" > "$DST/gstacks.stdout" 2>&1 || {
  echo "gstacks failed" >&2; tail -20 "$DST/gstacks.stdout" >&2; exit 1; }
grep -iE "duplicate|coverage" "$DST/gstacks.log" | head -8
populations -P "$DST" -M $PROJ/metadata/popmap_replicates_onepop.tsv \
        -r 0.30 --vcf --threads "$SLURM_CPUS_PER_TASK" > "$DST/populations.stdout" 2>&1 || {
  echo "populations failed" >&2; tail -10 "$DST/populations.stdout" >&2; exit 1; }
echo "r30 sites: $(grep -vc '^#' $DST/populations.snps.vcf)"
echo DEDUP_TEST_DONE
