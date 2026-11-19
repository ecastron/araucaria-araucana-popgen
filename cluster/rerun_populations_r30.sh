#!/bin/bash
#SBATCH --job-name=arau-repop
#SBATCH --partition=basic
#SBATCH --array=1-8%4
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/repop_%A_%a.out
#SBATCH --error=/home/ecastron/araucaria/logs/repop_%A_%a.err

# Re-run only the populations stage of the sweep with a single-population popmap.
#
# Stacks' populations -r is a WITHIN-population threshold. The sweep popmap grouped the 30
# libraries by site, giving 8 populations of 2-4, so a locus present in two libraries at one
# site satisfied -r 0.30 for that population and was exported. The result was tens of
# thousands of "r80" SNPs with a median per-library call rate of 123 sites, which made the
# replicate-concordance estimates noise.
#
# With one population, -r 0.30 means called in >=80% of all 30 libraries, which is the rule
# Paris et al. describe. ustacks/cstacks/sstacks/gstacks output is unaffected and reused.

set -uo pipefail
PROJ=/home/ecastron/araucaria
M=$SLURM_ARRAY_TASK_ID
D=$PROJ/rad/sweep/M${M}

module load Stacks/stacks-2.68
[[ -f "$D/catalog.calls" ]] || { echo "FAIL M=$M no catalog.calls" >&2; exit 1; }

if ! populations -P "$D" -M $PROJ/metadata/popmap_replicates_onepop.tsv \
      -r 0.30 --vcf --threads "$SLURM_CPUS_PER_TASK" > "$D/populations_r30.log" 2>&1; then
  echo "FAIL M=$M populations failed" >&2; tail -15 "$D/populations_r30.log" >&2; exit 1
fi

SNPS=$(grep -vc '^#' "$D/populations.snps.vcf")
echo "OK M=$M r80_snps=$SNPS"
