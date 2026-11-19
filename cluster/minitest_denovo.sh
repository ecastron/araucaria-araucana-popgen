#!/bin/bash
#SBATCH --job-name=arau-minitest
#SBATCH --partition=basic
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/minitest_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/minitest_%j.err

# End-to-end validation of the denovo_map chain on a tiny subset, BEFORE recommitting the
# 8-way parameter grid.
#
# The previous sweep failed at tsv2bam - the last stage - roughly two hours into each run,
# because SRA deflines carry a length= field that differs between mates. The deflines are now
# truncated, but that fix is unproven until tsv2bam actually pairs reads successfully. This
# job runs the complete chain (ustacks -> cstacks -> sstacks -> tsv2bam -> gstacks ->
# populations) on 4 libraries x 100k read pairs, which finishes in minutes rather than hours.
#
# Two replicate PAIRS are used, not four arbitrary libraries, so tsv2bam has genuine
# paired-end matching to do and the catalog has something to match across samples.

set -uo pipefail

PROJ=/home/ecastron/araucaria
SRC=$PROJ/rad/cleaned
OUT=$PROJ/rad/minitest
SUB=100000        # read pairs per library

module load Stacks/stacks-2.68

rm -rf "$OUT"; mkdir -p "$OUT/samples"

# CABP07 and CAHL16: two trees each sequenced as two independent libraries.
SAMPLES=(SRR11841691 SRR11841700 SRR11841648 SRR11841655)

for s in "${SAMPLES[@]}"; do
  for m in 1 2; do
    f=$SRC/${s}.${m}.fq.gz
    [[ -s "$f" ]] || { echo "FAIL missing $f" >&2; exit 1; }
    # Reads are in identical order in both mates, so taking the first N from each keeps
    # them synchronised without needing a subsampler.
    zcat "$f" | head -n $((SUB * 4)) | gzip > "$OUT/samples/${s}.${m}.fq.gz"
  done
  # Confirm the deflines actually match across mates before Stacks sees them.
  a=$(zcat "$OUT/samples/${s}.1.fq.gz" | awk 'NR%4==1{print $1; exit}')
  b=$(zcat "$OUT/samples/${s}.2.fq.gz" | awk 'NR%4==1{print $1; exit}')
  if [[ "${a%/1}" != "${b%/2}" ]]; then
    echo "FAIL  $s deflines still mismatched: R1=$a R2=$b" >&2
    exit 1
  fi
  echo "  prepared $s (R1=$a R2=$b)"
done

cat > "$OUT/popmap.tsv" <<'EOF'
SRR11841691	BP
SRR11841700	BP
SRR11841648	HL
SRR11841655	HL
EOF

echo "=== running denovo_map.pl ==="
if ! denovo_map.pl --samples "$OUT/samples" --popmap "$OUT/popmap.tsv" -o "$OUT" \
       -M 2 -n 2 -m 3 --paired -T "$SLURM_CPUS_PER_TASK" \
       -X "populations: -r 0.80 --vcf" > "$OUT/denovo_map.stdout" 2>&1; then
  echo "MINITEST FAILED" >&2
  tail -25 "$OUT/denovo_map.stdout" >&2
  echo "--- denovo_map.log tail ---" >&2
  tail -25 "$OUT/denovo_map.log" >&2
  exit 1
fi

if [[ ! -s "$OUT/populations.snps.vcf" ]]; then
  echo "MINITEST FAILED: no VCF produced" >&2
  exit 1
fi

snps=$(grep -vc '^#' "$OUT/populations.snps.vcf")
# The decisive line: tsv2bam reports how many paired reads it actually placed.
echo "=== tsv2bam summary ==="
grep -iE "paired|matching" "$OUT/denovo_map.log" | tail -6
echo "MINITEST PASSED  r80_snps=$snps"
