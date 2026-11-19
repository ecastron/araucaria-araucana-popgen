#!/bin/bash
#SBATCH --job-name=arau-Q2
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=150G
#SBATCH --time=72:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/Q2_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/Q2_%j.err

# Q2: at what spatial scale is kinship structured? Reports Sp and Nb.
#
# Kinship comes from GENOTYPE LIKELIHOODS, not called genotypes. This dataset has no technical
# replicates (the fifteen paired tree_ids are related trees), so the genotyping error rate cannot
# be measured directly, and low-coverage sites are demonstrably under-called. Likelihoods carry
# that uncertainty into the estimate rather than discarding it.
#
# Sigma is NOT computed. Converting Nb to sigma needs effective density (Nb = 4*pi*De*sigma^2),
# which was never measured. Sp is the quantity built for cross-study comparison and carries the
# result on its own.
#
# The AMOVA precondition is met: 84.2% of genetic variance sits WITHIN stands (p = 1e-4), so the
# within-stand contrasts Q2 depends on are identifiable.

set -uo pipefail
module load Popgen/popgen-1.0
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh

P=/home/ecastron/araucaria
V=$P/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf
OUT=$P/results/Q2
mkdir -p "$OUT"

command -v ngsRelate >/dev/null || { echo "FATAL: ngsRelate unavailable" >&2; exit 1; }

echo "=== 1. CATG read counts -> genotype likelihoods (Beagle) ==="
conda activate ipyrad
python $P/scripts/catg_to_beagle.py "$V" "$OUT/gl.beagle" 0.01 || exit 1
gzip -f "$OUT/gl.beagle"
NS=$(zcat "$OUT/gl.beagle.gz" | head -1 | awk '{print (NF-3)/3}')
NL=$(( $(zcat "$OUT/gl.beagle.gz" | wc -l) - 1 ))
echo "  $NL sites x $NS individuals"

echo
echo "=== 2. ngsRelate pairwise relatedness ==="
zcat "$OUT/gl.beagle.gz" | head -1 | cut -f4- | tr '\t' '\n' | awk 'NR%3==1' > "$OUT/samples.txt"
# -L is required: ngsRelate will not infer the site count from the beagle file itself.
ngsRelate -G "$OUT/gl.beagle.gz" -n "$NS" -L "$NL" -z "$OUT/samples.txt" \
  -O "$OUT/ngsrelate.tsv" -p 16 2>&1 | tail -20
[ -s "$OUT/ngsrelate.tsv" ] || { echo "FATAL: ngsRelate produced nothing" >&2; exit 1; }
echo "  pairs estimated: $(( $(wc -l < "$OUT/ngsrelate.tsv") - 1 ))"

echo
echo "=== 3. Sp / Nb and the kinship-distance regression ==="
conda activate dada2
Rscript $P/scripts/q2_sgs.R
echo "Q2 COMPLETE"
