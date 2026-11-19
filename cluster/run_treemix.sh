#!/bin/bash
#SBATCH --job-name=arau-treemix
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/treemix_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/treemix_%j.err

# Q4: is there gene flow between the Coastal and Andean cordilleras?
#
# Everything so far INFERS gene flow rather than testing it: FST of 0.10-0.16 is far too low for
# lineages isolated since the Pleistocene, ADMIXTURE returns K = 1, and between-cordillera pairs
# are only +0.044 more differentiated than distance predicts. All consistent with restricted but
# continuing exchange -- none of it a direct test.
#
# TreeMix (Pickrell & Pritchard 2012) is the direct test. It fits a population tree from the
# covariance of allele frequencies, then adds migration edges where the tree alone cannot explain
# the data. If gene flow crosses the Central Valley, an edge should appear between cordilleras
# and the variance explained should jump when it does.
#
# Read the output on three things:
#   1. Does an edge connect Costa (PN, TR, VA) to an Andean stand?
#   2. How much residual variance does each successive edge remove? A large gain at m=1 that
#      flattens afterwards means one real exchange event, not noise being fitted.
#   3. Are the migration weights large enough to matter biologically?
#
# CAUTION recorded in advance: TreeMix assumes populations are discrete. ADMIXTURE says they are
# not -- structure here is clinal. Under isolation by distance TreeMix can fit migration edges
# that are really just the cline. So an edge is suggestive, NOT proof; corroboration is that the
# edge should cross the Central Valley rather than link geographic neighbours.

set -uo pipefail
module load Popgen/popgen-1.0
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh

P=/home/ecastron/araucaria
V=$P/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf
OUT=$P/results/Q4/treemix
mkdir -p "$OUT"; cd "$OUT"

command -v treemix >/dev/null || { echo "FATAL: treemix unavailable" >&2; exit 1; }

echo "=== 1. VCF -> per-stand allele counts ==="
conda activate ipyrad
python $P/scripts/vcf_to_treemix.py "$V" "$P/metadata/popmap_177.tsv" "$OUT/araucaria.treemix" || exit 1
gzip -f "$OUT/araucaria.treemix"
echo "  sites: $(( $(zcat "$OUT/araucaria.treemix.gz" | wc -l) - 1 ))"

echo
echo "=== 2. TreeMix, m = 0..5 migration edges ==="
# -root uses an Andean stand as the outgroup: no true outgroup was sequenced, so the tree is
# unrooted in substance and the root only fixes a drawing orientation.
# -k blocks SNPs to account for linkage; 500 is conservative for RAD loci, which are unlinked
# between loci but correlated within one.
for m in 0 1 2 3 4 5; do
  treemix -i araucaria.treemix.gz -m $m -k 500 -root LR -o "tm_m${m}" >/dev/null 2>&1 \
    && echo "  m=$m done" || echo "  m=$m FAILED"
done

echo
echo "=== 3. variance explained by each model ==="
conda activate dada2
Rscript $P/scripts/treemix_summary.R
echo "TREEMIX COMPLETE"
