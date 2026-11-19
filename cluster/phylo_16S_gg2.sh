#!/bin/bash
#SBATCH --job-name=arau-gg2
#SBATCH --partition=basic
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G   # n002/n003 cap at 257G; anything above pins the job to n001
#SBATCH --time=48:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/phylo_gg2_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/phylo_gg2_%j.err

# Re-place the 16S ASVs against Greengenes2 2022.10 instead of Greengenes 13_8.
#
# Why 2022.10 and not 2024.09: only 2022.10 publishes a SEPP reference package. The 2024.09
# release ships phylogeny and taxonomy artefacts but no `backbone.sepp-reference.qza` (verified
# 404), so fragment insertion is not possible against it. 2022.10 is therefore the most recent
# reference compatible with placement, and is still nine years newer than GG 13_8.
#
# The upgrade is not merely a date. GG 13_8 is a 16S-only tree; Greengenes2 inserts 16S into a
# whole-genome (Web of Life) backbone, so the deep structure our 253 bp fragments cannot resolve
# comes from genomes rather than from the same short marker (McDonald et al. 2023).
#
# The GG 13_8 tree is kept, not discarded: comparing the two quantifies how much the choice of
# reference moves the result, which is the same sensitivity logic already applied to de novo.

set -euo pipefail
PROJ=/home/ecastron/araucaria
OUT=$PROJ/microbiome/16S/phylo
REFS=$PROJ/refs
cd "$OUT"

source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.5

FA=$PROJ/microbiome/16S/asv_16S.fasta
REF=$REFS/gg2-2022.10-sepp-reference.qza
if [[ ! -s "$REF" ]]; then
  echo "downloading Greengenes2 2022.10 SEPP reference"
  curl -fSL --retry 3 -o "$REF" \
    http://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.sepp-reference.qza
fi
echo "reference: $(du -h "$REF" | cut -f1)"

echo "=== SEPP insertion against Greengenes2 2022.10 ==="; date
qiime fragment-insertion sepp \
  --i-representative-sequences rep-seqs.qza \
  --i-reference-database "$REF" \
  --p-threads "$SLURM_CPUS_PER_TASK" \
  --o-tree insertion-tree-gg2.qza \
  --o-placements insertion-placements-gg2.qza
date

qiime tools export --input-path insertion-tree-gg2.qza --output-path exported-gg2
mv exported-gg2/tree.nwk asv_16S_sepp_gg2_raw.nwk

# Prune to our own features. The reference supplies topology and branch lengths; the tree that
# is used downstream must have exactly the tips the feature table has.
Rscript -e '
  suppressPackageStartupMessages(library(ape))
  tr  <- read.tree("asv_16S_sepp_gg2_raw.nwk")
  ids <- sub("^>", "", grep("^>", readLines("'"$FA"'"), value = TRUE))
  hit <- intersect(tr$tip.label, ids)
  cat(sprintf("raw tips: %d | our ASVs present: %d / %d | unplaced: %d\n",
              Ntip(tr), length(hit), length(ids), length(ids) - length(hit)))
  write.tree(keep.tip(tr, hit), "asv_16S_sepp_gg2_pruned.nwk")
  cat(sprintf("pruned tree tips: %d\n", Ntip(read.tree("asv_16S_sepp_gg2_pruned.nwk"))))'

echo "DONE"; ls -la "$OUT"/*.nwk
