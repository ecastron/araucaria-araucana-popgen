#!/bin/bash
#SBATCH --job-name=arau-phylo
#SBATCH --partition=basic
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --time=24:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/phylo16S_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/phylo16S_%j.err

# Phylogeny for the 16S ASVs, by SEPP fragment insertion into a curated reference tree.
#
# 253 bp of V4 does not carry enough signal to reconstruct its own tree. Janssen et al. (2018,
# mSystems) benchmarked exactly this case -- 16S V4 fragments from a denoiser -- and showed a
# de novo tree gives incorrect results, while inserting the fragments into a reference
# phylogeny does not. We therefore take the reference-insertion tree as primary and build the
# de novo tree only as a sensitivity check, so that any claim resting on UniFrac or Faith's PD
# can be shown not to hinge on which tree was used.
#
# ITS gets NO tree, deliberately. ITS is not alignable across fungi at the depth these ASVs
# span (121-461 bp here), and phylogenetic placement of fungal metabarcodes needs full-length
# ITS from long reads (Tedersoo 2018; Furneaux 2020). Q7 is CLR-primary and needs no tree.

set -euo pipefail
PROJ=/home/ecastron/araucaria
AMP=16S
OUT=$PROJ/microbiome/$AMP/phylo
REFS=$PROJ/refs
mkdir -p "$OUT" "$REFS"

source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.5

FA=$PROJ/microbiome/$AMP/asv_${AMP}.fasta
[[ -s "$FA" ]] || { echo "no ASV fasta at $FA" >&2; exit 1; }
echo "ASVs to place: $(grep -c '^>' "$FA")"

# ---- reference package -------------------------------------------------------------------
REF=$REFS/sepp-refs-gg-13-8.qza
if [[ ! -s "$REF" ]]; then
  echo "downloading SEPP reference"
  curl -fsSL -o "$REF" \
    https://data.qiime2.org/classifiers/sepp-ref-dbs/sepp-refs-gg-13-8.qza
fi

# ---- primary: fragment insertion ---------------------------------------------------------
cd "$OUT"
if [[ ! -s rep-seqs.qza ]]; then
  qiime tools import --input-path "$FA" --output-path rep-seqs.qza \
    --type 'FeatureData[Sequence]'
fi

if [[ ! -s insertion-tree.qza ]]; then
  echo "=== SEPP insertion ==="; date
  qiime fragment-insertion sepp \
    --i-representative-sequences rep-seqs.qza \
    --i-reference-database "$REF" \
    --p-threads "$SLURM_CPUS_PER_TASK" \
    --o-tree insertion-tree.qza \
    --o-placements insertion-placements.qza
  date
fi

# Fragments too divergent to place are silently absent from the tree. Recording which ones
# were dropped matters: they are real ASVs that phylogenetic metrics will ignore, so the
# fraction lost is a caveat on every UniFrac result, not a technicality.
TAB=$PROJ/microbiome/$AMP/seqtab_biom.qza
if [[ -s "$TAB" ]]; then
  qiime fragment-insertion filter-features \
    --i-table "$TAB" --i-tree insertion-tree.qza \
    --o-filtered-table placed-table.qza --o-removed-table unplaced-table.qza
fi

qiime tools export --input-path insertion-tree.qza --output-path exported-insertion
mv exported-insertion/tree.nwk asv_16S_sepp.nwk
PLACED=$(grep -o '[A-Za-z0-9]\{32\}' asv_16S_sepp.nwk | sort -u | wc -l)
TOTAL=$(grep -c '^>' "$FA")
# SEPP has no reject option: every fragment is placed somewhere. A count of 100% placed is
# therefore NOT evidence that placements are well supported, and must not be reported as such.
# It only rules out the opposite failure, features silently missing from the tree.
echo "PLACED: $PLACED / $TOTAL ASVs present as tips (placement quality is a separate question)"

# The insertion output carries the entire reference backbone, so the tree has ~200k tips while
# the feature table has ~5k. Prune to our own features: the topology and branch lengths still
# come from the reference alignment, but the tree that leaves this script has exactly the tips
# the table has. Everything downstream uses the pruned tree.
Rscript -e '
  suppressPackageStartupMessages(library(ape))
  tr <- read.tree("asv_16S_sepp.nwk")
  ids <- sub("^>", "", grep("^>", readLines("'"$FA"'"), value = TRUE))
  tr <- keep.tip(tr, intersect(tr$tip.label, ids))
  write.tree(tr, "asv_16S_sepp_pruned.nwk")
  cat(sprintf("pruned tree tips: %d\n", ape::Ntip(tr)))'

# ---- sensitivity: de novo tree -----------------------------------------------------------
# Built only to demonstrate robustness, never as the primary. MAFFT then FastTree under GTR,
# with the alignment masked for highly gapped / conserved columns as QIIME's standard pipeline
# does, since an unmasked V4 alignment is dominated by uninformative sites.
if [[ ! -s asv_16S_denovo.nwk ]]; then
  echo "=== de novo tree (sensitivity only) ==="; date
  qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences rep-seqs.qza \
    --p-n-threads "$SLURM_CPUS_PER_TASK" \
    --o-alignment aln.qza --o-masked-alignment aln-masked.qza \
    --o-tree denovo-tree.qza --o-rooted-tree denovo-tree-rooted.qza
  qiime tools export --input-path denovo-tree-rooted.qza --output-path exported-denovo
  mv exported-denovo/tree.nwk asv_16S_denovo.nwk
  date
fi

echo "DONE"
ls -la "$OUT"/*.nwk
