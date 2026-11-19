#!/bin/bash
#SBATCH --job-name=arau-Q1b
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/Q1b_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/Q1b_%j.err

# Q1, corrected run. Three problems in the first attempt:
#
#  1. FST came out with a NEGATIVE MEAN (-0.050), which is impossible given independently
#     demonstrated structure (isolation by distance r = 0.601; stand R2 significant in Q5). The
#     hand-rolled Weir & Cockerham dropped the unequal-sample-size correction and the
#     heterozygosity term. Replaced with vcftools --weir-fst-pop, a validated implementation.
#  2. PCA ran on only 767 loci because it required >=50% of libraries at >=20x. Use plink on the
#     full cr30 set instead -- more loci, and LD-pruned, which PCA needs and the first run lacked.
#  3. AMOVA failed outright: it demanded libraries with ZERO undefined pairs and only 1 of 177
#     qualified. 170 have under 10% missing, so drop the few bad ones and impute the rest.
#
# Per-stand diversity from the first run is kept -- that calculation was verified independently
# against the heterozygote-deficit diagnostic and is not affected by any of the above.

set -uo pipefail
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
P=/home/ecastron/araucaria
V=$P/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf
OUT=$P/results/Q1
mkdir -p "$OUT"

# Tools come from the Popgen module, not a conda activate: module load propagates correctly to
# non-interactive Slurm shells. ADMIXTURE, plink and vcftools are all REQUIRED -- fail loudly
# rather than skip, so a missing tool cannot yield a partial Q1 that looks complete in ten days.
module load Popgen/popgen-1.0
for t in vcftools plink admixture; do
  command -v $t >/dev/null || { echo "FATAL: $t unavailable. Q1 cannot complete." >&2; exit 1; }
done
if true; then
  echo "=== FST via vcftools (validated Weir & Cockerham) ==="
  mkdir -p "$OUT/fst"
  # ipyrad names samples from the symlink TARGET, so the VCF carries a .1 mate suffix
  # (SRR11841624.1) that the popmap does not. Build the population files from the VCF's own
  # sample names, or vcftools silently keeps 0 of 177 individuals.
  grep -m1 '^#CHROM' "$V" | cut -f10- | tr '\t' '\n' > "$OUT/fst/vcf_samples.txt"
  awk -F'\t' 'NR==FNR{stand[$1]=$2; next}
       {bare=$1; sub(/\.[12]$/,"",bare); if (bare in stand) print $1 > "'$OUT'/fst/pop_" stand[bare] ".txt"}' \
      $P/metadata/popmap_177.tsv "$OUT/fst/vcf_samples.txt"
  echo "population files: $(ls "$OUT"/fst/pop_*.txt | wc -l) stands, \
$(cat "$OUT"/fst/pop_*.txt | wc -l) libraries matched to VCF names"
  ls "$OUT"/fst/pop_*.txt | sed 's/.*pop_//;s/\.txt//' > "$OUT/fst/stands.txt"
  while read -r A; do while read -r B; do
    [ "$A" \< "$B" ] || continue
    r=$(vcftools --vcf "$V" --weir-fst-pop "$OUT/fst/pop_$A.txt" \
          --weir-fst-pop "$OUT/fst/pop_$B.txt" --stdout 2>&1 |
        grep -i "weighted Fst estimate" | awk '{print $NF}')
    echo -e "$A\t$B\t$r"
  done < "$OUT/fst/stands.txt"; done < "$OUT/fst/stands.txt" | tee "$OUT/q1_fst_vcftools.tsv"

  echo; echo "=== PCA via plink (LD-pruned) ==="
  plink --vcf "$V" --allow-extra-chr --double-id --set-missing-var-ids @:# \
        --indep-pairwise 50 5 0.2 --out "$OUT/prune" >/dev/null 2>&1
  plink --vcf "$V" --allow-extra-chr --double-id --set-missing-var-ids @:# \
        --extract "$OUT/prune.prune.in" --pca 10 --out "$OUT/q1_pca_plink" >/dev/null 2>&1
  if [ -s "$OUT/q1_pca_plink.eigenval" ]; then
    echo "eigenvalues (% of variance):"
    awk '{s+=$1; v[NR]=$1} END{for(i=1;i<=10;i++) printf "  PC%-2d %6.2f%%\n", i, 100*v[i]/s}' \
      "$OUT/q1_pca_plink.eigenval"
    echo "LD-pruned SNPs used: $(wc -l < "$OUT/prune.prune.in")"
  fi

  if true; then
    echo; echo "=== ADMIXTURE K=1..8 (cross-validation) ==="
    cd "$OUT"
    plink --vcf "$V" --allow-extra-chr --double-id --set-missing-var-ids @:# \
          --extract "$OUT/prune.prune.in" --make-bed --out admix >/dev/null 2>&1
    cut -f1 -d' ' admix.bim | sort -u | head -1 >/dev/null
    awk '{$1=1; print}' admix.bim > tmp.bim && mv tmp.bim admix.bim   # ADMIXTURE needs numeric chr
    for K in 1 2 3 4 5 6 7 8; do
      admixture --cv admix.bed $K -j8 2>&1 | grep -E "^CV" || echo "  K=$K failed"
    done | tee "$OUT/q1_admixture_cv.txt"
  fi
fi

echo; echo "=== AMOVA (NA-tolerant) ==="
conda activate dada2
Rscript $P/scripts/q1_amova.R
echo "Q1b COMPLETE"
