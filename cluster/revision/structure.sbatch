#!/bin/bash
#SBATCH -J struct_rev
#SBATCH --partition=basic
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH -o /home/ecastron/araucaria/logs/struct__rev%j.out
# Revision (25 Nov 2026). Runs structure.py, then PCA and ADMIXTURE
# (plink 1.9, ADMIXTURE 1.3.0 from the Popgen module, as q1_fixed.sh) on the 163 analyzed trees:
#   ld   : LD-pruned as in q1_fixed.sh (--indep-pairwise 50 5 0.2), now on 163 trees, not 177
#   thin1: one SNP per RAD locus
# and sNMF (snmf.R) on all SNPs and on the thinned set.
# Outputs in results/revision/: pca_{ld,thin1}.eigen{vec,val}, admix_{ld,thin1}.2.Q, *_cv.txt, snmf_*
set -euo pipefail
P=/home/ecastron/araucaria; S=$P/scripts/revision; OUT=$P/results/revision
V=$P/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf
mkdir -p $OUT; cd $OUT
/home/ecastron/envs/feems/bin/python -u $S/structure.py

module load Popgen/popgen-1.0
COMMON="--vcf $V --allow-extra-chr --double-id --set-missing-var-ids @:# --keep $OUT/keep163.plink --biallelic-only strict --snps-only just-acgt"
plink $COMMON --indep-pairwise 50 5 0.2 --out ld_prune
plink $COMMON --extract ld_prune.prune.in --pca 10 --out pca_ld
plink $COMMON --extract thin1.extract --pca 10 --out pca_thin1
echo "LD-pruned SNPs: $(wc -l < ld_prune.prune.in)   thinned SNPs: $(wc -l < thin1.extract)"
for t in ld thin1; do
  ex=$([ $t = ld ] && echo ld_prune.prune.in || echo thin1.extract)
  plink $COMMON --extract $ex --make-bed --out admix_$t
  # ADMIXTURE needs integer chromosome codes; RAD loci are unplaced, so code all as 1
  awk '{$1=1; print}' admix_$t.bim > tmp.bim && mv tmp.bim admix_$t.bim
  for K in 1 2 3; do admixture --cv=10 -j8 -s 2026 admix_$t.bed $K | grep -i "CV error" >> admix_${t}_cv.txt; done
done

module load R/R-4.4.0
Rscript $S/snmf.R
echo DONE
