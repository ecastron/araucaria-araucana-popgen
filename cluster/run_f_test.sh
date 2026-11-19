#!/bin/bash
#SBATCH --job-name=arau-Ftest
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/Ftest_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/Ftest_%j.err
set -uo pipefail
module load Popgen/popgen-1.0
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
P=/home/ecastron/araucaria; V=$P/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf
OUT=$P/results/Q2; cd $OUT
for DP in 10 20; do
  conda activate ipyrad
  python $P/scripts/gl_highdepth.py "$V" "$OUT/gl_dp${DP}.beagle" $DP || continue
  gzip -f "$OUT/gl_dp${DP}.beagle"
  NS=$(zcat "$OUT/gl_dp${DP}.beagle.gz" | head -1 | awk '{print (NF-3)/3}')
  NL=$(( $(zcat "$OUT/gl_dp${DP}.beagle.gz" | wc -l) - 1 ))
  ngsRelate -G "$OUT/gl_dp${DP}.beagle.gz" -n "$NS" -L "$NL" -z "$OUT/samples.txt" \
    -O "$OUT/ngsrelate_dp${DP}.tsv" -p 16 >/dev/null 2>&1
  conda activate dada2
  Rscript -e "
    r <- read.delim('$OUT/ngsrelate_dp${DP}.tsv', check.names=FALSE)
    cat(sprintf('depth floor ${DP}x: median F = %.3f | median theta = %.4f | median nSites = %d | pairs = %d\n',
        median(c(r\$Fa, r\$Fb), na.rm=TRUE), median(r\$theta, na.rm=TRUE),
        median(r\$nSites), nrow(r)))"
done
echo "baseline (no depth floor): median F = 0.700"
