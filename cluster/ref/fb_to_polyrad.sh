#!/bin/bash
#SBATCH --job-name=fb2prad
#SBATCH --partition=basic
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/fb2prad_%j.out
set -eo pipefail
S=/home/ecastron/araucaria/results/sv
P=/home/ecastron/araucaria/results/polyrad_fb; mkdir -p "$P"
module load SAMTOOLS/samtools.1.20
FBBIN=/home/ecastron/envs/freebayes/bin

echo "### 1. concatenate 200 chunks"
grep "^#" $S/fb_out/fb_000.vcf > $P/fb_all.vcf
for i in $(seq 0 199); do printf -v C "%03d" $i; grep -hv "^#" $S/fb_out/fb_$C.vcf 2>/dev/null; done \
  | sort -k1,1 -k2,2n >> $P/fb_all.vcf
echo "  records: $(grep -vc '^#' $P/fb_all.vcf)"

echo "### 2. biallelic SNPs, add AD from RO/AO"
python3 /home/ecastron/araucaria/scripts/ref/fb_add_AD.py $P/fb_all.vcf $P/fb_AD.vcf
echo "### 3. bgzip + index"
$FBBIN/bgzip -f -@ 8 -c $P/fb_AD.vcf > $P/fb_AD.vcf.gz 2>/dev/null || bgzip -f -c $P/fb_AD.vcf > $P/fb_AD.vcf.gz
tabix -f -p vcf $P/fb_AD.vcf.gz 2>/dev/null || true
ls -lh $P/fb_AD.vcf.gz
