#!/bin/bash
#SBATCH --job-name=fbconv
#SBATCH --partition=basic
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/fbconv_%j.out
set -eo pipefail
P=/home/ecastron/araucaria/results/polyrad_fb
python3 /home/ecastron/araucaria/scripts/ref/fb_add_AD.py $P/fb_all.vcf $P/fb_AD.vcf
/opt/samtools/bgzip -f -c $P/fb_AD.vcf > $P/fb_AD.vcf.gz
/opt/samtools/tabix -f -p vcf $P/fb_AD.vcf.gz 2>/dev/null || true
echo "--- spot check ---"
grep -m1 -v "^#" $P/fb_AD.vcf | cut -f9,10,11
grep -v "^#" $P/fb_AD.vcf | head -5000 | awk -F'\t' '{for(i=10;i<=NF;i++) if($i!="."){n++; if($i ~ /:[0-9]+,[0-9]+$/) ok++}} END{printf "cells with data: %d | ending in AD: %d (%.1f%%)\n",n,ok,100*ok/n}'
ls -lh $P/fb_AD.vcf.gz
