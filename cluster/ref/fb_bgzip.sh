#!/bin/bash
#SBATCH --job-name=fbgz
#SBATCH --partition=basic
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/fbgz_%j.out
set -eo pipefail
P=/home/ecastron/araucaria/results/polyrad_fb
/opt/samtools/bgzip -f -@ 8 -c $P/fb_AD.vcf > $P/fb_AD.vcf.gz
/opt/samtools/tabix -f -p vcf $P/fb_AD.vcf.gz
echo "verify from the .gz:"
/opt/samtools/bgzip -dc $P/fb_AD.vcf.gz | grep -m1 -v "^#" | cut -f9,10,11
/opt/samtools/bgzip -dc $P/fb_AD.vcf.gz | grep -vc "^#" | xargs echo "  SNPs:"
ls -lh $P/fb_AD.vcf.gz $P/fb_AD.vcf.gz.tbi
