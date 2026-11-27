#!/bin/bash
#SBATCH -J export_rev
#SBATCH --partition=gpu,basic
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH -o /home/ecastron/araucaria/logs/export_rev_%j.out
# Revision (27 Nov 2026): genotype and ASV tables for archiving (export_genotypes.py, export_asv_table.R)
set -euo pipefail
cd /home/ecastron/araucaria
/home/ecastron/envs/feems/bin/python -u scripts/revision/export_genotypes.py
module load R/R-4.4.0
Rscript scripts/revision/export_asv_table.R
ls -la results/revision/araucaria_163trees_78026snps.vcf.gz results/revision/its1_asv_*.csv.gz
