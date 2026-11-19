#!/bin/bash
#SBATCH -J mibd
#SBATCH -p gpu
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 03:00:00
#SBATCH -o /home/ecastron/araucaria/logs/mibd_%j.out
#SBATCH -e /home/ecastron/araucaria/logs/mibd_%j.err
module load R
cd /home/ecastron/araucaria
Rscript scripts/q9_microbiome_ibd.R
