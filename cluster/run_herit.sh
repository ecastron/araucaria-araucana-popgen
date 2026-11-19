#!/bin/bash
#SBATCH -J herit
#SBATCH -p gpu
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 04:00:00
#SBATCH -o /home/ecastron/araucaria/logs/herit_%j.out
#SBATCH -e /home/ecastron/araucaria/logs/herit_%j.err
module load R
cd /home/ecastron/araucaria
Rscript scripts/q6_heritability.R
