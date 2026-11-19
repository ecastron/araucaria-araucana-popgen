#!/bin/bash
#SBATCH -J f1cv
#SBATCH -p gpu
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o /home/ecastron/araucaria/results/F1C/f1c_%j.out
#SBATCH -e /home/ecastron/araucaria/results/F1C/f1c_%j.err
module load R/R-4.4.0
cd /home/ecastron/araucaria
Rscript scripts/f1c_varpart.R
