#!/bin/bash
#SBATCH -J q6life
#SBATCH -p gpu
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o /home/ecastron/araucaria/results/Qlife/q6life_%j.out
#SBATCH -e /home/ecastron/araucaria/results/Qlife/q6life_%j.err
module load R/R-4.4.0
cd /home/ecastron/araucaria
Rscript scripts/q6_lifestyle.R
