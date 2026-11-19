#!/bin/bash
#SBATCH -J q6pow
#SBATCH -p gpu
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 06:00:00
#SBATCH -o /home/ecastron/araucaria/logs/pow_%j.out
#SBATCH -e /home/ecastron/araucaria/logs/pow_%j.err
module load R
export LD_LIBRARY_PATH=/home/ecastron/araucaria/lib:$LD_LIBRARY_PATH
cd /home/ecastron/araucaria
Rscript scripts/q6_power.R
