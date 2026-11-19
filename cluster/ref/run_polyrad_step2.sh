#!/bin/bash
#SBATCH --job-name=polyrad2
#SBATCH --partition=basic
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/polyrad2_%j.out
set -uo pipefail
module load R/R-4.4.0
export R_LIBS_USER=/home/ecastron/R/libs
export LD_LIBRARY_PATH=/usr/lib64:/lib64:${LD_LIBRARY_PATH:-}
Rscript --vanilla /home/ecastron/araucaria/scripts/ref/run_polyrad_step2.R
