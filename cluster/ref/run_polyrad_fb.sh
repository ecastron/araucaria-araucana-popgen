#!/bin/bash
#SBATCH --job-name=pradfb
#SBATCH --partition=basic
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --time=72:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/pradfb_%j.out
set -eo pipefail
module load R/R-4.4.0
export R_LIBS_USER=/home/ecastron/R/libs
export LD_LIBRARY_PATH=/usr/lib64:/lib64:${LD_LIBRARY_PATH:-}
Rscript --vanilla /home/ecastron/araucaria/scripts/ref/polyrad_fb.R
