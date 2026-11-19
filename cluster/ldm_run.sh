#!/bin/bash
#SBATCH -J ldm
#SBATCH -p gpu
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 04:00:00
#SBATCH -o /home/ecastron/araucaria/results/fungfig/ldm_%j.out
#SBATCH -e /home/ecastron/araucaria/results/fungfig/ldm_%j.err
module load R/R-4.4.0
cd /home/ecastron/araucaria
# ITS_TABLE and OUT_DIR pass through from the submitting shell (sbatch exports the
# environment); a table name given as an argument is forwarded. No argument = merged table.
Rscript scripts/ldm_run.R "$@"
