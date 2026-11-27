#!/bin/bash
#SBATCH -J its_rev
#SBATCH --partition=basic
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH -o /home/ecastron/araucaria/logs/its__rev%j.out
# Revision (25 Nov 2026): runs its.R (stand-level range test, rarefaction sensitivity)
module load R/R-4.4.0
Rscript /home/ecastron/araucaria/scripts/revision/its.R
