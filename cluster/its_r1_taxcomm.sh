#!/bin/bash
#SBATCH --job-name=arau-r1tax
#SBATCH --partition=basic
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G   # n002/n003 cap at 257G; anything above pins the job to n001
#SBATCH --time=12:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/r1tax_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/r1tax_%j.err

# The interactive R1-vs-merged comparison died at assignTaxonomy with "Memory allocation
# failed" -- UNITE against two ASV sets at once exceeded the login node. Retention and ASV
# counts were already obtained; this rerun completes the two outstanding pieces, taxonomy
# resolution and community concordance, with memory to spare.

set -euo pipefail
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate dada2
Rscript /home/ecastron/araucaria/scripts/its_r1_taxcomm.R
