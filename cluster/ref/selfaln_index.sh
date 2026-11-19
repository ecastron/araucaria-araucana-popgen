#!/bin/bash
#SBATCH --job-name=sa_idx
#SBATCH --partition=basic
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --time=6:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/sa_idx_%j.out
set -euo pipefail
module load Minimap2/minimap2
D=/home/ecastron/araucaria/refs/araucana
# asm20 (not asm5): segmental duplications in a repeat-rich conifer can be
# appreciably diverged, and we would rather over-collect and filter on identity.
minimap2 -x asm20 -t 16 -I 25G -d "$D/v1.asm20.mmi" "$D/v1.fasta"
ls -lh "$D/v1.asm20.mmi"
