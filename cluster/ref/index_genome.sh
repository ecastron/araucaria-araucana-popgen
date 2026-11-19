#!/bin/bash
#SBATCH --job-name=mm2idx
#SBATCH --partition=basic
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G          # >257G pins to the 2 TB fat node n001
#SBATCH --time=24:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/mm2idx_%j.out
set -euo pipefail
module load Minimap2/minimap2
D=/home/ecastron/araucaria/refs/araucana
cd "$D"
echo "node: $(hostname)  start: $(date)"
# -I 25G keeps the whole 20 Gb genome in ONE index part, so every read is compared
# against the entire genome in a single pass. Multi-part indexes split the reference,
# which would let a read map "uniquely" within a part while its true paralog sits in
# another part -- exactly the artifact this experiment is trying to detect.
minimap2 -t 16 -x sr -I 25G -d v1.sr.mmi v1.fasta.gz
echo "index done: $(date)"; ls -lh v1.sr.mmi
