#!/bin/bash
#SBATCH --job-name=bt2build
#SBATCH --partition=basic
#SBATCH --cpus-per-task=32
#SBATCH --mem=1200G          # >257G pins to n001 (2 TB, 40 cores)
#SBATCH --time=96:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/bt2build_%j.out
# the Bowtie2 module sources a conda activate.d script that dereferences JAVA_HOME;
# under `set -u` that unbound variable aborts the job before bowtie2-build runs.
export JAVA_HOME="${JAVA_HOME:-}"
set -eo pipefail
module load Bowtie2/bowtie2
D=/home/ecastron/araucaria/refs/araucana
cd "$D"
echo "node $(hostname)  start $(date)"
bowtie2-build --version | head -2
# 19.3 Gb exceeds the 4 Gb small-index limit, so --large-index is required
bowtie2-build --threads 32 --large-index v1.fasta v1_bt2
echo "done $(date)"; ls -lh v1_bt2*
