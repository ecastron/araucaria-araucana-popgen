#!/bin/bash
#SBATCH --job-name=arau-Q1
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/Q1_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/Q1_%j.err
set -uo pipefail
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate ipyrad
P=/home/ecastron/araucaria
python $P/scripts/q1_popgen.py || echo "Q1 python stage FAILED"

# AMOVA in R on the distance matrix Q1 wrote. Separate stage so a failure here cannot lose the
# diversity, FST and PCA results already on disk.
conda activate dada2
Rscript $P/scripts/q1_amova.R || echo "Q1 AMOVA stage FAILED"

# ADMIXTURE, only if the popgen env finished installing. Optional by design: PCA already gives
# the structure picture, and ancestry coefficients are a complement rather than a dependency.
if [ -x /opt/miniconda/miniconda-py312_24.3.0-0/envs/popgen/bin/admixture ]; then
  bash $P/scripts/q1_admixture.sh || echo "Q1 ADMIXTURE stage FAILED"
else
  echo "ADMIXTURE not installed; skipped (PCA covers structure)"
fi
echo "Q1 COMPLETE"
