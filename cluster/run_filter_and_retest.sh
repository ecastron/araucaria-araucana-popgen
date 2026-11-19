#!/bin/bash
#SBATCH --job-name=arau-filt
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --time=12:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/filt_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/filt_%j.err
set -euo pipefail
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate ipyrad
P=/home/ecastron/araucaria
V=$P/rad/ipyrad_full/full_ct091_outfiles/full_ct091.vcf

python $P/scripts/filter_vcf_callrate.py "$V" 0.3 0.5 0.7

# Re-run the structure test at every threshold. Reporting all three matters: a positive that
# appears at only one call-rate floor out of three is threshold-shopping, not a result.
for c in 30 50 70; do
  F="${V%.vcf}.cr${c}.vcf"
  [[ -s "$F" ]] || { echo "missing $F"; continue; }
  echo; echo "############ call-rate floor ${c}% ############"
  python $P/scripts/pop_structure_full.py "$F" 2>&1 | tail -22
done
