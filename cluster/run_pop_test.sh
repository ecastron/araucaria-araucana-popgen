#!/bin/bash
#SBATCH --job-name=arau-poptest
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/poptest_%j.out
#SBATCH --error=/home/ecastron/araucaria/logs/poptest_%j.err
set -euo pipefail
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate ipyrad
V=$(ls /home/ecastron/araucaria/rad/ipyrad_full/*_outfiles/*.vcf 2>/dev/null | head -1)
[[ -s "$V" ]] || { echo "no full-assembly VCF found" >&2; exit 1; }
echo "using $V"
python /home/ecastron/araucaria/scripts/pop_structure_full.py "$V"
