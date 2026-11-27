#!/bin/bash
# Runs the diagnostic scripts whose stdout supplies numbers to figdata/het_diagnosis.json,
# saving each stdout under results/derived/het_diag/ so het_diagnosis_assemble.py can read
# every value from a file.
#SBATCH -J hetdiag
#SBATCH -p basic
#SBATCH -c 8
#SBATCH --mem 24G
#SBATCH -t 04:00:00
#SBATCH -o /home/ecastron/araucaria/logs/hetdiag_%j.out
#SBATCH -e /home/ecastron/araucaria/logs/hetdiag_%j.err
set -euo pipefail
source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
conda activate /home/ecastron/envs/feems
P=/home/ecastron/araucaria
V=$P/rad/ipyrad_full/full_ct091_outfiles/full_ct091.cr30.vcf
POP=$P/metadata/popmap_177.tsv
O=$P/results/derived/het_diag
mkdir -p $O
# independent single-threaded scans of the same VCF: run them side by side
python $P/scripts/ref/two_phenomena.py      $V $POP    > $O/two_phenomena.txt      &
python $P/scripts/ref/rarefy_and_cv.py      $V $POP 30 > $O/rarefy_and_cv.txt      &
python $P/scripts/ref/locus_depth_filter.py $V $POP    > $O/locus_depth_filter.txt &
python $P/scripts/ref/cnv_filter.py         $V $POP    > $O/cnv_filter.txt         &
python $P/scripts/ref/hwe_per_stand.py      $V $POP    > $O/hwe_per_stand.txt      &
python $P/scripts/ref/hwe_by_maf.py         $V         > $O/hwe_by_maf.txt         &
wait
echo done
