#!/bin/bash
#SBATCH -J q6dp20
#SBATCH -p gpu
#SBATCH --nodelist=g001
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 8:00:00
#SBATCH -o /home/ecastron/araucaria/logs/q6dp20_%j.out
#SBATCH -e /home/ecastron/araucaria/logs/q6dp20_%j.err
set -uo pipefail
P=/home/ecastron/araucaria
module load R 2>/dev/null
export LD_LIBRARY_PATH=$P/lib:${LD_LIBRARY_PATH:-}
for s in q6_rerun_clean.R q6_ksens_mrm.R q6_power.R; do
  echo "######################## $s ########################"
  Rscript $P/scripts/dp20/$s 2>&1 || echo "!! $s FAILED"
done
echo "######################## DONE $(date) ########################"
