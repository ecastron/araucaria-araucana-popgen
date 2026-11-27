#!/bin/bash
#SBATCH --job-name=ancnull
#SBATCH --partition=gpu,basic
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --array=1-10
#SBATCH --output=/home/ecastron/araucaria/logs/ancnull_%A_%a.out
# Revision (26 Nov 2026): sNMF null under zero gene flow (ancestry_null.py, .R)
set -euo pipefail
cd /home/ecastron/araucaria
PY=/home/ecastron/envs/feems/bin/python
mkdir -p results/revision/ancnull
$PY - <<'PYEOF'
import pickle, numpy as np
D = pickle.load(open("results/Q1/abc/observed2.pkl", "rb"))
n = np.bincount(D["deme"], minlength=len(D["stands"]))
import os
tmp = f"results/revision/ancnull/stand_order.{os.getpid()}.tmp"
with open(tmp, "w") as fh:
    fh.write("stand\tn\tcoastal\n")
    for s, k, c in zip(D["stands"], n, D["coastal"]): fh.write(f"{s}\t{k}\t{'TRUE' if c else 'FALSE'}\n")
os.replace(tmp, "results/revision/ancnull/stand_order.tsv")
PYEOF
$PY -u scripts/revision/ancestry_null.py ${SLURM_ARRAY_TASK_ID} 10
module load R/R-4.4.0
Rscript scripts/revision/ancestry_null.R ${SLURM_ARRAY_TASK_ID} 10
