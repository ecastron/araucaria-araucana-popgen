#!/bin/bash
#SBATCH -J mglm
#SBATCH -p gpu
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 06:00:00
#SBATCH -o /home/ecastron/araucaria/logs/mglm_%j.out
#SBATCH -e /home/ecastron/araucaria/logs/mglm_%j.err
module load R
# GSL is installed on the login node only; mvabund.so needs libgsl at runtime
export LD_LIBRARY_PATH=/home/ecastron/araucaria/lib:$LD_LIBRARY_PATH
cd /home/ecastron/araucaria
Rscript scripts/manyglm_q5q6.R
