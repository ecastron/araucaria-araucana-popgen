#!/bin/bash
#SBATCH --job-name=abcfinal
#SBATCH --partition=basic
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=4:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/abcfinal_%j.out
set -euo pipefail
~/envs/feems/bin/python /home/ecastron/araucaria/scripts/abc/final_abc.py
