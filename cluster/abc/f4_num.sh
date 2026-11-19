#!/bin/bash
#SBATCH --job-name=f4num
#SBATCH --partition=basic
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=4:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/f4num_%j.out
set -euo pipefail
~/envs/feems/bin/python /home/ecastron/araucaria/scripts/abc/f4_numbers.py
