#!/bin/bash
#SBATCH -J clockda
#SBATCH --partition=basic
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH -o /home/ecastron/araucaria/logs/clockda_%j.out
/home/ecastron/envs/feems/bin/python -u /home/ecastron/araucaria/scripts/clock_da.py
