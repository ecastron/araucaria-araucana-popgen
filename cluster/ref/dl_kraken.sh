#!/bin/bash
#SBATCH --job-name=dlk2
#SBATCH --partition=basic
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/dlk2_%j.out
set -euo pipefail
S=https://genome-idx.s3.amazonaws.com/kraken
D=/home/ecastron/araucaria/refs/kraken2
cd "$D"
# PlusPFP: archaea bacteria viral plasmid human UniVec_Core protozoa FUNGI PLANT.
# Plant is the reason for this collection over PlusPF -- without Viridiplantae our
# own Araucaria reads return "unclassified", which is indistinguishable from a
# genuinely unknown read and would make the whole test uninterpretable.
echo ">>> capped PlusPFP-16 first (fast, gives a preliminary read)"
mkdir -p pluspfp16 && cd pluspfp16
curl -sS -L -C - --retry 20 --retry-delay 15 -O "$S/k2_pluspfp_16_GB_20260626.tar.gz"
tar xzf k2_pluspfp_16_GB_20260626.tar.gz && rm -f k2_pluspfp_16_GB_20260626.tar.gz
echo "pluspfp16 ready: $(ls -lh hash.k2d 2>/dev/null | awk '{print $5}')"
cd "$D"
echo ">>> full PlusPFP (171.8 GB archive, 231.5 GB index)"
mkdir -p pluspfp && cd pluspfp
curl -sS -L -C - --retry 30 --retry-delay 20 -O "$S/k2_pluspfp_20260626.tar.gz"
echo "archive down: $(ls -lh k2_pluspfp_20260626.tar.gz | awk '{print $5}')"
tar xzf k2_pluspfp_20260626.tar.gz && rm -f k2_pluspfp_20260626.tar.gz
echo "pluspfp ready:"; ls -lh
