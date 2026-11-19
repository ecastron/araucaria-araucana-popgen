#!/bin/bash
#SBATCH --job-name=arau-sra
#SBATCH --partition=basic
#SBATCH --array=1-774%10
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=/home/ecastron/araucaria/logs/sra_%A_%a.out
#SBATCH --error=/home/ecastron/araucaria/logs/sra_%A_%a.err

# Download one SRA run per array task from NCBI, convert to gzipped FASTQ, validate.
#
# NCBI only. Chile routes to NCBI far better than to EBI: a single prefetch stream
# sustains ~5 MB/s here, against ~0.65 MB/s per aria2 connection from ENA.
#
# .sra accession and fasterq-dump scratch both live on node-local NVMe (/scratch),
# which is per-node and invisible from the login node. Only the final gzipped FASTQ
# is written to the NFS home. The .sra is deleted after conversion.
#
# Idempotent: a run whose gzipped mates already exist is skipped, so the array can be
# resubmitted over a partially complete download.

set -uo pipefail

PROJ=/home/ecastron/araucaria
TASKS=$PROJ/metadata/sra_tasks.tsv

module load SRAtoolkit/sratoolkit-3.1.1

read -r idx dataset run < <(awk -v i="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==i' "$TASKS" | tr -d '\r')

if [[ -z "${run:-}" ]]; then
  echo "FATAL: no task row for array index $SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi

DEST=$PROJ/raw/$dataset
mkdir -p "$DEST"

if [[ -s "$DEST/${run}_1.fastq.gz" && -s "$DEST/${run}_2.fastq.gz" ]]; then
  echo "SKIP  $run (both mates already present)"
  exit 0
fi

WORK=/scratch/$USER/arau_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$WORK"
# Always clear node-local scratch, including on failure, so a node does not fill up.
trap 'rm -rf "$WORK"' EXIT

# --max-size u lifts the default 20G cap; the largest RAD run is well under it but the
# cap would otherwise silently skip a big run.
if ! prefetch --max-size u -O "$WORK" "$run" >"$WORK/prefetch.log" 2>&1; then
  echo "FAIL  $run prefetch failed" >&2
  tail -5 "$WORK/prefetch.log" >&2
  exit 1
fi

# prefetch verifies the download against NCBI's own checksums; vdb-validate re-checks
# the archive's internal consistency before we spend CPU converting it.
if ! vdb-validate "$WORK/$run" >"$WORK/validate.log" 2>&1; then
  echo "FAIL  $run vdb-validate failed" >&2
  tail -5 "$WORK/validate.log" >&2
  exit 1
fi

if ! fasterq-dump --split-files --threads "$SLURM_CPUS_PER_TASK" \
       --temp "$WORK/tmp" --outdir "$WORK" "$run" >"$WORK/dump.log" 2>&1; then
  echo "FAIL  $run fasterq-dump failed" >&2
  tail -5 "$WORK/dump.log" >&2
  exit 1
fi

if [[ ! -s "$WORK/${run}_1.fastq" || ! -s "$WORK/${run}_2.fastq" ]]; then
  echo "FAIL  $run expected paired output missing after fasterq-dump" >&2
  ls -la "$WORK" >&2
  exit 1
fi

pigz -p "$SLURM_CPUS_PER_TASK" "$WORK/${run}_1.fastq" "$WORK/${run}_2.fastq"

# Move into place only after both mates are fully written, so an interrupted task never
# leaves a half-file that the skip check above would mistake for complete.
mv "$WORK/${run}_1.fastq.gz" "$WORK/${run}_2.fastq.gz" "$DEST/"

n1=$(stat -c%s "$DEST/${run}_1.fastq.gz")
n2=$(stat -c%s "$DEST/${run}_2.fastq.gz")
echo "OK    $run $dataset ${n1} ${n2} bytes"
